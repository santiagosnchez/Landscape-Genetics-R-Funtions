dist_spline_res <- function(r, sp){
    require(raster)
    res <- r
    dat <- extract(r, sp)
    for (i in 1:dim(dat)[2]){
        d <- hist(dat[,i], plot=F, breaks=50)
        x <- d$breaks[-1]
        y <- d$density
        fitsp <- smooth.spline(x, y, df=5)
        xx <- r[[i]][ !is.na(r[[i]][]) ]
        yy <- predict(fitsp, xx)$y
        res[[i]][ !is.na(res[[i]][]) ] <- 1-yy
    }
    return(res)
}

rescale_raster <- function(r){
    require(raster)
    require(scales)
    res <- r
    res[ !is.na(res[]) ] <- rescale(res[ !is.na(res[]) ], c(1,10))
    return(res)
}

get_cos_mat <- function(r, sp){
    require(raster)
    require(gdistance)
    n <- dim(sp@coords)[1]
    res <- matrix(0, ncol=dim(r)[3], nrow=n*(n-1) / 2)
    for (i in 1:dim(r)[3]){
        tr <- transition(r[[i]], function(x) 1/mean(x), directions=8)
        cst <- costDistance(tr, sp)
        res[,i] <- as.vector(cst)
    }
    return(res)
}

# needs to have circuitscape installed

runCS <- function(layer, sites){
    require(raster)
    cls <- cellFromXY(layer, sites)
    if (length(unique(cls)) != length(cls)){
        stop("sites are do not fall within unique cells")
    }
    sites <- rasterize(x = sites,y = layer)
    dir.create("CS",showWarnings=FALSE)
    writeRaster(sites,"CS/sites_rast.asc",overwrite=TRUE)
    writeRaster(layer,"CS/resis_rast.asc",overwrite=TRUE)
    CS_ini <- c("[circuitscape options]",            
            "data_type = raster",
            "scenario = pairwise",
            paste(c("point_file =",
                    "habitat_file =",
                    "output_file ="),
            paste(c("CS/sites_rast.asc",
                            "CS/resis_rast.asc",
                            "CS/CS.out"))))
    cat(CS_ini, sep="\n", file="CS/myini.ini")
    CS_run <- paste("csrun.py", paste("CS/myini.ini"), paste("&> /dev/null"))
    system(CS_run)
    rdist <- as.dist(read.csv("CS/CS_resistances.out",sep=" ",row.names=1,header=1))
    return(rdist)
}

get_res_mat <- function(r, sp){
    require(raster)
    n <- dim(sp@coords)[1]
    res <- matrix(0, ncol=dim(r)[3], nrow=n*(n-1) / 2)
    for (i in 1:dim(r)[3]){
        rd <- runCS(r[[i]], sp)
        res[,i] <- as.vector(rd)
    }
    return(res)
}

MLPE <- function(x, responseCol=1, resistanceCol=2:dim(x)[2], pops=NULL, ID=NULL, scale=TRUE){
    # This is largely based on MLPE.lmm from ResistanceGA
    require(lme4)
    get_spr_matrix <- function(ID){
        # from ResistanceGA's ZZ.mat function
        Zl <- lapply(c("pop1", "pop2"), function(nm) Matrix::fac2sparse(ID[[nm]], 
        "d", drop = FALSE))
        ZZ <- Reduce("+", Zl[-1], Zl[[1]])
        return(ZZ)
    }
    if (all(sapply(list(pops,ID), is.null))){
        stop("At least provide a numeric vector for population assignments (pops) for each sample, 
            ordered in the same way as in the distance matrices rows/columns")
    } else if (is.null(ID) & !is.null(pops)){
        n <- length(pops)
        pw1 <- (n^2 - n)/2
        pw2 <- dim(x)[1]
        if (pw1 != pw2){
            stop("pops and dim(x)[1] are not the same length. Needs same number of individuals.")
        }
        ID <- matrix(NA, ncol=2, nrow=pw1)
        nr = 0
        message("Building \"ID\" matrix...")
        for (i in 1:n){
            if (i != n)
                for (j in (i+1):n){
                    nr = nr+1
                    ID[nr,] <- c(pops[i],pops[j])
                }
        }
        ID <- as.data.frame(ID)
        colnames(ID) <- c("pop1","pop2")
        message("Building \"ZZ\" sparse matrix...")
        ZZ <- get_spr_matrix(ID)
    } else if (!is.null(ID) & is.null(pops)){
        message("Building \"ZZ\" sparse matrix...")
        ZZ <- get_spr_matrix(ID)
    } else if (!is.null(ID) & !is.null(pops)){
        message("Building \"ZZ\" sparse matrix...")
        ZZ <- get_spr_matrix(ID)
    }
    res <- data.frame(matrix(NA, ncol=5, nrow=length(resistanceCol)))
    var <- colnames(x)[resistanceCol]
    for (i in resistanceCol){
        if (scale == TRUE){
            cat("Rescaling...\n")
            resistance = scale(x[,i], center=T, scale=T)
        } else {
            resistance = x[,i]
        }
        dat <- data.frame(ID, resistance, response = x[,responseCol])
        mod <- lFormula(response ~ resistance + (1 | pop1),
               data = dat,
               REML = FALSE)
        mod$reTrms$Zt <- ZZ
        message("Fitting MLPE model, column ",i,"...",sep="")
        dfun <- do.call(mkLmerDevfun, mod)
        opt <- optimizeLmer(dfun)
        MOD <- (mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
        res[(i-resistanceCol[1])+1,] <- summary(MOD)[14]$AICtab
    }
    colnames(res) <- names(summary(MOD)[14]$AICtab)
    rownames(res) <- var
    return(res)
}

RCM <- function(x, responseCol=1, resistanceCol=2:dim(x)[2]){
    require(ecodist)
    x <- as.data.frame(x)
    vars <- colnames(x)
    forms <- list()
    for (i in resistanceCol)
        for (j in resistanceCol)
            forms <- c(forms, as.formula(paste(vars[responseCol],"~",vars[i], "+", vars[j])))
    message("Running mantel tests...")
    res.man <- lapply(forms, mantel, data=x)
    mr <- sapply(res.man, function(x) x[1] )
    pv <- sapply(res.man, function(x) x[4] )
    rmat <- matrix(mr, ncol=length(vars)-1, nrow=length(vars)-1)
    pmat <- matrix(pv, ncol=length(vars)-1, nrow=length(vars)-1)
    rcm <- rmat
    dimnames(rcm) <- list(vars[-responseCol],vars[-responseCol])
    diag(rcm) <- 0
    rcm[ lower.tri(rcm, diag=F) ] <- rmat[ lower.tri(rmat, diag=F) ] - rmat[ upper.tri(rmat, diag=F) ]
    rcm[ upper.tri(rcm, diag=F) ] <- rmat[ upper.tri(rmat, diag=F) ] - rmat[ lower.tri(rmat, diag=F) ]
    return(rcm)
}

plot_RCM <- function(x, pal=NULL, width=1000, height=800, filename="plot"){
    require(plotly)
    ax <- rownames(x)
    if (is.null(pal))
        z <- plot_ly(z = x, type= "heatmap", x=ax, y=ax)
    else
        z <- plot_ly(z = x, colors=pal, type= "heatmap", x=ax, y=ax)
    # https://stackoverflow.com/a/46692547/1706987
    text <- paste0("function(el, x) {
  var gd = document.getElementById(el.id); 
  Plotly.downloadImage(gd, {format: 'svg', width: ",width,", height: ",height,", filename: '",file,"'});
  }")
    z %>% htmlwidgets::onRender(text)
}



