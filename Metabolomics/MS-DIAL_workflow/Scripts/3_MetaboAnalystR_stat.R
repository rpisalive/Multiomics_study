library(MetaboAnalystR)

setwd('path_to/working_directory')
#Prepare .csv file for MetaboAnalyst statistics analysis
raw_int <- read.csv('path_to/raw_intensities_matrix.csv', row.names = 1)
metab_SE <- readRDS('path_to/annotated_qp.RDS')
matrix <- as.data.frame(metab_SE@assays@data$counts)
raw_int <- raw_int[rownames(raw_int) %in% rownames(matrix), ]
rownames(raw_int)<- metab_SE@elementMetadata$shortname
raw_int <- rbind(metab_SE@metadata$metadata$group, raw_int)
rownames(raw_int)[1] <- "Label"
#Removing the QC columns
raw_int <- raw_int[, raw_int[1, ] != "QC"]
#Please modify the code to define labels according to your experiment, refer to https://dev.metaboanalyst.ca/docs/Format.xhtml
#The followoing is for a paired example
label <- c(1,2,3,-4,-5,-6,-7,-8,-9,4,5,6,7,8,9,-1,-2,-3)
raw_int[1,] <- label
write.csv(raw_int, "meta_stat.csv", row.names = TRUE)



#MetboAnalyst statistics
mSet<-InitDataObjects("pktable", "stat", paired = TRUE)
mSet<-Read.TextData(mSet, "meta_stat.csv", "colp", "disc")
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "MeanCenter", "S10T0", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", format ="png", dpi=72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "png", dpi=72, width=NA)
# Perform fold-change analysis on uploaded data, paired
mSet<-FC.Anal(mSet, 2.0, 0, TRUE)
# Plot fold-change analysis
mSet<-PlotFC(mSet, "fc_0_", "png", 72, width=NA)
# To view fold-change 
mSet$analSet$fc$fc.log

# Perform T-test (parametric)
mSet<-Ttests.Anal(mSet, nonpar=F, threshp=0.05, paired=TRUE, equal.var=TRUE, "fdr", TRUE)
# Plot of the T-test results
mSet<-PlotTT(mSet, imgName = "tt_0_", format = "png", dpi = 72, width=NA)

# Perform the volcano analysis
mSet<-Volcano.Anal(mSet, paired = TRUE, fcthresh = 2.0, cmpType = 0, nonpar = F, threshp = 0.1, equal.var = TRUE, pval.type = "fdr")
# Create the volcano plot
mSet<-PlotVolcano(mSet, "volcano_0_", 1, 0, format ="png", dpi=72, width=NA)

#Heatmap
mSet<-PlotStaticCorrHeatMap(mSetObj = mSet, imgName = "corr_0_", format = "png", dpi = 300,
                            width=NA, target = "col", cor.method = "pearson", colors = "bwm",
                            viewOpt = "detail", fix.col = FALSE, no.clst = FALSE, corrCutoff = 0)

#Pattern searching
# Perform correlation analysis on a pattern (a feature of interest in this case)
mSet<-FeatureCorrelation(mSet, "pearson", "Bctc_pos")

# Plot the correlation analysis on a pattern
mSet<-PlotCorr(mSet, "ptn_3_", format="png", dpi=300, width=NA)

# Perform PCA analysis
mSet<-PCA.Anal(mSet)

# Create PCA overview
mSet<-PlotPCAPairSummary(mSet, "pca_pair_0_", format = "png", dpi = 300, width=NA, 5)

# Create PCA scree plot
mSet<-PlotPCAScree(mSet, "pca_scree_0_", "png", dpi = 300, width=NA, 5)

# Create a 2D PCA score plot
mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", format = "png", dpi=300, width=NA, 1, 2, 0.95, 1, 0)


#######################################################################################################
#This part is established due to a bug in the function PlotPCA3DScoreImg, already reported to the author
.on.public.web <- FALSE
.get.mSet <- function(mSetObj=NA){
  if(.on.public.web){
    return(mSet)
  }else{
    return(mSetObj);
  }
}

GetColorSchema <- function(my.cls, grayscale=F){
  
  lvs <- levels(my.cls); 
  grp.num <- length(lvs);
  if(grayscale){
    dist.cols <- colorRampPalette(c("grey90", "grey30"))(grp.num);
    names(dist.cols) <- lvs;
  }else if(exists("colVec") && length(colVec) >0 && !any(colVec =="#NA")){
    
    # make sure to sync with cls in case user exclude some groups
    dist.cols <- colVec[names(colVec) %in% levels(my.cls)];
    
  }else{             
    if(grp.num <= 18){ # update color and respect default
      dist.cols <- pal_18[1:grp.num];
    }else{
      dist.cols <- colorRampPalette(pal_18)(grp.num);
    }
    names(dist.cols) <- lvs;
  }
  
  return (dist.cols);
}

pal_18 <- c("#e6194B", "#3cb44b", "#4363d8", "#42d4f4", "#f032e6", "#ffe119", "#911eb4", "#f58231", "#bfef45",
            "#fabebe", "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075")

ExpandSchema<-function(my.cls, schema){
  my.vec <- rep(0, length=length(my.cls));
  clsVec <- as.character(my.cls)
  grpnms <- levels(my.cls);
  for(i in 1:length(grpnms)){
    nm <- grpnms[i];
    my.vec[clsVec == nm] <- schema[nm];
  }
  
  return(my.vec);
}

Plot3D <- function(x, y, z, xlab= xlabel, ylab=ylabel, zlab=zlabel, 
                   angle =angl, color=cols, pch=pchs){
  
  if(.on.public.web){
    # make this lazy load
    if(!exists("my.plot.scatter3d")){ # public web on same user dir
      .load.scripts.on.demand("util_plot3d.Rc");    
    }
    return(my.plot.scatter3d(x, y, z, xlab=xlab, ylab=ylab, 
                             zlab=zlab, angle =angle, color=color, pch=pch));
  }else{
    return(my.plot.scatter3d(x, y, z, xlab=xlab, ylab=ylab, 
                             zlab=zlab, angle =angle, color=color, pch=pch));
  }
}

my.plot.scatter3d <- function(x, y = NULL, z = NULL, color = par("col"), pch = NULL,
                              main = NULL, sub = NULL, xlim = NULL, ylim = NULL, zlim = NULL,
                              xlab = NULL, ylab = NULL, zlab = NULL, scale.y = 1, angle = 40,
                              axis = TRUE, tick.marks = TRUE, label.tick.marks = TRUE,
                              x.ticklabs = NULL, y.ticklabs = NULL, z.ticklabs = NULL,
                              y.margin.add = 0, grid = TRUE, box = FALSE, lab = par("lab"),
                              lab.z = mean(lab[1:2]), type = "p", highlight.3d = FALSE,
                              mar = c(5, 3, 4, 3) + 0.1, col.axis = par("col.axis"),
                              col.grid = "grey", col.lab = par("col.lab"), cex.symbols = par("cex"),
                              cex.axis = 0.8 * par("cex.axis"), cex.lab = par("cex.lab"),
                              font.axis = par("font.axis"), font.lab = par("font.lab"),
                              lty.axis = par("lty"), lty.grid = 2, lty.hide = 1,
                              lty.hplot = par("lty"), log = "", ...)
{
  ## Uwe Ligges <ligges@statistik.tu-dortmund.de>,
  ## http://www.statistik.tu-dortmund.de/~ligges
  ##
  ## For MANY ideas and improvements thanks to Martin Maechler!!!
  ## Parts of the help files are stolen from the standard plotting functions in R.
  
  mem.par <- par(mar = mar)
  x.scal <- y.scal <- z.scal <- 1
  xlabel <- if (!missing(x)) deparse(substitute(x))
  ylabel <- if (!missing(y)) deparse(substitute(y))
  zlabel <- if (!missing(z)) deparse(substitute(z))
  
  ## color as part of `x' (data.frame or list):
  if(!is.null(d <- dim(x)) && (length(d) == 2) && (d[2] >= 4))
    color <- x[,4]
  else if(is.list(x) && !is.null(x$color))
    color <- x$color
  
  ## convert 'anything' -> vector
  xyz <- xyz.coords(x=x, y=y, z=z, xlab=xlabel, ylab=ylabel, zlab=zlabel,
                    log=log)
  if(is.null(xlab)) { xlab <- xyz$xlab; if(is.null(xlab)) xlab <- "" }
  if(is.null(ylab)) { ylab <- xyz$ylab; if(is.null(ylab)) ylab <- "" }
  if(is.null(zlab)) { zlab <- xyz$zlab; if(is.null(zlab)) zlab <- "" }
  
  if(length(color) == 1)
    color <- rep(color, length(xyz$x))
  else if(length(color) != length(xyz$x))
    stop("length(color) ", "must be equal length(x) or 1")
  
  angle <- (angle %% 360) / 90
  yz.f <- scale.y * abs(if(angle < 1) angle else if(angle > 3) angle - 4 else 2 - angle)
  yx.f <- scale.y * (if(angle < 2) 1 - angle else angle - 3)
  if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
    temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
    temp <- xlab;  xlab <- ylab;   ylab <- temp
    temp <- xlim;  xlim <- ylim;   ylim <- temp
  }
  angle.1 <- (1 < angle && angle < 2) || angle > 3
  angle.2 <- 1 <= angle && angle <= 3
  dat <- cbind(as.data.frame(xyz[c("x","y","z")]), col = color)
  
  n <- nrow(dat);
  y.range <- range(dat$y[is.finite(dat$y)])
  
  ### 3D-highlighting / colors / sort by y
  if(type == "p" || type == "h") {
    y.ord <- rev(order(dat$y))
    dat <- dat[y.ord, ]
    if(length(pch) > 1)
      if(length(pch) != length(y.ord))
        stop("length(pch) ", "must be equal length(x) or 1")
    else pch <- pch[y.ord]
    daty <- dat$y
    daty[!is.finite(daty)] <- mean(daty[is.finite(daty)])
    if(highlight.3d && !(all(diff(daty) == 0)))
      dat$col <- rgb(seq(0, 1, length = n) * (y.range[2] - daty) / diff(y.range), g=0, b=0)
  }
  
  ### optim. axis scaling
  p.lab <- par("lab")
  ## Y
  y.range <- range(dat$y[is.finite(dat$y)], ylim)
  y.prty <- pretty(y.range, n = lab[2],
                   min.n = max(1, min(.5 * lab[2], p.lab[2])))
  y.scal <- round(diff(y.prty[1:2]), digits = 12)
  y.add <- min(y.prty)
  dat$y <- (dat$y - y.add) / y.scal
  y.max <- (max(y.prty) - y.add) / y.scal
  
  x.range <- range(dat$x[is.finite(dat$x)], xlim)
  x.prty <- pretty(x.range, n = lab[1],
                   min.n = max(1, min(.5 * lab[1], p.lab[1])))
  x.scal <- round(diff(x.prty[1:2]), digits = 12)
  dat$x <- dat$x / x.scal
  x.range <- range(x.prty) / x.scal
  x.max <- ceiling(x.range[2])
  x.min <-   floor(x.range[1])
  if(!is.null(xlim)) {
    x.max <- max(x.max, ceiling(xlim[2] / x.scal))
    x.min <- min(x.min,   floor(xlim[1] / x.scal))
  }
  x.range <- range(x.min, x.max)
  ## Z
  z.range <- range(dat$z[is.finite(dat$z)], zlim)
  z.prty <- pretty(z.range, n = lab.z,
                   min.n = max(1, min(.5 * lab.z, p.lab[2])))
  z.scal <- round(diff(z.prty[1:2]), digits = 12)
  dat$z <- dat$z / z.scal
  z.range <- range(z.prty) / z.scal
  z.max <- ceiling(z.range[2])
  z.min <-   floor(z.range[1])
  if(!is.null(zlim)) {
    z.max <- max(z.max, ceiling(zlim[2] / z.scal))
    z.min <- min(z.min,   floor(zlim[1] / z.scal))
  }
  z.range <- range(z.min, z.max)
  
  ### init graphics
  plot.new()
  if(angle.2) {x1 <- x.min + yx.f * y.max; x2 <- x.max}
  else        {x1 <- x.min; x2 <- x.max + yx.f * y.max}
  plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max))
  temp <- strwidth(format(rev(y.prty))[1], cex = cex.axis/par("cex"))
  if(angle.2) x1 <- x1 - temp - y.margin.add
  else        x2 <- x2 + temp + y.margin.add
  plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max))
  if(angle > 2) par("usr" = par("usr")[c(2, 1, 3:4)])
  usr <- par("usr") # we have to remind it for use in closures
  title(main, sub, ...)
  
  ### draw axis, tick marks, labels, grid, ...
  xx <- if(angle.2) c(x.min, x.max) else c(x.max, x.min)
  if(grid) {
    ## grids
    ###################
    # XY wall
    i <- x.min:x.max;
    segments(i, z.min, i + (yx.f * y.max), yz.f * y.max + z.min,
             col = col.grid, lty = lty.grid);
    
    i <- 0:y.max;
    segments(x.min + (i * yx.f), i * yz.f + z.min,
             x.max + (i * yx.f), i * yz.f + z.min,
             col = col.grid, lty = lty.grid);
    
    ######################
    # XZ wall
    # verticle lines
    temp <- yx.f * y.max;
    temp1 <- yz.f * y.max;
    i <- (x.min + temp):(x.max + temp);
    segments(i, z.min + temp1, i, z.max + temp1,
             col = col.grid, lty = lty.grid);
    
    # horizontal lines
    i <- (z.min + temp1):(z.max + temp1);
    segments(x.min + temp, i, x.max + temp, i,
             col = col.grid, lty = lty.grid)
    
    
    ##################
    # YZ wall
    # horizontal lines
    i <- xx[2]:x.min;
    mm <- z.min:z.max;
    segments(i, mm, i + temp, mm + temp1,
             col = col.grid, lty = lty.grid);
    # verticle lines
    i <- 0:y.max;
    segments(x.min + (i * yx.f), i * yz.f + z.min,
             xx[2] + (i * yx.f), i * yz.f + z.max,
             col = col.grid, lty = lty.grid)
    
    
    # make the axis into solid line
    segments(x.min, z.min, x.min + (yx.f * y.max), yz.f * y.max + z.min,
             col = col.grid, lty = lty.hide);
    segments(x.max, z.min, x.max + (yx.f * y.max), yz.f * y.max + z.min,
             col = col.axis, lty = lty.hide);
    segments(x.min + (y.max * yx.f), y.max * yz.f + z.min,
             x.max + (y.max* yx.f), y.max * yz.f + z.min,
             col = col.grid, lty = lty.hide);
    segments(x.min + temp, z.min + temp1, x.min + temp, z.max + temp1,
             col = col.grid, lty = lty.hide);
    segments(x.max + temp, z.min + temp1, x.max + temp, z.max + temp1,
             col = col.axis, lty = lty.hide);
    segments(x.min + temp, z.max + temp1, x.max + temp, z.max + temp1,
             col = col.axis, lty = lty.hide);
    segments(xx[2], z.max, xx[2] + temp, z.max + temp1,
             col = col.axis, lty = lty.hide);
  }
  if(axis) {
    if(tick.marks) { ## tick marks
      xtl <- (z.max - z.min) * (tcl <- -par("tcl")) / 50
      ztl <- (x.max - x.min) * tcl / 50
      mysegs <- function(x0,y0, x1,y1)
        segments(x0,y0, x1,y1, col=col.axis, lty=lty.axis)
      ## Y
      i.y <- 0:y.max
      mysegs(yx.f * i.y - ztl + xx[1], yz.f * i.y + z.min,
             yx.f * i.y + ztl + xx[1], yz.f * i.y + z.min)
      ## X
      i.x <- x.min:x.max
      mysegs(i.x, -xtl + z.min, i.x, xtl + z.min)
      ## Z
      i.z <- z.min:z.max
      mysegs(-ztl + xx[2], i.z, ztl + xx[2], i.z)
      
      if(label.tick.marks) { ## label tick marks
        las <- par("las")
        mytext <- function(labels, side, at, ...)
          mtext(text = labels, side = side, at = at, line = -.5,
                col=col.lab, cex=cex.axis, font=font.lab, ...)
        ## X
        if(is.null(x.ticklabs))
          x.ticklabs <- format(i.x * x.scal)
        mytext(x.ticklabs, side = 1, at = i.x)
        ## Z
        if(is.null(z.ticklabs))
          z.ticklabs <- format(i.z * z.scal)
        mytext(z.ticklabs, side = if(angle.1) 4 else 2, at = i.z,
               adj = if(0 < las && las < 3) 1 else NA)
        ## Y
        temp <- if(angle > 2) rev(i.y) else i.y ## turn y-labels around
        if(is.null(y.ticklabs))
          y.ticklabs <- format(y.prty)
        else if (angle > 2)
          y.ticklabs <- rev(y.ticklabs)
        text(i.y * yx.f + xx[1],
             i.y * yz.f + z.min, y.ticklabs,
             pos=if(angle.1) 2 else 4, offset=1,
             col=col.lab, cex=cex.axis/par("cex"), font=font.lab)
      }
    }
    
    ## axis and labels
    
    mytext2 <- function(lab, side, line, at)
      mtext(lab, side = side, line = line, at = at, col = col.lab,
            cex = cex.lab, font = font.axis, las = 0)
    ## X
    lines(c(x.min, x.max), c(z.min, z.min), col = col.axis, lty = lty.axis)
    mytext2(xlab, 1, line = 1.5, at = mean(x.range))
    ## Y
    lines(xx[1] + c(0, y.max * yx.f), c(z.min, y.max * yz.f + z.min),
          col = col.axis, lty = lty.axis)
    mytext2(ylab, if(angle.1) 2 else 4, line= 0.5, at = z.min + y.max * yz.f)
    
    ## Z
    lines(xx[c(2,2)], c(z.min, z.max), col = col.axis, lty = lty.axis)
    mytext2(zlab, if(angle.1) 4 else 2, line= 1.5, at = mean(z.range))
    
  }
  
  ### plot points
  x <- dat$x + (dat$y * yx.f)
  z <- dat$z + (dat$y * yz.f)
  col <- as.character(dat$col)
  if(type == "h") {
    z2 <- dat$y * yz.f + z.min
    segments(x, z, x, z2, col = col, cex = cex.symbols, lty = lty.hplot, ...)
    points(x, z, type = "p", col = col, pch = pch, cex = cex.symbols, ...)
  }
  else points(x, z, type = type, col = col, pch = pch, cex = cex.symbols, ...)
  
  ### box-lines in front of points (overlay)
  if(axis && box) {
    lines(c(x.min, x.max), c(z.max, z.max),
          col = col.axis, lty = lty.axis)
    lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + z.max,
          col = col.axis, lty = lty.axis)
    lines(xx[c(1,1)], c(z.min, z.max), col = col.axis, lty = lty.axis)
  }
  
  
  # par(mem.par) # we MUST NOT set the margins back
  ### Return Function Object
  ob <- ls() ## remove all unused objects from the result's enviroment:
  rm(list = ob[!ob %in% c("angle", "mar", "usr", "x.scal", "y.scal", "z.scal", "yx.f",
                          "yz.f", "y.add", "z.min", "z.max", "x.min", "x.max", "y.max",
                          "x.prty", "y.prty", "z.prty")])
  rm(ob)
  invisible(list(
    xyz.convert = function(x, y=NULL, z=NULL) {
      xyz <- xyz.coords(x, y, z)
      if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
        temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
      }
      y <- (xyz$y - y.add) / y.scal
      return(list(x = xyz$x / x.scal + yx.f * y,
                  y = xyz$z / z.scal + yz.f * y))
    },
    points3d = function(x, y = NULL, z = NULL, type = "p", ...) {
      xyz <- xyz.coords(x, y, z)
      if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
        temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
      }
      y2 <- (xyz$y - y.add) / y.scal
      x <- xyz$x / x.scal + yx.f * y2
      y <- xyz$z / z.scal + yz.f * y2
      mem.par <- par(mar = mar, usr = usr)
      on.exit(par(mem.par))
      if(type == "h") {
        y2 <- z.min + yz.f * y2
        segments(x, y, x, y2, ...)
        points(x, y, type = "p", ...)
      }
      else points(x, y, type = type, ...)
    },
    plane3d = function(Intercept, x.coef = NULL, y.coef = NULL,
                       lty = "dashed", lty.box = NULL, ...){
      if(!is.atomic(Intercept) && !is.null(coef(Intercept))) Intercept <- coef(Intercept)
      if(is.null(lty.box)) lty.box <- lty
      if(is.null(x.coef) && length(Intercept) == 3){
        x.coef <- Intercept[if(angle > 2) 3 else 2]
        y.coef <- Intercept[if(angle > 2) 2 else 3]
        Intercept <- Intercept[1]
      }
      mem.par <- par(mar = mar, usr = usr)
      on.exit(par(mem.par))
      x <- x.min:x.max
      ltya <- c(lty.box, rep(lty, length(x)-2), lty.box)
      x.coef <- x.coef * x.scal
      z1 <- (Intercept + x * x.coef + y.add * y.coef) / z.scal
      z2 <- (Intercept + x * x.coef +
               (y.max * y.scal + y.add) * y.coef) / z.scal
      segments(x, z1, x + y.max * yx.f, z2 + yz.f * y.max, lty = ltya, ...)
      y <- 0:y.max
      ltya <- c(lty.box, rep(lty, length(y)-2), lty.box)
      y.coef <- (y * y.scal + y.add) * y.coef
      z1 <- (Intercept + x.min * x.coef + y.coef) / z.scal
      z2 <- (Intercept + x.max * x.coef + y.coef) / z.scal
      segments(x.min + y * yx.f, z1 + y * yz.f,
               x.max + y * yx.f, z2 + y * yz.f, lty = ltya, ...)
    },
    
    wall3d = function(Intercept, x.coef = NULL, y.coef = NULL,
                      lty = "dashed", lty.box = NULL, ...){
      if(!is.atomic(Intercept) && !is.null(coef(Intercept))) Intercept <- coef(Intercept)
      if(is.null(lty.box)) lty.box <- lty
      if(is.null(x.coef) && length(Intercept) == 3){
        x.coef <- Intercept[if(angle > 2) 3 else 2]
        y.coef <- Intercept[if(angle > 2) 2 else 3]
        Intercept <- Intercept[1]
      }
      mem.par <- par(mar = mar, usr = usr)
      on.exit(par(mem.par))
      x <- x.min:x.max
      ltya <- c(lty.box, rep(lty, length(x)-2), lty.box)
      x.coef <- x.coef * x.scal
      z1 <- (Intercept + x * x.coef + y.add * y.coef) / z.scal
      z2 <- (Intercept + x * x.coef +
               (y.max * y.scal + y.add) * y.coef) / z.scal
      segments(x, z1, x + y.max * yx.f, z2 + yz.f * y.max, lty = ltya, ...)
      y <- 0:y.max
      ltya <- c(lty.box, rep(lty, length(y)-2), lty.box)
      y.coef <- (y * y.scal + y.add) * y.coef
      z1 <- (Intercept + x.min * x.coef + y.coef) / z.scal
      z2 <- (Intercept + x.max * x.coef + y.coef) / z.scal
      segments(x.min + y * yx.f, z1 + y * yz.f,
               x.max + y * yx.f, z2 + y * yz.f, lty = ltya, ...)
    },
    box3d = function(...){
      mem.par <- par(mar = mar, usr = usr)
      on.exit(par(mem.par))
      lines(c(x.min, x.max), c(z.max, z.max), ...)
      lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + z.max, ...)
      lines(c(0, y.max * yx.f) + x.min, c(0, y.max * yz.f) + z.max, ...)
      lines(c(x.max, x.max), c(z.min, z.max), ...)
      lines(c(x.min, x.min), c(z.min, z.max), ...)
      lines(c(x.min, x.max), c(z.min, z.min), ...)
    }
  ))
}

.set.mSet <- function(mSetObj=NA){
  if(.on.public.web){
    mSet <<- mSetObj;
    return (1);
  }
  return(mSetObj);
}

# Copy the function to a new object
PlotPCA3DScoreImg <- MetaboAnalystR::PlotPCA3DScoreImg

# Now edit the function
fix(PlotPCA3DScoreImg)
#Change Plot3D(..., color = cols, ...) to Plot3D(..., color = all.cols, ...) 

###################################################################################################

# Create a 3D PCA score plot, currently not working
mSet<-PlotPCA3DScoreImg(mSet, "pca_score3d_0_", "png", 300, width=NA, 1,2,3, 40)

# Create a PCA loadings Plots
mSet<-PlotPCALoading(mSet, "pca_loading_0_", "png", 300, width=NA, 1,2);

# Create a PCA Biplot
mSet<-PlotPCABiplot(mSet, "pca_biplot_0_", format = "png", dpi = 300, width=NA, 1, 2)

#Partial Least Squares - Discriminant Analysis (PLS-DA)
mSet<-PLSR.Anal(mSet, reg=TRUE)
mSet<-PlotPLSPairSummary(mSet, "pls_pair_0_", "png", 300, width=NA, 5)
mSet<-PlotPLS2DScore(mSet, "pls_score2d_0_", "png", 300, width=NA, 1,2,0.95,1,0)
mSet<-PlotPLS3DScoreImg(mSet, "pls_score3d_0_", "png", 300, width=NA, 1,2,3, 40)
mSet<-PlotPLSLoading(mSet, "pls_loading_0_", "png", 300, width=NA, 1, 2);
mSet<-PLSDA.CV(mSet, "5", 5,5, "Q2")
mSet<-PlotPLS.Classification(mSet, "pls_cv_0_", "png", 300, width=NA)
mSet<-PlotPLS.Imp(mSet, "pls_imp_0_", "png", 300, width=NA, "vip", "Comp. 1", 15, FALSE)
mSet<-PLSDA.Permut(mSet, 100, "accu")
mSet<-PlotPLS.Permutation(mSet, "pls_perm_1_", "png", 300, width=NA)

#Hierarchical CLustering: Dendogram
mSet<-PlotHCTree(mSet, "tree_0_", format = "png", dpi=300, width=NA, "euclidean", "complete")
?hclust

#Hierarchical CLustering: Heatmaps
mSet<-PlotStaticHeatMap(mSet, "heatmap_0_", "png", 300, width=NA, dataOpt = "norm", scaleOpt = "row",
                  smplDist = "euclidean", clstDist = "complete", palette = "bwm", fzCol = 8,
                  fzRow = 8, viewOpt = "detail", rowV = T, colV = T, var.inx = NULL,
                  border = T, grp.ave = F, show.legend = T, show.annot.legend = T, includeRowNames = T)
