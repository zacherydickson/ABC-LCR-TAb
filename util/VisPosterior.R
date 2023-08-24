args <- commandArgs(trailingOnly = T)

if(length(args) < 1){
    stop("Usage: Vis....R inFile.res [alpha = 0.05] (parses inFile.log also creates inFile.pdf and inFile.eval")
}

resFile <- args[1]
#threads <- as.numeric(args[2])
#if(is.na(threads)){
#    threads <- 1
#}
CIalpha <- as.numeric(args[2])
if(is.na(CIalpha)){
    CIalpha=0.05
}
logFile <- sub("res$","log",resFile)
pdfFile <- sub("res$","pdf",resFile)
evalFile <- sub("res$","eval",resFile)
ciFile <- sub("res$","ci",resFile)

if(any(resFile == c(logFile,pdfFile,evalFile))){
    stop("Could not generate unique log, pdf, and eval files from the inFile, does it end in 'res'?")
}

require("vioplot")
require("Rcpp")
CPPSource <- "/home/zac/scripts/LCR/TemporalOrder/Primates/ABC2/util/VisPosterior.cpp"
sourceCpp(CPPSource)

ColourPalette <- c(SpanishGrey="#989898",SeaGreen="#388659",Coral="#FF7F50",SteelBlue="#4682B4",Veronica="#A020F0",LightBlue="#ADD8E6")
#CIPalette = setNames(c("grey","#7A9CC6","#B3D2B2","#FFFD98"),c("0","90","95","99"))
CIPalette = setNames(ColourPalette[c("Coral","Veronica","SteelBlue","LightBlue","SpanishGrey")],c("0","50","90","95","99"))

######## FUNCTIONS ##############

parseLog <- function(file){
    lines <- scan(file,what="character",sep="\n")
    #Reduce to iteration, and swap indication and reult lines
    lines <- lines[grepl("===Iteration|Swap|with p=",lines)]
    #Remove any indication of non-main chain swaps
    nonMainSwaps <- grep("Swap between chains [^0]",lines)
    if(length(nonMainSwaps) > 0){
        lines <- lines[-c(nonMainSwaps,nonMainSwaps+1)]
    }
    #Remove Swap indicating lines, no longer needed
    lines <- lines[!grepl("Swap",lines)]
    AccSwaps <- grep("Accept with",lines)
    InitialAcceptCount <- sapply(strsplit(lines[1],"\\(|/|\\)",perl=T),function(x){as.numeric(x[length(x)-1])})
    AcceptCountAtSwap <- sapply(strsplit(lines[AccSwaps-1],"\\(|/|\\)",perl=T),function(x){as.numeric(x[length(x)-1])})
    return (AcceptCountAtSwap + InitialAcceptCount)
}

#Effective Sample Size = n / (1 + 2*Σ_k cor@lag(k))
lagcor <- function(x,k){
    idx <- (1+k):length(x)
    cor(x[idx],x[idx-k])
}


windowedSumDeviation <- function(x,w){
    z<-sapply(seq(1,length(x),by=w),function(i){
                  y <- x[i:(i+w-1)];
                  m <- mean(y);
                  sum(y - m)
        });
    z[!is.na(z)]
}

kneedle <- function(x,guess=length(x),bPlot=FALSE,...){
    if(is.na(guess) || guess > length(x) / 2){
        guess = length(x)/2
    }
    x = x[1:(guess*2)]
    tmp <- lowess(x,...)
    y = tmp$y
    y0 = y[1]
    yn= y[length(y)]
    x0=0
    xn=length(y)-1
    m = (yn-y0)/(xn-x0)
    b = y0
    d = abs(y - (m*(x0:xn) + b))
    B <- which.max(d) + 1
    #Get the equation of line for the perpendicular at point B
    #m2 = -1/m
    #b2 = tmp$y[B] - m2*tmp$x[B]
    #x_int = (b - b2) / (m2-m)
    #y_int = x_int*m + b
    #segments(tmp$x[B],tmp$y[B],x_int,y_int,col=guidecol)
    if(bPlot){
        guidecol <- rgb(0.5,0.5,0.5,0.5)
        plot(x,main="Kneedle Point Estimation",xlab="",ylab="",type="l")
        lines(tmp$x,tmp$y,col=guidecol)
        segments(tmp$x[1],tmp$y[1],tmp$x[nrow(df)],tmp$y[nrow(df)],col=guidecol)
        abline(v=B,lwd=3,col="red")
    }
    return(B)
}

densityJitter <- function(x,a=0,b=1){
    jitter <- rep(0,length(x))
    bNA <- is.na(x)
    x <- x[!bNA]
    dens <- density(x)
    densIdx <- sapply(x,function(z){which.min(abs(dens$x - z))})
    jitter[!bNA] <- dens$y[densIdx] / max(dens$y) * runif(length(x),a,b)
    jitter
}


mESSBurninEst <- function(p,a=0.05,e=0.05){
    2^(2/p) * pi / (p*gamma(p/2))^(2/p) * qchisq(1-a,p) / e^2
}

mBM <- function(df){
    n <- nrow(df)
    a_n <- floor(sqrt(n));
    b_n <- n / a_n
    batchFactor <- cut(1:n,breaks=a_n)
    batches <- split(df,batchFactor)
    M = apply(df,2,mean)
    BatchMeanDeviates <- function(x){
        mat <- apply(x,2,mean) - M;
        mat %*% t(mat)
    }
    result = Reduce("+",lapply(batches,BatchMeanDeviates))
    result * b_n / (a_n - 1)
}

eSS <- function(df){
    sampleCov <- cov(df)
    mBMEst <- mBM(df)
    nrow(df) * (det(sampleCov)/det(mBMEst))^(1/ncol(df))
}

modeEst <- function(x){
   dens <- density(x[!is.na(x)])
   dens$x[which.max(dens$y)]
}

#### OPTIMIZATION FUNCTIONS

gradientEst <- function(point,dens,mat,bw,factor=10^-2){
    p <- length(point)
    jitPoint = point + apply(mat,2,function(x){min(abs(diff(sort(unique(x)))))}) * factor
    testPoints <- matrix(point,nrow=p,ncol=p,byrow=T)
    diag(testPoints) <- jitPoint
    testDens <- apply(testPoints,1,calculateMVDensity,data=mat,bw=bw)
    (testDens - dens) / diag(testPoints - point)
}

findUpperMult <- function(point,val,jitSize,gradient,FUN,...){
    mult <- 1
    repeat {
        point2 <- mult*jitSize*gradient + point 
        val2 <- FUN(point2,...)
        if(val2 <= val){
            break
        }
        mult <- mult * 2
    }
    mult
}

generalizedGoldenSearch <- function(l,r,FUN,tolProp=10^-6,...){
    phi <- (1+sqrt(5))/2
    lVal <- FUN(l,...)
    rVal <- FUN(r,...)
    p1 <- l + (r-l) / (1+phi)
    p1Val <- FUN(p1,...)
    tol <- (r-l) * tolProp
    repeat {
        p2 <- l + (r - p1)
        p2Val <- FUN(p2,...)
        #Define p1 to always be the left of the two probe points
        if(p1 > p2){
            tmp <- p1
            p1 <- p2
            p2 <- tmp
            tmp <- p1Val
            p1Val <- p2Val
            p2Val <- tmp
        }
        if(p2 - l < tol){
            break;
        }
        if(p2Val > p1Val){
            l <- p1
            lVal <- p1Val
            p1 <- p2
            p1Val <- p2Val
        } else {
            r <- p2
            rVal <- p2Val
        }
    }
    c(p1,p1Val)
}

optimizeMVMode <- function(point,val,mat,bw,tolProp=10^-4,gssTolProp=10^-6,jitFactor=10^-2){
    jitSize <- apply(mat,2,function(x){min(abs(diff(sort(unique(x)))))}) * jitFactor
    tol <- sqrt(tolProp^2*length(point))
    repeat {
        gradient <-  gradientEst(point,val,mat=mat,bw=bw)
        gradMag <- sqrt(sum(gradient^2))
        if(gradMag < tol){
            break
        }
        r <- findUpperMult(point,val,jitSize,gradient,calculateMVDensity,data=mat,bw=bw) 
        if(r <= 1){
            break
        }
        res <- generalizedGoldenSearch(1,r,function(mult,...){p <- point + mult*jitSize*gradient; calculateMVDensity(p,...)},data=mat,bw=bw)
        if(res[2] < val){
            break
        }
        inc <- res[2] - val
        val <- res[2]
        point <- point + res[1]*jitSize*gradient
        if(inc < tolProp){
            break
        }
    }
    list(point=point,val=val)
}

##MULTIVARIATE MODE AND CREDIBILITY REGION FUNCTIONS

getStandardizedMatrix <- function(df){
    stdInfo <- apply(df,2,function(x){c(mean(x),sd(x))})
    mat <- sapply(names(df),function(i){(df[,i]- stdInfo[1,i])/stdInfo[2,i]})
    #Take an initial guess of the mode at the point which is the mode for each individual
    #component, also get the bandwidth estimate
    densInfo <- apply(mat,2,function(x){dens <- density(x); c(dens$x[which.max(dens$y)],dens$bw)})
    bw <- mean(densInfo[2,])
    uniMode <- densInfo[1,]
    #The general info
    obj <- list(stdInfo=stdInfo,mat=mat,bw=bw,uniMode=uniMode)
    #Add some metadata
    attr(obj,"class") <- "StandardizedMatrix"
    obj$n = nrow(df)
    obj$names <- colnames(stdInfo)
    obj$m = length(obj$names)
    #Add Fields to be added
    obj$mode <- numeric(0)
    obj$modalDensity <- numeric(0)
    obj$density <- numeric(0)
    obj$modalDist <- numeric(0)
    obj$credRadii <- numeric(0)
    obj
}

getDensityEst <- function(stdMat) {
    if(attr(stdMat,"class") != "StandardizedMatrix"){
        stop("Attempt to getDensityEst from non-StandardizedMatrix object")
    }
    stdMat$density <- apply(stdMat$mat,1,calculateMVDensity,data=stdMat$mat,bw=stdMat$bw)
    stdMat
}

estimateMvMode <- function(stdMat,...){
    if(attr(stdMat,"class") != "StandardizedMatrix"){
        stop("Attempt to estimateMvMode from non-StandardizedMatrix object")
    }
    if(length(stdMat$density) == 0){
        stop("Attempt to estimateMvMode to improperly initialized StandardizedMatrix")
    }
    res <- optimizeMVMode(stdMat$mat[which.max(stdMat$density),],max(stdMat$density),stdMat$mat,stdMat$bw,...)
    stdMat$mode <- res$point
    stdMat$modalDensity <- res$val
    stdMat$modalDist <- apply(stdMat$mat,1,function(x){sqrt(sum((x - stdMat$mode)^2))})
    stdMat
}

#subsetStdMat can alos be used to reorder columns
subsetStdMat <- function(stdMat,names){
    if(attr(stdMat,"class") != "StandardizedMatrix"){
        stop("Attempt to subsetStdMat from non-StandardizedMatrix object")
    }
    isValid = prod(sapply(c("uniMode","mode","stdInfo","mat"),function(i){length(stdMat[[i]])}))
    if(!isValid){
        stop("Attempt to subsetStdMat before finding the mode")
    }
    if(length(names) == stdMat$m && all(stdMat$names == names)){
        #All names present and in order
        return(stdMat)
    }
    if(length(names) > stdMat$m || any(! names %in% stdMat$names)){
        stop("Attempt to subsetStdMat with non-existant names")
    }
    stdMat$names <- names
    stdMat$m <- length(stdMat$names)
    for(member in c("uniMode","mode")){
        stdMat[[member]] <- stdMat[[member]][stdMat$names]
    }
    for(member in c("stdInfo","mat")){
        stdMat[[member]] <- stdMat[[member]][,stdMat$names]
    }
    stdMat
}

addFixedNames <- function(stdMat,isFixed,mu){
    if(attr(stdMat,"class") != "StandardizedMatrix"){
        stop("Attempt to addFixedNames from non-StandardizedMatrix object")
    }
    if(length(stdMat$mode) == 0){
        stop("Attempt to addFixedNames to StandardizedMatrix before finding mode")
    }
    if(!sum(isFixed)){
        return(stdMat)
    }
    #Add info as necessary
    fixedMu <- mu[isFixed]
    stdMat$uniMode <- c(stdMat$uniMode,fixedMu)
    stdMat$stdInfo <- cbind(stdMat$stdInfo,sapply(fixedMu,c,0))
    stdMat$mat <- cbind(stdMat$mat,sapply(fixedMu,function(i){rep(0,stdMat$n)}))
    stdMat$mode <- c(stdMat$mode,fixedMu)
    #Reorder
    stdMat$names <- names(stdMat$mode)
    stdMat$m <- length(stdMat$names)
    subsetStdMat(stdMat,names(isFixed)) 
}

getRescaledMode <- function(stdMat){
    if(attr(stdMat,"class") != "StandardizedMatrix"){
        stop("Attempt to get RescaledMode from non-StandardizedMatrix object")
    }
    if(length(stdMat$mode) == 0){
        stop("Attempt to get RescaledMode from improperly initialized StandardizedMatrix")
    }
    stdMat$mode * stdMat$stdInfo[2,] + stdMat$stdInfo[1,]
}

ColourByModalDistance <- function(stdMat,ciPalette,default="black"){
    if(attr(stdMat,"class") != "StandardizedMatrix"){
        stop("Attempt to ColourByModalDistance with a non-StandardizedMatrix object")
    }
    if(length(stdMat$modalDist) == 0){
        stop("Attempt to ColourByModalDistance with improperly initialized StandardizedMatrix")
    }
    pointCol = rep("black",stdMat$n)
    ciPalette <- ciPalette[order(as.numeric(names(ciPalette)))]
    q = as.numeric(names(ciPalette))/100
    colIdx <- rowSums(sapply(quantile(stdMat$modalDist,q),"<",stdMat$modalDist))
    pointCol[colIdx > 0] <- ciPalette[colIdx]
    pointCol
}

estimateCredibilityRadii <- function(stdMat,cutoffs = c(0.5,0.9,0.95,0.99)){
    if(attr(stdMat,"class") != "StandardizedMatrix"){
        stop("Attempt to EstimateCredibilityRadii with a non-StandardizedMatrix object")
    }
    if(length(stdMat$modalDist) == 0 || length(stdMat$density) == 0){
        stop("Attempt to EstimateCredibilityRadii with improperly initialized StandardizedMatrix")
    }
    distOrder <- order(stdMat$modalDist)
    cumProp <- cumsum(stdMat$density[distOrder]) / sum(stdMat$density)
    idx <- distOrder[sapply(cutoffs,function(x){min(which(cumProp > x))})]
    stdMat$credRadii <- setNames(c(0,stdMat$modalDist[idx]),c(0,cutoffs*100))
    stdMat
}


ColourByCredibility <- function(stdMat,ciPalette,defaultCol="black"){
    if(attr(stdMat,"class") != "StandardizedMatrix"){
        stop("Attempt to ColourByCredibility  with a non-StandardizedMatrix object")
    }
    if(length(stdMat$credRadii) == 0){
        stop("Attempt to ColourByCredibility with improperly initialized StandardizedMatrix")
    }
    #Calculate cumulative densities and colour the points
    pointCol = rep("black",stdMat$n)
    ciPalette <- ciPalette[order(as.numeric(names(ciPalette)))]
    colIdx <- rowSums(sapply(stdMat$credRadii[names(ciPalette)],"<",stdMat$modalDist))
    pointCol[colIdx > 0] <- ciPalette[colIdx]
    pointCol
}

euclidianDist <- function(mat,x){
    mat |> sweep(2,x,"-") |> (`^`)(2) |> rowSums() |> sqrt()
}

#Assume the point most extreme in the CR for one dim is near the edge of CR
#Assume the nearest point outside the CR, with a more extreme dim value is still close
#Linearaly interpolate a dim value at the CR limit given a line
findCIExtremePoint <- function(dim,mat, nmat, cpI, cpO, alpha = 0.05, maximum = T){
    which.ext = which.max
    compFunc = (`>`)
    if(!maximum){
        which.ext <- which.min
        compFunc <- (`<`)
    }
    innerIdx <- which.ext(mat[,dim])
    point <- mat[innerIdx,]
    extIdx <- which(compFunc(nmat[,dim],point[dim]))
    if(length(extIdx) == 0){
        return(point)
    }
    nmat <- nmat[extIdx,,drop=FALSE]
    idx <- nmat |> euclidianDist(point) |> which.min()
    outerIdx <- extIdx[idx]
    weight <- (cpO[outerIdx] - 1 + alpha) / (cpO[outerIdx] - cpI[innerIdx])
    limitPoint <- point
    limitPoint[dim] <- nmat[idx,dim,drop=FALSE] * (1-weight) + point[dim] * weight
    limitPoint
}

getRescaledCI <- function(stdMat,alpha=0.05){
    if(attr(stdMat,"class") != "StandardizedMatrix"){
        stop("Attempt to getRescaledCI from a non-StandardizedMatrix object")
    }
    if(length(stdMat$modalDist) == 0 || length(stdMat$density) == 0){
        stop("Attempt to getRescaledCI from improperly initialized StandardizedMatrix")
    }
    if(alpha > 0){
        distOrder <- order(stdMat$modalDist)
        cumProp <- cumsum(stdMat$density[distOrder]) / sum(stdMat$density)
        #Identitify the first point outside of the credibility region
        idx <- min(which(cumProp >=  1 - alpha))
        #Subset the matrix to be only points inside the credibility region
        if(idx == 1){ #The closest point to the mode is outside of the credibility region
            mat <- rbind(stdMat$mode,stdMat$mode)
            nmat <- stdMat$mat
            cpI <- c(0,0)
            cpO <- cumProp
        } else {
            inReg <- 1:(idx-1)
            mat <- stdMat$mat[distOrder[inReg],,drop=F]
            nmat <- stdMat$mat[distOrder[-inReg],,drop=F]
            cpI <- cumProp[inReg]
            cpO <- cumProp[-inReg]
        }
        #For each dimension, find the point at the extremes of the credibility region
        #Then find the point outside the region, and which has a more extreme value for that dim
        #Add a new point to the credibility interval which has the same non-dim values as
        # the original, but has dim value weighted by distance beyond the credibility
        # region
        mat <- mapply(findCIExtremePoint,
                      rep(1:stdMat$m,each=2),rep(c(T,F),stdMat$m),
                      MoreArgs=list(mat=mat,nmat =nmat,cpI=cpI,cpO=cpO,alpha=alpha)) |>
            t() |> rbind(mat)
    } else {
        mat <- stdMat$mat
    }
    mat |> apply(2,range) |>
        sweep(2,stdMat$stdInfo[2,],"*") |> sweep(2,stdMat$stdInfo[1,],"+")
}


### PLOTTING FUNCTIONS

vioplotWPoints <- function(data,pointCol=NULL,mode=NULL,order=1:(data),...){
    vioplot(data,...,side="left",plotCentre="line")
    if(is.null(pointCol)){
        pointCol = rep("black",nrow(data))
    }
    if(is.null(ncol(data))){ #input is a vector
        dJit <- 1.05 + densityJitter(data,0,0.4)
        points(dJit,data,col=pointCol,pch=19)
        if(!is.null(mode)){
            segments(1,mode,1.5,mode,col="black",lwd=2)
        }
        return(NULL)
    }
    dJit <- lapply(data,densityJitter,b=0.4)
    dJit <- mapply("+",dJit,1:ncol(data)+0.05,SIMPLIFY=F)
    mapply(function(x,y){points(x,y,col=pointCol,pch=19)},dJit,data)
    if(!is.null(mode)){
        segments(1:length(data),mode,1:length(data)+ 0.5,mode,col="black",lwd=2)
    }

    return(NULL)
}

#by default no axis are shown
#axis can be a vector containing the values 1:4, an axis which be shown on each side
#present
scatterPlotWiCred <- function(stdMat,axes=NULL,col=rep("black",stdMat$n),...){
    if(attr(stdMat,"class") != "StandardizedMatrix"){
        stop("Attempt to scatterPlotWiCred with a non-StandardizedMatrix object")
    }
    if(length(stdMat$modalDist) == 0 || length(stdMat$credRadii) == 0){
        stop("Attempt to scatterPlotWiCred before estimating the mode or estimating credibility")
    }
    if(stdMat$m > 2){
        stop("Attempt to build a scatterplot for a Standardized Matrix which hasn't been subset to two columns")
    }
    mat <- stdMat$mat |> sweep(2,stdMat$stdInfo[2,],"*") |> sweep(2,stdMat$stdInfo[1,],"+") 
    order <- order(stdMat$modalDist,decreasing=T)

    mat <- mat[order,]
    col <- col[(1:stdMat$n -1) %% length(col) +1]
    col <- col[order]
    plot(mat,ann="F",xaxt="n",yaxt="n",pch=20,col=col,...)
    garbage <- lapply(axes,axis)
}

specialHist <- function(x,horiz=F,...){
    obj <- hist(x,plot=F)
    if(length(obj$density) == 1){
        val <- x[1]
        plot(val,val,type="n",...)
        if(horiz){
            abline(h=val)
        } else {
            abline(v=val)
        }
        return(invisible())
    }
    if(horiz){
        w <- diff(obj$breaks)[1] 
        args <- list(...)
        if(hasArg("xlim")){
            args$xlim <- rev(args$xlim)
        } else {
            args$xlim <- rev(range(obj$density))
        }
        args$x=range(obj$density)
        args$y=range(obj$mids)+c(-w/2,w/2)
        args$type="n"
        do.call(plot,args)
        col <- "gray"
        if(hasArg("col")){
            col <- args$col
        }
        rect(0,obj$mids -w/2,obj$density,obj$mids +w/2,col=col)
    } else {
        obj <- hist(x,...)
    }
    invisible(obj)
}

scatterPlotMatrix <- function(stdMat,...){
    pointCol <- ColourByCredibility(StdMat,CIPalette)
    oldPar <- par()
    par(mfrow = c(stdMat$m,stdMat$m),mar=c(0.5,0.5,0.1,0.1),oma=c(4.1,4.1,2.1,0.1),las=2,xpd=T)
        #,cex=1.3)
    plotW <- (dev.size()[2] - sum(par()$omi[c(2,4)]))/stdMat$m -sum(par()$mai[c(2,4)])
    plotH <- (dev.size()[1] - sum(par()$omi[c(1,3)]))/stdMat$m -sum(par()$mai[c(1,3)])
    mat <- stdMat$mat |> sweep(2,stdMat$stdInfo[2,],"*") |> sweep(2,stdMat$stdInfo[1,],"+") 
    limits <- sweep(stdMat$mat,2,stdMat$stdInfo[2,],"*") |> sweep(2,stdMat$stdInfo[1,],"+") |>
              apply(2,pretty,min.n=3) |> sapply(range)
    mode <- stdMat$mode * stdMat$stdInfo[2,] + stdMat$stdInfo[1,]
    for(row in 1:stdMat$m){
        for(col in 0:(stdMat$m-1)){
            if(row == col) {
                specialHist(mat[,row],ann=F,axes=F,xlim=limits[,row])
                if(row == 1){
                    pos <- (plotW + plotW/2 + 2*par()$mai[2] + par()$mai[4]) / (plotW + sum(par()$mai[c(2,4)])) / stdMat$m
                    mtext(stdMat$names[1L],side=3,at=pos,outer=T,cex=par()$cex,las=1,adj=0.5)
                }
                next
            }
            if(col == 0){
                if(row==1){
                    plot.new()
                    next
                }
                specialHist(mat[,row],ann=F,axes=F,ylim=limits[,row],horiz=T)
                axis(side=2)
                irow <- stdMat$m - row + 1
                pos <- (plotH *(irow-1 + 1/2) + irow*par()$mai[2] + (irow-1)*par()$mai[4]) /
                        (plotH + sum(par()$mai[c(1,3)])) / stdMat$m
                mtext(stdMat$names[row],side=2,at=pos,outer=T,cex=par()$cex,adj=0.5,las=3,line=3)
                next
            }
            if(row == 1 & col==stdMat$m-2){
                plot.new()
#                text(0.5,0.9,"Credibility")
                legend("topright",legend=paste("<",names(CIPalette)[2:3],"%"),fill=CIPalette[1:2],bty="n",cex=par()$cex*3/2)
                pos <- (plotW * (stdMat$m-1) + (stdMat$m-0.5)*(par()$mai[2] + par()$mai[4])) /
                        (plotW + sum(par()$mai[c(1,3)])) / stdMat$m
                mtext("Credibility",side=3,at=pos,outer=T,cex=par()$cex*3/2,adj=0.5,las=1,line=0)
                next
            }
            if(row == 1 & col==stdMat$m-1){
                plot.new()
#                text(0.5,0.9,"Credibility")
                legend("topleft",legend=paste("<",names(CIPalette)[4:5],"%"),fill=CIPalette[3:4],bty="n",cex=par()$cex*3/2)
                next
            }
            if(col == row + 1)
            {
                plot.new()
                text(0.5,0.1,stdMat$names[col])
                next;
            }
            if(col > row){
                plot.new()
                next
            }
            axes <- (1:2)[c(row==stdMat$m,col==-1)]
            scatterPlotWiCred(subsetStdMat(stdMat,stdMat$names[c(col,row)]),axes,
                              col=pointCol,xlim=limits[,col],ylim=limits[,row])
            points(mode[col],mode[row],pch=19)
        }
    }
    par(mfrow=oldPar$mfrow,mar=oldPar$mar,oma=oldPar$oma,las=oldPar$las,xpd=oldPar$xpd,cex=oldPar$cex)
}

######## MAIN ##############


#stop("Here")

df <- read.table(resFile,sep="\t",stringsAsFactors=F,header=T,check.names=F)
df <- df[-1,]
nProt <- max(grep("Prot",names(df)))
n = nProt + sum(names(df)=="") + 1;
col.names = names(df)[-(1:n)];
modelName <- names(df)[nProt+1]
for (cn in col.names) {
    df[,cn] <- as.numeric(sapply(strsplit(df[,cn]," "),"[",2))
}

RowstoKeep = seq(1,nrow(df));
isFixed = sapply(setNames(nm=col.names[-1]),function(cn){length(unique(df[,cn])) == 1})

truncAway <- function(x,center=0){
    y <- ifelse(x < center, floor(x), ceiling(x))
    y[x == center] <- round(center,0)
    y
}

signLog <- function(x,base=exp(1)){
    sign(x) * log(abs(x)+1,base)
}

OoM = apply(df[,col.names[-1]][,!isFixed],2,function(x){diff(range(truncAway(signLog(x[!is.na(x)]))))})
ymin = apply(df[,col.names[-1]][,!isFixed],2,function(x){range(truncAway(signLog(x[!is.na(x)])))[1]})
#message(paste0(OoM,collapse=" "))
#message(paste0(ymin,collapse=" "))
#
#stop("There")

SwapIdx <- parseLog(logFile)

pdf(pdfFile,title=paste("ABC2 Results",resFile, sep= " - "))


burnin = kneedle(df$nLogP,mESSBurninEst(sum(!isFixed)),bPlot=TRUE)
RowstoKeep = RowstoKeep[RowstoKeep > burnin]

message("Calculating Multivariate Stats ...")
StdMat <- getStandardizedMatrix(df[RowstoKeep,col.names[-1][!isFixed]])
message("\tSmoothed Density ...")
StdMat <- getDensityEst(StdMat)
message("\tMode ...")
StdMat <- estimateMvMode(StdMat)
message("\tCredibility ...")
StdMat <- estimateCredibilityRadii(StdMat)
StdMat <- addFixedNames(StdMat,isFixed,unlist(df[1,names(isFixed)]))

pointCol <- ColourByCredibility(StdMat,CIPalette)
pointOrder <- order(StdMat$modalDist,decreasing=T)

message("Plotting ...")
SwapIdx = SwapIdx - burnin
SwapIdx = SwapIdx[SwapIdx > 0]
layout(matrix(c(rep(1,5),2),ncol=1))
mar <- par()$mar
for(cn in col.names){
    if(cn == col.names[1] | (cn %in% names(isFixed) & !isFixed[cn])){
        tmp <- df[-(1:burnin),cn]
        mar1 = c(0,mar[-1])
        mar1[3] = 6.1
        par(mar = mar1)
        x <- (1:length(tmp))[!is.na(tmp)]
        y <- tmp[!is.na(tmp)]
        plot(x,y,type="l",main=cn,ylab="",xaxt="n",xlab="")
        mar2 = mar; mar2[0] = 2.1; mar2[3]=0
        par(mar = mar2)
        plot(x,y,type="n",yaxt="n",ylab="swaps",xlab = "")
        abline(v=SwapIdx,col=rgb(0.5,0.5,0.5,0.5))
    }
}
layout(matrix(1,ncol=1))
par(mar = mar)
garbage <- lapply(split(col.names[-1][!isFixed],interaction(OoM,ymin,drop=T)),function(cn){vioplotWPoints(df[RowstoKeep[pointOrder],cn],names=cn,pointCol[pointOrder],getRescaledMode(StdMat)[cn])});

scatterPlotMatrix(StdMat);
garbage <- dev.off()

message("Outputting eval-prior ...")
lines <- paste0(">",modelName)
for(cn in col.names[-1]){
    mode <- df[1,cn];
    if(!isFixed[cn]){
        mode <- round(getRescaledMode(StdMat)[cn],7)
    }
    lines <- c(lines,paste0(cn,"\tFixed\tmu:",mode))
}
str <- unlist(df[1,1:nProt])
parts <- strsplit(str,"[;,:[\\]]",perl=T)
initLines <- lapply(parts,function(x){paste0("Init",x[c(3,5)],"_",x[2],"\tFixed\tmu:",x[c(4,6)])})
lines <- c(lines,sort(unlist(initLines)))

fileConn<-file(evalFile)
writeLines(lines,fileConn)
close(fileConn)

message("Outputting Credibility intervals ... ")
ci <- getRescaledCI(StdMat,CIalpha)
write.table(ci,ciFile,quote=F,sep=",",row.names=F,col.names=T)

message("Done")
