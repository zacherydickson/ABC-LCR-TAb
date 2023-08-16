args <- commandArgs(trailingOnly = T)

if(length(args) < 1){
    stop("Usage: Vis....R inFile.res (parses inFile.log also creates inFile.pdf and inFile.eval")
}

resFile <- args[1]
#threads <- as.numeric(args[2])
#if(is.na(threads)){
#    threads <- 1
#}
logFile <- sub("res$","log",resFile)
pdfFile <- sub("res$","pdf",resFile)
evalFile <- sub("res$","eval",resFile)
#burnin <- as.numeric(args[2])

if(any(resFile == c(logFile,pdfFile,evalFile))){
    stop("Could not generate unique log, pdf, and eval files from the inFile, does it end in 'res'?")
}

require("vioplot")
require("Rcpp")
CPPSource <- "/home/zac/scripts/LCR/TemporalOrder/Primates/ABC2/util/VisPosterior.cpp"
sourceCpp(CPPSource)


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

vioplotWPoints <- function(data,pointCol=NULL,mode=NULL,...){
    vioplot(data,...,side="left",plotCentre="line")
    if(is.null(pointCol)){
        pointCol = rep("black",nrow(data))
    }
    if(is.null(ncol(data))){ #input is a vector
        dJit <- 1.05 + densityJitter(data,0,0.4)
        points(dJit,data,col=pointCol)
        if(!is.null(mode)){
            segments(1,mode,1.5,mode,col="black",lwd=2)
        }
        return(NULL)
    }
    dJit <- lapply(data,densityJitter,b=0.4)
    dJit <- mapply("+",dJit,1:ncol(data)+0.05,SIMPLIFY=F)
    mapply(function(x,y){points(x,y,col=pointCol)},dJit,data)
    if(!is.null(mode)){
        segments(1:length(data),mode,1:length(data)+ 0.5,mode,col="black",lwd=2)
    }

    return(NULL)
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

K_H <- function(p1,p2,mu,sigma){
    dmvnorm(p1-p2,mu,sigma)
}


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
        if(sqrt(sum(gradient^2)) < tol){
            break
        }
        r <- findUpperMult(point,val,jitSize,gradient,calculateMVDensity,data=mat,bw=bw) 
        if(r <= 1){
            break
        }
        res <- generalizedGoldenSearch(1,r,function(mult,...){p <- point + mult*jitSize*gradient; calculateMVDensity(point,...)},data=mat,bw=bw)
        if(res[2] < val){
            break
        }
        val <- res[2]
        point <- point + res[1]*jitSize*gradient
    }
    list(point=point,val=val)
}

##Golden search assumes unimodal
#goldenSearch <- function(point,dens,mat,bw,gradient,tolProp=10^-6){
#    phi = (1 + sqrt(5)) / 2
#    #Determine the extreme point; find the distance to the boundary in the gradient
#    #direction for each dimension, the extreme is the point at the nearest boundary
#    extrema <- mapply("[",asplit(apply(mat,2,range),2),(sign(gradient)+1)/2 + 1)
#    #The maximum number of steps in the gradient Direction to take
#    #The golden search will find the number of gradient steps to take to find a maximum
#    #density
#    limitDist <- min(abs(point - extrema) / abs(gradient))
#    #The stopping point; when the smallest distance being considered is some fraction of
#    #the total length
#    tol = limitDist * tolProp
#    #initialize probe points; f1(lower), f2(lower-mid), f3(upper), f4(upper-mid)
#    pts <- list(f1=list(x=0,y=dens),f2 = list(x=limitDist / (1+phi),y=-1), f3 = list(x=limitDist,y=-1), f4 = list(x=-1,y=-1));
#    mvDensityAtStep <- function(step){
#        calculateMVDensity(point + step * gradient,mat,bw)
#    }
#    #Reduce the search area until a lower-mid point is found which is greater than the
#    #lower point; Esentially reducing the problem to finding the nearest local maximum to
#    #the lower point
#    repeat {
#        pts$f2$y = mvDensityAtStep(pts$f2$x)
#        if(pts$f2$y < pts$f1$y){
#            pts$f3 = pts$f2
#            pts$f2$x = (pts$f3$x - pts$f1$x) / (1+phi);
#        } else {
#            break;
#        }
#    }
#    if(pts$f3$y < 0){
#        pts$f3$y = mvDensityAtStep(pts$f3$x)
#    }
#    repeat {
#        pts$f4$x = pts$f1$x + (pts$f3$x - pts$f2$x)
#        pts$f4$y = mvDensityAtStep(pts$f4$x)
#        #Set f4 always to be the righmost of f2 and f4
#        if(pts$f4$x > pts$f2$x){
#            tmp <- pts$f2
#            pts$f2 = pts$f4
#            pts$f4 = tmp
#        }
#        #Check stopping rule
#        if(pts$f4$x - pts$f1$x < tol){
#            break;
#        }
#        #If upper-mid is greater than lower-mid search lower-mid to upper
#        if(pts$f4$y > pts$f2$y){
#            pts$f1 = pts$f2
#            pts$f2 = pts$f4
#        } else { #search lower to upper-mid
#            pts$f3 = pts$f4
#        }
#    }
#    list(point = point + pts$f1$x * gradient, dens = pts$f1$y)
#}
#
#sampleDensity <- function(x,n){
#    d <- density(x)
#    sample(d$x,n,replace=T,p=d$y)
#}

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
    #Add Fields to be added
    obj$mode <- numeric(0)
    obj$modalDensity <- numeric(0)
    obj$density <- numeric(0)
    obj$modalDist <- numeric(0)
    obj$credRadii <- numeric(0)
    obj
}


###Finds the multivariate Mode via a gradient ascent
###   Start at the point defined by the univariate modes
###   Until Sufficiently close to a peak
###       Evaluate the gradient at the current point
###       Use golden search to find the maximum density along the line segment from the
###         current point to the furthest in direction of the gradient
###       Set the current point to be the point at the maximum density
#estimateMvMode <- function(stdMat,tol=10^-4,gssTolProp=10^-6){
#    if(attr(stdMat,"class") != "StandardizedMatrix"){
#        stop("Attempt to estimate multivariate mode with non-StadardizedMatrix object")
#    }
#    #Standardize the matrix to make density estimation easier
#    mode <- stdMat$uniMode
#    bw <- stdMat$bw
#    modalDens <- mvDensity(mode,stdMat$mat,bw)
#    res <- list(point = mode, dens = modalDens)
#    repeat {
#        gradient <- gradientEst(res$point,res$dens,stdMat$mat,bw)
#        res <- goldenSearch(res$point,res$dens,stdMat$mat,bw,gradient,gssTolProp)
#        if(res$dens > modalDens){
#            mode = res$point
#            modalDens = res$dens
#        } else {
#            break
#        }
#        gradMag = sqrt(sum(gradient ^2))
#        #Stopping Point: At a Peak the slope be zero
#        if(gradMag < sqrt(length(mode)) * tol){
#            break;
#        }
#    }
#    stdMat$mode <- mode
#    #Calculate the euclidian distance for each point from the mode for every point
#    stdMat$modalDist <- apply(stdMat$mat,1,function(x){sqrt(sum((x - mode)^2))})
#    stdMat
#    #Return to the original scales
#    #mode * stdMat$stdInfo[2,] + stdMat$stdInfo[1,]
#}


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

addFixedNames <- function(stdMat,isFixed,mu){
    if(attr(stdMat,"class") != "StandardizedMatrix"){
        stop("Attempt to addFixedNames from non-StandardizedMatrix object")
    }
    stdMat$names <- names(isFixed)
    #Add info as necessary
    fixedMu <- mu[isFixed]
    stdMat$uniMode <- c(stdMat$uniMode,fixedMu)
    stdMat$stdInfo <- cbind(stdMat$stdInfo,sapply(fixedMu,c,1))
    stdMat$mat <- cbind(stdMat$mat,sapply(fixedMu,function(i){rep(0,stdMat$n)}))
    stdMat$mode <- c(stdMat$mode,fixedMu)
    #Reorder
    for(member in c("uniMode","mode")){
        stdMat[[member]] <- stdMat[[member]][stdMat$names]
    }
    for(member in c("stdInfo","mat")){
        stdMat[[member]] <- stdMat[[member]][,stdMat$names]
    }
    stdMat
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

#EstimateMVDensity <- function(stdMat,nProp=0.001){
#    if(attr(stdMat,"class") != "StandardizedMatrix"){
#        stop("Attempt to EstimateMVDensity with a non-StandardizedMatrix object")
#    }
#    if(length(stdMat$modalDist) == 0){
#        stop("Attempt to EstimateMVDensity with improperly initialized StandardizedMatrix")
#    }
#    n <- stdMat$n
#    modalDist <- stdMat$modalDist
#    #Process densities in order fom close to the mode to far
#    densOrder <- order(modalDist)
#    #Take a sample of all the densities
#    nPoints = max(2,as.integer(nProp*n))
#    sampleIdx <- as.integer(seq(1,n,l=nPoints))
#    densitySample <- apply(stdMat$mat[densOrder[sampleIdx],],1,calculateMVDensity,mat=stdMat$mat,bw=stdMat$bw)
#    #linearly interpolate densities between pairs of points
#    density <- densitySample[1]
#    for(i in 2:nPoints){
#        slope <- (densitySample[i] - densitySample[i-1]) / (sampleIdx[i] - sampleIdx[i-1])
#        density <- c(density,1:(diff(sampleIdx[(i-1):i]) * slope + densitySample[i-1]))
#    }
#    #Put densities back into the order of the points
#    stdMat$density <- density[match(1:n,densOrder)]
#    stdMat
#}

estimateCredibilityRadii <- function(stdMat,cutoffs = c(0.9,0.95,0.99)){
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

scatterPlotMatrix <- function(stdMat,...){

}

######## MAIN ##############


stop("Here")

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
isFixed = setNames(rep(FALSE,length(col.names[-1])),col.names[-1])

#for(cn in col.names[-1]){
#    tmp <- rle(df[,cn])
#    if(length(tmp$values) == 1){
#        isFixed[cn]=TRUE
#    }
#    tmp <- lapply(1:length(tmp$values),function(i){c(tmp$values[i],rep(NA,tmp$lengths[i]-1))})
#    df[,cn] = unlist(tmp)
#}

isFixed = sapply(names(isFixed),function(cn){length(unique(df[,cn])) == 1})

OoM = apply(df[,col.names[-1]][,!isFixed],2,function(x){round(log10(diff(range(x[!is.na(x)]))),0)})
ymin = apply(df[,col.names[-1]][,!isFixed],2,function(x){y <- median(x[!is.na(x)]); ceiling(log10(abs(y)))*sign(y)})
#message(paste0(OoM,collapse=" "))
#message(paste0(ymin,collapse=" "))

SwapIdx <- parseLog(logFile)
#lowessFactor <- kneedle(abs(sapply(1:(2*nrow(df)/3),function(l){tmp <- windowedSumDeviation(df$nLogP,l); sum(tmp)*length(tmp)})))/nrow(df)

CIPalette = setNames(c("grey","#7A9CC6","#B3D2B2","#FFFD98"),c("0","90","95","99"))


pdf(pdfFile,title=paste("ABC2 Results",resFile, sep= " - "))


#burnin = kneedle(df$nLogP,burnin,bPlot=TRUE,f=lowessFactor)
burnin = kneedle(df$nLogP,mESSBurninEst(sum(!isFixed)),bPlot=TRUE)
RowstoKeep = RowstoKeep[RowstoKeep > burnin]

message("Calculating Multivariate Stats ...")
StdMat <- getStandardizedMatrix(df[RowstoKeep,col.names[-1][!isFixed]])
StdMat <- getDensityEst(StdMat)
StdMat <- estimateMvMode(StdMat)
StdMat <- estimateCredibilityRadii(StdMat)

pointCol <- ColourByCredibility(StdMat,CIPalette)
#pointCol <- ColourByModalDistance(StdMat,CIPalette)
#pointCol = NULL


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
garbage <- lapply(split(col.names[-1][!isFixed],interaction(OoM,ymin,drop=T)),function(cn){vioplotWPoints(df[RowstoKeep,cn],names=cn,pointCol,getRescaledMode(StdMat)[cn])});
scatterPlotMatrix(StdMat);

garbage <- dev.off()


message("Outputting eval-prior ...")

StdMat <- addFixedNames(StdMat,isFixed,unlist(df[1,names(isFixed)]))

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

message("Done")
