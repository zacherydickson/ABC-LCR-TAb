args <- commandArgs(trailingOnly = T)

if(length(args) < 1){
    stop("Usage: Vis....R inFile.res (parses inFile.log also creates inFile.pdf and inFile.eval")
}

resFile <- args[1]
logFile <- sub("res$","log",resFile)
pdfFile <- sub("res$","pdf",resFile)
evalFile <- sub("res$","eval",resFile)
#burnin <- as.numeric(args[2])

if(any(resFile == c(logFile,pdfFile,evalFile))){
    stop("Could not generate unique log, pdf, and eval files from the inFile, does it end in 'res'?")
}

require("vioplot")
require("mvtnorm")


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

modeEst <- function(x){
   dens <- density(x[!is.na(x)])
   dens$x[which.max(dens$y)]
}

#point - vector defining cordiates in multivariate space
#mat - the matrix of values
#bw - the bandwidth of the density estimate
mvDensity <- function(point,mat,bw){
    p <- length(point)
    mu <- rep(0,p)
    sigma <- diag(bw,p,p)
    K_H <- function(x){
        dmvnorm(point - x,mu,sigma)
    }
    sum(apply(mat,1,K_H))
}

gradientEst <- function(point,dens,mat,bw,factor=10^-2){
    p <- length(point)
    jitPoint = point + apply(mat,2,function(x){min(abs(diff(sort(unique(x)))))}) * factor
    testPoints <- matrix(point,nrow=p,ncol=p,byrow=T)
    diag(testPoints) <- jitPoint
    testDens <- apply(testPoints,1,mvDensity,mat=mat,bw=bw)
    (testDens - dens) / diag(testPoints - point)
}

#Golden search assumes unimodal
goldenSearch <- function(point,dens,mat,bw,gradient,tolProp=10^-6){
    phi = (1 + sqrt(5)) / 2
    #Determine the extreme point; find the distance to the boundary in the gradient
    #direction for each dimension, the extreme is the point at the nearest boundary
    extrema <- mapply("[",asplit(apply(mat,2,range),2),(sign(gradient)+1)/2 + 1)
    #The maximum number of steps in the gradient Direction to take
    #The golden search will find the number of gradient steps to take to find a maximum
    #density
    limitDist <- min(abs(point - extrema) / abs(gradient))
    #The stopping point; when the smallest distance being considered is some fraction of
    #the total length
    tol = limitDist * tolProp
    #initialize probe points; f1(lower), f2(lower-mid), f3(upper), f4(upper-mid)
    pts <- list(f1=list(x=0,y=dens),f2 = list(x=limitDist / (1+phi),y=-1), f3 = list(x=limitDist,y=-1), f4 = list(x=-1,y=-1));
    mvDensityAtStep <- function(step){
        mvDensity(point + step * gradient,mat,bw)
    }
    #Reduce the search area until a lower-mid point is found which is greater than the
    #lower point; Esentially reducing the problem to finding the nearest local maximum to
    #the lower point
    repeat {
        pts$f2$y = mvDensityAtStep(pts$f2$x)
        if(pts$f2$y < pts$f1$y){
            pts$f3 = pts$f2
            pts$f2$x = (pts$f3$x - pts$f1$x) / (1+phi);
        } else {
            break;
        }
    }
    if(pts$f3$y < 0){
        pts$f3$y = mvDensityAtStep(pts$f3$x)
    }
    repeat {
        pts$f4$x = pts$f1$x + (pts$f3$x - pts$f2$x)
        pts$f4$y = mvDensityAtStep(pts$f4$x)
        #Set f4 always to be the righmost of f2 and f4
        if(pts$f4$x > pts$f2$x){
            tmp <- pts$f2
            pts$f2 = pts$f4
            pts$f4 = tmp
        }
        #Check stopping rule
        if(pts$f4$x - pts$f1$x < tol){
            break;
        }
        #If upper-mid is greater than lower-mid search lower-mid to upper
        if(pts$f4$y > pts$f2$y){
            pts$f1 = pts$f2
            pts$f2 = pts$f4
        } else { #search lower to upper-mid
            pts$f3 = pts$f4
        }
    }
    list(point = point + pts$f1$x * gradient, dens = pts$f1$y)
}

##Finds the multivariate Mode via a gradient ascent
##   Start at the point defined by the univariate modes
##   Until Sufficiently close to a peak
##       Evaluate the gradient at the current point
##       Use golden search to find the maximum density along the line segment from the
##         current point to the furthest in direction of the gradient
##       Set the current point to be the point at the maximum density
mvMode <- function(df,tol=10^-4,gssTolProp=10^-6){
    #Standardize the matrix to make density estimation easier
    stdInfo <- apply(df,2,function(x){c(mean(x),sd(x))})
    mat <- sapply(names(df),function(i){(df[,i]- stdInfo[1,i])/stdInfo[2,i]})
    #Take an initial guess of the mode at the point which is the mode for each individual
    #component, also get the bandwidth estimate
    densInfo <- apply(mat,2,function(x){dens <- density(x); c(dens$x[which.max(dens$y)],dens$bw)})
    mode <- densInfo[1,]
    bw <- mean(densInfo[2,])
    modalDens <- mvDensity(mode,mat,bw)
    res <- list(point = mode, dens = modalDens)
    repeat {
        gradient = gradientEst(res$point,res$dens,mat,bw)
        res <- goldenSearch(res$point,res$dens,mat,bw,gradient,gssTolProp)
        if(res$dens > modalDens){
            mode = res$point
            modalDens = res$dens
        } else {
            break
        }
        gradMag = sqrt(sum(gradient ^2))
        #Stopping Point: At a Peak the slope be zero
        if(gradMag < sqrt(length(mode)) * tol){
            break;
        }
    }
    #Return to the original scales
    mode * stdInfo[2,] + stdInfo[1,]
}

getScaledModalDist <- function(df,mode){
    stdInfo <- apply(df,2,function(x){c(mean(x),sd(x))})
    mat <- sapply(names(df),function(i){(df[,i]- stdInfo[1,i])/stdInfo[2,i]})
    mode <- (mode - stdInfo[1,]) / stdInfo[2,]
    apply(mat,1,function(x){sqrt(sum(abs(x - mode)^2))})
}

#stop("Here")

######## MAIN ##############

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

CIPalette = setNames(c("#FFFD98","#B3D2B2","#7A9CC6"),c("99","95","90"))


pdf(pdfFile,title=paste("ABC2 Results",resFile, sep= " - "))


#burnin = kneedle(df$nLogP,burnin,bPlot=TRUE,f=lowessFactor)
burnin = kneedle(df$nLogP,mESSBurninEst(sum(!isFixed)),bPlot=TRUE)
RowstoKeep = RowstoKeep[RowstoKeep > burnin]
message("Calculating Multivariate Mode ...")
multiVarMode <- mvMode(df[RowstoKeep,col.names[-c(1,which(isFixed)+1)]])
multiVarMode[names(isFixed)[isFixed]] = df[1,names(isFixed)[isFixed]]
multiVarMode <- multiVarMode[col.names[-1]]
multiVarMode = unlist(multiVarMode)
pointCol = rep("grey",length(RowstoKeep))
modalDist <- getScaledModalDist(df[RowstoKeep,col.names[-1]],multiVarMode)
pointCol[modalDist < quantile(modalDist,0.99)] = CIPalette["99"]
pointCol[modalDist < quantile(modalDist,0.95)] = CIPalette["95"]
pointCol[modalDist < quantile(modalDist,0.90)] = CIPalette["90"]
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
garbage <- lapply(split(col.names[-1][!isFixed],interaction(OoM,ymin,drop=T)),function(cn){vioplotWPoints(df[RowstoKeep,cn],names=cn,pointCol,multiVarMode[cn])});

garbage <- dev.off()


message("Outputting eval-prior ...")
#Mode Estimation

lines <- paste0(">",modelName)
for(cn in col.names[-1]){
    mode <- df[1,cn];
    if(!isFixed[cn]){
        mode <- round(multiVarMode[cn],7)
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
