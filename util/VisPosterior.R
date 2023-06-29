args <- commandArgs(trailingOnly = T)

if(length(args) < 1){
    stop("Usage: Vis....R inFile.res [burnin] (parses inFile.log also creates inFile.pdf and inFile.eval")
}

resFile <- args[1]
logFile <- sub("res$","log",resFile)
pdfFile <- sub("res$","pdf",resFile)
evalFile <- sub("res$","eval",resFile)
burnin <- as.numeric(args[2])

if(any(resFile == c(logFile,pdfFile,evalFile))){
    stop("Could not generate unique log, pdf, and eval files from the inFile, does it end in 'res'?")
}

require("vioplot")


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

for(cn in col.names[-1]){
    tmp <- rle(df[,cn])
    if(length(tmp$values) == 1){
        isFixed[cn]=TRUE
    }
    tmp <- lapply(1:length(tmp$values),function(i){c(tmp$values[i],rep(NA,tmp$lengths[i]-1))})
    df[,cn] = unlist(tmp)
}


OoM = apply(df[,col.names[-1]][,!isFixed],2,function(x){round(log10(diff(range(x[!is.na(x)]))),0)})
ymin = apply(df[,col.names[-1]][,!isFixed],2,function(x){y <- min(x[!is.na(x)]); ceiling(log10(abs(y)))*sign(y)})
message(paste0(OoM,collapse=" "))
message(paste0(ymin,collapse=" "))

SwapIdx <- parseLog(logFile)
lowessFactor <- kneedle(abs(sapply(1:(2*nrow(df)/3),function(l){tmp <- windowedSumDeviation(df$nLogP,l); sum(tmp)*length(tmp)})))/nrow(df)

pdf(pdfFile,title=paste("ABC2 Results",resFile, sep= " - "))

burnin = kneedle(df$nLogP,burnin,bPlot=TRUE,f=lowessFactor)
RowstoKeep = RowstoKeep[RowstoKeep > burnin]
SwapIdx = SwapIdx - burnin
SwapIdx = SwapIdx[SwapIdx > 0]
layout(matrix(c(rep(1,5),2),ncol=1))
for(cn in col.names){
    if(cn == col.names[1] | (cn %in% names(isFixed) & !isFixed[cn])){
        tmp <- df[-(1:burnin),cn]
        mar <- par()$mar
        mar1 = c(0,mar[-1])
        par(mar = mar1)
        x <- (1:length(tmp))[!is.na(tmp)]
        y <- tmp[!is.na(tmp)]
        plot(x,y,type="l",main=cn,ylab="",xaxt="n",xlab="")
        mar2 = mar; mar2[3]=0
        par(mar = mar2)
        plot(x,y,type="n",yaxt="n",ylab="swaps",xlab = "")
        abline(v=SwapIdx,col=rgb(0.5,0.5,0.5,0.5))
    }
}
layout(matrix(1,ncol=1))
garbage <- lapply(split(col.names[-1][!isFixed],interaction(OoM,ymin,drop=T)),function(cn){vioplot(df[RowstoKeep,cn],names=cn)});

garbage <- dev.off()


#Mode Estimation

lines <- paste0(">",modelName)
for(cn in col.names[-1]){
    mode <- df[1,cn];
    if(!isFixed[cn]){
        tmp <- df[RowstoKeep,cn]
        dens <- density(tmp[!is.na(tmp)])
        mode <- dens$x[which.max(dens$y)]
        mode <- round(mode,7)
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
