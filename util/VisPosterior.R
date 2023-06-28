args <- commandArgs(trailingOnly = T)

if(length(args) < 2){
    stop("Usage: Vis....R inFile outFile [burnin]")
}

resFile <- args[1]
logFile <- sub("res$","log",resFile)
pdfFile <- args[2]
burnin <- as.numeric(args[3])


require("vioplot")

df <- read.table(resFile,sep="\t",stringsAsFactors=F,header=T,check.names=F)
df <- df[-1,]
nProt <- max(grep("Prot",names(df)))
n = nProt + sum(names(df)=="") + 1;
col.names = names(df)[-(1:n)];
modelName <- names(df)[nProt+1]
for (cn in col.names) {
    df[,cn] <- as.numeric(sapply(strsplit(df[,cn]," "),"[",2))
}

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

SwapIdx <- parseLog(logFile)


#On the assumption that The leftmost paramter is proposed first, and on each acceptance the
#next pramter is proposed
#If every column, except the first is shifted up by its paramter index, then all actual
#acceptances will be in the same row, and every nParamth row can be kept
#index = 1;
#for(i in 1:(n-1)){
#    cn = col.names[i+2]
#    df[,cn] = df[c((i+1):nrow(df),1:i),cn]
#}

#Effective Sample Size = n / (1 + 2*Σ_k cor@lag(k))
lagcor <- function(x,k){
    idx <- (1+k):length(x)
    cor(x[idx],x[idx-k])
}
#RowstoKeep = seq(1,nrow(df),by=n);
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
ymin = apply(df[,col.names[-1]][,!isFixed],2,function(x){y <- min(x[!is.na(x)]); round(log10(abs(y)),0)*sign(y)})


##Subsample to every seventh entry for each parameter, offset by when it was being adjusted
#df2 = df[1:(nrow(df)/n-n+1),];
#for (cn in col.names){
#    df2[,cn] <- df[seq(1,nrow(df),by=n)+match(cn,col.names)-2,cn]
#}
#
#df = df2;

#kneedle <- function(x,guess=length(x),...){
#    if(is.na(guess) || guess > length(x) / 2){
#        guess = length(x)/2
#    }
#    x = x[1:(guess*2)]
#    y = lowess(x,...)$y
#    y0 = y[1]
#    yn= y[length(y)]
#    x0=0
#    xn=length(y)-1
#    m = (yn-y0)/(xn-x0)
#    theta = atan(abs(1/m))
#    b = y0
#    d = abs(y - (m*(x0:xn) + b)) * sin(theta)
#    return(which.max(d) + 1)
#}

plotkneedle <- function(x,guess=length(x),...){
    if(is.na(guess) || guess > length(x) / 2){
        guess = length(x)/2
    }
    guidecol <- rgb(0.5,0.5,0.5,0.5)
    plot(x,main="Keedle Point Estimation",xlab="",ylab="",type="l")
    x = x[1:(guess*2)]
    tmp <- lowess(x,...)
    y = tmp$y
    lines(tmp$x,tmp$y,col=guidecol)
    segments(tmp$x[1],tmp$y[1],tmp$x[nrow(df)],tmp$y[nrow(df)],col=guidecol)
    y0 = y[1]
    yn= y[length(y)]
    x0=0
    xn=length(y)-1
    m = (yn-y0)/(xn-x0)
    theta = atan(abs(1/m))
    b = y0
    d = abs(y - (m*(x0:xn) + b)) * sin(theta)
    B <- which.max(d) + 1
    segments(tmp$x[B],tmp$y[B],tmp$x[B]+max(d)*cos(theta),tmp$y[B]+max(d)*sin(theta),col=guidecol)
    abline(v=B,lwd=3,col="red")
    return(B)
}


pdf(pdfFile,title="BC Posterior")

burnin = plotkneedle(df$nLogP,burnin)
RowstoKeep = RowstoKeep[RowstoKeep > burnin]
SwapIdx = SwapIdx - burnin
SwapIdx[SwapIdx > 0]
for(cn in col.names){
    if(cn == col.names[1] | (cn %in% names(isFixed) & !isFixed[cn])){
        tmp <- df[-(1:burnin),cn]
        plot((1:length(tmp))[!is.na(tmp)],tmp[!is.na(tmp)],type="l",main=cn,ylab="")
        abline(v=SwapIdx,col=rgb(0.5,0.5,0.5,0.5))
    }
}
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
cat(paste(lines,collapse="\n"))
cat("\n")

