args <- commandArgs(trailingOnly = T)

if(length(args) < 2){
    stop("Usage: Vis....R inFile outFile [burnin]")
}

resFile <- args[1]
pdfFile <- args[2]
burnin <- as.numeric(args[3])


require("vioplot")

df <- read.table(resFile,sep="\t",stringsAsFactors=F,header=T,check.names=F)
nProt <- max(grep("Prot",names(df)))
n = nProt + sum(names(df)=="") + 1;
col.names = names(df)[-(1:n)];
modelName <- names(df)[nProt+1]
for (cn in col.names) {
    df[,cn] <- as.numeric(sapply(strsplit(df[,cn]," "),"[",2))
}

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
ymin = apply(df[,col.names[-1]][,!isFixed],2,function(x){sign(min(x[!is.na(x)]))})



##Subsample to every seventh entry for each parameter, offset by when it was being adjusted
#df2 = df[1:(nrow(df)/n-n+1),];
#for (cn in col.names){
#    df2[,cn] <- df[seq(1,nrow(df),by=n)+match(cn,col.names)-2,cn]
#}
#
#df = df2;

kneedle <- function(x,guess=length(x)){
    if(is.na(guess) || guess > length(x) / 2){
        guess = length(x)/2
    }
    x = x[1:(guess*2)]
    y = lowess(x)$y
    y0 = y[1]
    yn= y[length(y)]
    x0=0
    xn=length(y)-1
    m = (yn-y0)/(xn-x0)
    b = y0
    d = abs(y - (m*(x0:xn) + b)) * asin(pi/4)
    return(which.max(d) + 1)
}

burnin = kneedle(df$nLogP,burnin)
RowstoKeep = RowstoKeep[RowstoKeep > burnin]

pdf(pdfFile,title="BC Posterior")

plot(df$nLogP,type="l")
abline(v=burnin,lwd=3,col="red")
for(cn in col.names){
    if(cn == col.names[1] | (cn %in% names(isFixed) & !isFixed[cn])){
        tmp <- df[-(1:burnin),cn]
        plot((1:length(tmp))[!is.na(tmp)],tmp[!is.na(tmp)],type="l",main=cn,ylab="")
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

