require("truncnorm")

##PMF_Lt <- function(l,t,L0,lambda1,lambda2){
##    lambda = lambda1*(L0+1)+lambda2*L0
##    gamma = lambda1*(L0+1)
##    #Not complete
##}
##
##
###A0 and L0 are lists with two vector elements
###Values - value
###Probabilitiy Mass - prop
##PMF_At <- function(at,t,delta,sigma,tau,LR,AR,A0,L0){
##    AE <- function(l,alpha){
##        sel <- exp(-delta*tau)
##        alpha*sel + AR*(1-sel)*exp(tau * (l - LR))
##    }
##    PMF_Ae <- function(l,alpha){
##        p1 <- ptruncnorm(at+0.5,0,Inf,ae,sigma*t)
##        p3 <- 1-ptruncnorm(at-0.5,0,Inf,ae,sigma*t)
##        p5 <- sum(PMF_alpha(ae))
##        p <- p1 * p3 * p5
##        p
##    }
##    innerLoop <- function(alpha){
##        means <- AE(L0$value,alpha) 
##        p1 <- ptruncnorm(at+0.5,0,Inf,means,sigma*t)
##        p2 <- 1-ptruncnorm(at-0.5,0,Inf,means,sigma*t)
##        sum(p1 * p2 * L0$prob)
##    }
##    sum(sapply(A0$value,innerLoop)*A0$prob)
##}



simedge <- function(Len,Ab,time,lambda,kappa,delta,sigma,tau,A0,L0){
    outLen <- Len + rpois(1,(Len+1)*time*lambda) - rpois(1,Len*time*kappa)
    if(outLen < 0){
        outLen = 0;
    }
    selCoef <- exp(-delta*time)
    driftCoef <- sigma*time
    ae <- Ab * selCoef + A0 *(1-selCoef) * exp(tau * (Len - L0))
    outAb = round(rtruncnorm(1,0,Inf,ae,driftCoef),2)
    c(outLen,outAb)
}

sim <- function(tree,rootLen,rootAb,lambda,kappa,delta,sigma,tau){
    tree$lengths = c(rep(0,tree$Ntip),rootLen,rep(0,tree$Nnode-1))
    tree$abundances = c(rep(0,tree$Ntip),rootAb,rep(0,tree$Nnode-1))
    for(i in 1:nrow(tree$edge)){
        parentIdx <- tree$edge[i,1]
        vals <- simedge(tree$lengths[parentIdx],tree$abundances[parentIdx],tree$edge.length[i],lambda,kappa,delta,sigma,tau,rootAb,rootLen)
        childIdx <- tree$edge[i,2]
        tree$lengths[childIdx] = vals[1]
        tree$abundances[childIdx] = vals[2]
    }
    #matrix(c(tree$lengths[1:tree$Ntip],tree$abundances[1:tree$Ntip]),nrow=2,byrow=T,ncol =tree$Ntip,dimnames=list(c("Len","Ab"),tree$tip.label))
    matrix(c(tree$lengths,tree$abundances),nrow=2,byrow=T,ncol =length(tree$labels),dimnames=list(c("Len","Ab"),tree$labels))
}

runsim <- function(n,tree,rootLen,rootAb,lambda,kappa,delta,sigma,tau){
    tree$labels = c(tree$tip.label,tree$node.label)
    tree$Ntip = length(tree$tip.label)
    tree$tip.label = NULL
    tree$node.label = NULL
    res <- lapply(1:n,function(i){sim(tree,rootLen,rootAb,lambda,kappa,delta,sigma,tau)})
    res <- lapply(c("Len","Ab"),function(type){sapply(tree$labels,function(label){sapply(res,"[",type,label)})})
    setNames(res,c("Len","Ab"))
}

eval <- function(simres,obs){
    n <- nrow(simres$Len)
    logP = 0
    for(i in 1:nrow(obs)){
        label <- obs$label[i]
        Ab <- obs$Ab[i]
        Len <- obs$Len[i]
        p = (sum(simres$Ab[,label] == Ab)+1) * (sum(simres$Len[,label] == Len)+1) / (n+2)^2
        logP = logP + log(p)
    } 
    logP
}

kneedle <- function(x,guess=length(x),smooth=TRUE){
    if(is.na(guess) || guess > length(x) / 2){
        guess = length(x)/2
    }
    x = x[1:(guess*2)]
    y <- x;
    if(smooth){
        y <- lowess(x)$y
    }
    y0 = y[1]
    yn= y[length(y)]
    x0=0
    xn=length(y)-1
    m = (yn-y0)/(xn-x0)
    b = y0
    d = abs(y - (m*(x0:xn) + b)) * asin(pi/4)
    return(which.max(d) + 1)
}
