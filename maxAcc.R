##jmd
##may 2012
##maxAcc.R

require("lpSolve")
require("quadprog")

#########################################################################################
#### maxAUC
#########################################################################################
## construct Q matrix
makeQ <- function(n.bins){
    mat <- matrix(0,nrow=n.bins,ncol=n.bins)
    for (i in 1:nrow(mat)){
        for (j in 1:ncol(mat)){
            if (i>j){ mat[i,j] <- -(n.bins+i)*j/2 } else { mat[i,j] <- -(n.bins+j)*i/2 }
        }
    }
    return(mat)
}

## max AUC given k=prevalence, pve=proportion of variance explained (e.g. heritability)
maxAUC <- function(k, pve, n.bins=100){
    Q <- makeQ(n.bins=n.bins)
    #check concavity to ensure that Q is negative semidefinite
    eigen.vals <- eigen(x=Q, symmetric=TRUE, only.values=TRUE)$values
    if (!all(eigen.vals<=0)){ warning('The matrix Q is not concave') }
    #constraints: prevalence, pve, 0<=p, p<=1, sum(p)<=1
    Amat <- rbind((1:n.bins)/n.bins, (1:n.bins)^2/(n.bins^2), diag(n.bins), -diag(n.bins), rep(-1, n.bins))
    rhs <- c(k, k*(1-k)*pve+k^2, numeric(n.bins), rep(-1, n.bins), -1)
    #run quadratic program over p
    qp.res <- solve.QP(Dmat=-Q, dvec=numeric(n.bins), Amat=t(Amat), bvec=rhs, meq=2)
    auc <- (-2*qp.res$value+n.bins^2*k)/(n.bins^2*k*(1-k))
    return(auc)
}

#########################################################################################
#### optSe
#########################################################################################
## opt (Se|thresh)
optSeGivenThresh <- function(thresh, k, pve, sp, n.bins, direction="max"){
    obj <- c(rep(0, thresh), thresh:n.bins)/(n.bins*k)
    sp.coeff <- c(n.bins-0:(thresh-1), rep(0, n.bins-thresh+1))/(n.bins*(1-k))
    #constraints: avg risk=k, pve, sp, sum(p)=1, 0<=p, p<=1
    Amat <- rbind((0:n.bins)/n.bins, (0:n.bins)^2/(n.bins^2), sp.coeff, rep(1,n.bins+1), diag(n.bins+1), -diag(n.bins+1))
    rhs <- c(k, k*(1-k)*pve+k^2, sp, 1, numeric(n.bins+1), rep(-1,n.bins+1))
    lp.res <- lp(direction=direction, objective.in=obj, const.mat=Amat, const.dir=c(rep("=",4), 
    rep(">=", 2*n.bins+2)), const.rhs=rhs)
    return(lp.res$objval)
}

## opt Se over thresholds
optSe <- function(k, pve, sp, n.bins=1000, thresh.vector=seq(from=1,to=n.bins,by=10), direction="max"){
    #vector of sensitivity at each threshold
    se.tmp.v <- apply(X=as.matrix(thresh.vector), MARGIN=1, FUN=optSeGivenThresh, k=k, pve=pve, n.bins=n.bins,
    sp=sp, direction=direction)
    se <- max(se.tmp.v)
    return(se)
}

#########################################################################################
#### optSp
#########################################################################################
## opt(Sp|thresh)
optSpGivenThresh <- function(thresh, k, pve, se, n.bins, direction="max"){
    obj <- c(n.bins-0:(thresh-1), rep(0, n.bins-thresh+1))/(n.bins*(1-k))
    se.coeff <- c(rep(0, thresh), thresh:n.bins)/(n.bins*k)
    #constraints: avg risk=k, pve, sp, sum(p)=1, 0<=p, p<=1
    Amat <- rbind((0:n.bins)/n.bins, (0:n.bins)^2/(n.bins^2), se.coeff, rep(1,n.bins+1), diag(n.bins+1), -diag(n.bins+1))
    rhs <- c(k, k*(1-k)*pve+k^2, se, 1, numeric(n.bins+1), rep(-1,n.bins+1))
    lp.res <- lp(direction=direction, objective.in=obj, const.mat=Amat, const.dir=c(rep("=",4), 
    rep(">=", 2*n.bins+2)), const.rhs=rhs)
    return(lp.res$objval)
}

## opt Sp over thresholds
optSp <- function(k, pve, se, n.bins=1000, thresh.vector=seq(from=1,to=n.bins,by=10), direction="max"){
    #vector of specificity at each threshold
    sp.tmp.v <- apply(X=as.matrix(thresh.vector), MARGIN=1, FUN=optSpGivenThresh, k=k, pve=pve, n.bins=n.bins,
    se=se, direction=direction)
    sp <- max(sp.tmp.v)
    return(sp)
}
