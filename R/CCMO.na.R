CCMO.na <- function(Y,gm,gc,Xo,Xm,Xc,f,ind){
  fo = "~"
  nXo <- ifelse(is.null(Xo),0,max(1,ncol(Xo)))
  if(!is.null(Xo)) fo = paste0(fo,"+Xo")
  nXm <- ifelse(is.null(Xm),0,max(1,ncol(Xm)))
  if(!is.null(Xm)) fo = paste0(fo,"+Xm")
  nXc <- ifelse(is.null(Xc),0,max(1,ncol(Xc)))
  if(!is.null(Xc)) fo = paste0(fo,"+Xc")
  nbeta <- 3 + nXo + nXm + nXc
  fo = paste(fo,"-1")
  X <- model.matrix(formula(fo))
  nX <- c(nXo,nXm,nXc)
  
  keep <- !is.na(gm)
  Xt <- as.matrix(X[keep,])
  Yt <- Y[keep]
  gmt <- gm[keep]
  gct <- gc[keep]
  
  n1 <- sum(Yt == 1)
  n0 <- sum(Yt == 0)
  n <- n1+n0
  lambda <- n1 / (n*f) - n0 / (n * (1-f))
  
  group2 <- is.na(gct)
  X1 <- as.matrix(Xt[!group2,])
  Y1 <- Yt[!group2]
  gm1 <- gmt[!group2]
  gc1 <- gct[!group2]
  X2 <- as.matrix(Xt[group2,])
  Y2 <- Yt[group2]
  gm2 <- gmt[group2]
  gc2 <- gct[group2]
  
  Z <- cbind(1,gm1,gc1,X1[,1:nXo],gm1*X1[,(nXo+1):(nXo+nXm)],gc1*X1[,(nXo+nXm+1):(nXo+nXm+nXc)])
  fit <- glm(Y1 ~ 0 + Z,family = binomial)
  res <- summary(fit)$coef
  est.log <- as.vector(res[,1])
  sd.log <- as.vector(res[,2])
  
  beta0 <- est.log
  theta0 <- (2*sum(gm1 == 2)+sum(gm1 == 1))/(2*length(gm1))
  theta0 <- log(theta0/(1-theta0))
  para0 <- c(beta0,theta0)
  
  if(ind == TRUE){
    llik <- function(para)
      -likeli.ccmo.ind(para,Y1,X1,gm1,gc1,f,lambda,n,nX) - likeli.ccmo.ind.na(para,Y2,X2,gm2,gc2,f,lambda,n,nX)
    fit <- optim(par = para0,fn = llik,method = 'L-BFGS-B',hessian = TRUE)
    est <- fit$par #beta 6, theta 1
    est <- c(est[1:nbeta],exp(est[nbeta+1])/(1 + exp(est[nbeta+1])))
    Matv <- solve(fit$hessian)[1:nbeta,1:nbeta]
    sd <- sqrt(diag(Matv))
  }
  
  if(ind == FALSE){
    llik <- function(para)
      -likeli.ccmo.rob(para,Y1,X1,gm1,gc1,f,lambda,n,nX) - likeli.ccmo.rob.na(para,Y2,X2,gm2,gc2,f,lambda,n,nX)
    fit <- optim(par = para0,fn = llik,method = 'L-BFGS-B',hessian = TRUE)
    est <- fit$par #beta 6, theta 1
    est <- c(est[1:nbeta],exp(est[nbeta+1])/(1 + exp(est[nbeta+1])))
    Matv <- solve(fit$hessian)[1:nbeta,1:nbeta]
    sd <- sqrt(diag(Matv))
  }
  return(list(est = est,sd = sd,Matv = Matv,est.log = est.log,sd.log = sd.log,logL = -fit$value))
}

likeli.ccmo.rob <- function(para,Y,X,gm,gc,f,lambda,n,nX){
  res <- 0
  nbeta <- 3 + nX[1] + nX[2] + nX[3]
  beta <- para[1:nbeta]
  theta <- exp(para[nbeta+1]) / (1 + exp(para[nbeta+1]))
  for(i in 1:length(gm)){
    res <- res + log(P_Y.ccmo(Y[i],X[i,],gm[i],gc[i],beta,nX)) + log(P_gc.ccmo(gm[i],gc[i],theta))
    res <- res - log(n * (1 + lambda * (H.ccmo.rob(gm[i],X[i,],beta,theta,nX) - f)))
  }
  return(res)
}

likeli.ccmo.ind <- function(para,Y,X,gm,gc,f,lambda,n,nX){
  res <- 0
  nbeta <- 3 + nX[1] + nX[2] + nX[3]
  beta <- para[1:nbeta]
  theta <- exp(para[nbeta+1]) / (1 + exp(para[nbeta+1]))
  for(i in 1:length(gm)){
    res <- res + log(P_Y.ccmo(Y[i],X[i,],gm[i],gc[i],beta,nX)) + log(P_gc.ccmo(gm[i],gc[i],theta)) + log(P_gm.ccmo.ind(gm[i],theta))
    res <- res - log(n * (1 + lambda * (H.ccmo.ind(X[i,],beta,theta,nX) - f)))
  }
  return(res)
}

P_Y.ccmo <- function(Y,X,gm,gc,beta,nX){
  Xo <- X[1:nX[1]]
  Xm <- X[(nX[1]+1):(nX[1]+nX[2])]
  Xc <- X[(nX[1]+nX[2]+1):(nX[1]+nX[2]+nX[3])]
  tmp <- exp(sum(c(1,gm,gc,Xo,gm*Xm,gc*Xc) * beta))
  return(1-Y+(2*Y-1)*tmp/(1+tmp))
}

P_gc.ccmo <- function(gm,gc,theta){
  Pg <- matrix(c(1-theta,theta,0,(1-theta)/2,0.5,theta/2,0,1-theta,theta),3,3)
  return(Pg[gc+1,gm+1])
}

P_gm.ccmo.ind <- function(gm,theta){
  Pgm <- c((1-theta)^2,2*theta*(1-theta),theta^2)
  return(Pgm[gm+1])
}

H.ccmo.ind <- function(X,beta,theta,nX){
  res <- P_Y.ccmo(1,X,0,0,beta,nX) * P_gc.ccmo(0,0,theta) * P_gm.ccmo.ind(0,theta)
  res <- res + P_Y.ccmo(1,X,0,1,beta,nX) * P_gc.ccmo(0,1,theta) * P_gm.ccmo.ind(0,theta)
  res <- res + P_Y.ccmo(1,X,1,0,beta,nX) * P_gc.ccmo(1,0,theta) * P_gm.ccmo.ind(1,theta)
  res <- res + P_Y.ccmo(1,X,1,1,beta,nX) * P_gc.ccmo(1,1,theta) * P_gm.ccmo.ind(1,theta)
  res <- res + P_Y.ccmo(1,X,1,2,beta,nX) * P_gc.ccmo(1,2,theta) * P_gm.ccmo.ind(1,theta)
  res <- res + P_Y.ccmo(1,X,2,1,beta,nX) * P_gc.ccmo(2,1,theta) * P_gm.ccmo.ind(2,theta)
  res <- res + P_Y.ccmo(1,X,2,2,beta,nX) * P_gc.ccmo(2,2,theta) * P_gm.ccmo.ind(2,theta)
  return(res)
}

H.ccmo.rob <- function(gm,X,beta,theta,nX){
  res <- P_Y.ccmo(1,X,gm,0,beta,nX) * P_gc.ccmo(gm,0,theta)
  res <- res + P_Y.ccmo(1,X,gm,1,beta,nX) * P_gc.ccmo(gm,1,theta)
  res <- res + P_Y.ccmo(1,X,gm,2,beta,nX) * P_gc.ccmo(gm,2,theta)
  return(res)
}

likeli.ccmo.rob.na <- function(para,Y,X,gm,gc,f,lambda,n,nX){
  res <- 0
  nbeta <- 3 + nX[1] + nX[2] + nX[3]
  beta <- para[1:nbeta]
  theta <- exp(para[nbeta+1]) / (1 + exp(para[nbeta+1]))
  for(i in 1:length(gm)){
    res <- res + log(likeli0.ccmo.rob.na(Y[i],X[i,],gm[i],beta,theta,nX))
    res <- res - log(n * (1 + lambda * (H.ccmo.rob(gm[i],X[i,],beta,theta,nX) - f)))
  }
  return(res)
}

likeli0.ccmo.rob.na <- function(Y,X,gm,beta,theta,nX){
  res <- 0
  for(gc in 0:2)
    res <- res + P_Y.ccmo(Y,X,gm,gc,beta,nX) * P_gc.ccmo(gm,gc,theta)
  return(res)
}

likeli.ccmo.ind.na <- function(para,Y,X,gm,gc,f,lambda,n,nX){
  res <- 0
  nbeta <- 3 + nX[1] + nX[2] + nX[3]
  beta <- para[1:nbeta]
  theta <- exp(para[nbeta+1]) / (1 + exp(para[nbeta+1]))
  for(i in 1:length(gm)){
    res <- res + log(likeli0.ccmo.ind.na(Y[i],X[i,],gm[i],beta,theta,nX))
    res <- res - log(n * (1 + lambda * (H.ccmo.ind(X[i,],beta,theta,nX) - f)))
  }
  return(res)
}

likeli0.ccmo.ind.na <- function(Y,X,gm,beta,theta,nX){
  res <- 0
  for(gc in 0:2)
    res <- res + P_Y.ccmo(Y,X,gm,gc,beta,nX) * P_gc.ccmo(gm,gc,theta) * P_gm.ccmo.ind(gm,theta)
  return(res)
}
