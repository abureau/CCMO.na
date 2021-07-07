CCMO.na <- function(Y,gm,gc,Xo,Xm,Xc,Xgm,f){
  n1 <- sum(Y == 1)
  n0 <- sum(Y == 0)
  n <- n1+n0
  lambda <- n1 / (n*f) - n0 / (n * (1-f))
  fo = "~"
  nXo <- ifelse(is.null(Xo),0,max(1,ncol(Xo)))
  if(!is.null(Xo)) fo = paste0(fo,"+Xo")
  nXm <- ifelse(is.null(Xm),0,max(1,ncol(Xm)))
  if(!is.null(Xm)) fo = paste0(fo,"+Xm")
  nXc <- ifelse(is.null(Xc),0,max(1,ncol(Xc)))
  if(!is.null(Xc)) fo = paste0(fo,"+Xc")
  nXgm <- ifelse(is.null(Xgm),0,max(1,ncol(Xgm)))
  if(!is.null(Xgm)) fo = paste0(fo,"+Xgm")
  nbeta <- 3 + nXo + nXm + nXc
  fo = paste(fo,"-1")
  X <- model.matrix(formula(fo))
  nX <- c(nXo,nXm,nXc,nXgm)
  
  group2 <- apply(is.na(cbind(gm,gc)),1,any)
  X1 <- as.matrix(X[!group2,])
  Y1 <- Y[!group2]
  gm1 <- gm[!group2]
  gc1 <- gc[!group2]
  X2 <- as.matrix(X[group2,])
  Y2 <- Y[group2]
  gm2 <- gm[group2]
  gc2 <- gc[group2]
  
  Z <- cbind(1,gm1,gc1,X1[,1:nXo],gm1*X1[,(nXo+1):(nXo+nXm)],gc1*X1[,(nXo+nXm+1):(nXo+nXm+nXc)])
  fit <- glm(Y1 ~ 0 + Z,family = binomial)
  res <- summary(fit)$coef
  est.log <- as.vector(res[,1])
  sd.log <- as.vector(res[,2])
  
  beta0 <- est.log
  theta0 <- (2*sum(gm1 == 2)+sum(gm1 == 1))/(2*length(gm1))
  theta0 <- log(theta0/(1-theta0))
  eta0 <- rep(0,nXgm)
  para0 <- c(beta0,theta0,eta0)
  
  llik <- function(para)
    -likeli.ccmo(para,Y1,X1,gm1,gc1,f,lambda,n,nX) - likeli.ccmo.na(para,Y2,X2,gm2,gc2,f,lambda,n,nX)
  fit <- optim(par = para0,fn = llik,method = 'L-BFGS-B',hessian = TRUE)
  est <- fit$par #beta 6, theta 1, eta 1
  if(nXgm>0) est <- c(est[1:nbeta],exp(est[nbeta+1])/(1 + exp(est[nbeta+1])),est[(nbeta+2):(nbeta+1+nXgm)])
  else est <- c(est[1:nbeta],exp(est[nbeta+1])/(1 + exp(est[nbeta+1])))
  Matv <- solve(fit$hessian)[1:nbeta,1:nbeta]
  sd <- sqrt(diag(Matv))
  return(list(est = est,sd = sd,Matv=Matv,est.log = est.log,sd.log = sd.log,logL = -fit$value))
}

likeli.ccmo <- function(para,Y,X,gm,gc,f,lambda,n,nX){
  res <- 0
  nbeta <- 3 + nX[1] + nX[2] + nX[3]
  beta <- para[1:nbeta]
  theta <- exp(para[nbeta+1]) / (1 + exp(para[nbeta+1]))
  eta <- para[(nbeta+2):(nbeta+1+nX[4])]
  Xgm <- as.matrix(X[,-(1:(nX[1]+nX[2]+nX[3]))])
  for(i in 1:length(gm)){
    res <- res + log(P_Y.ccmo(Y[i],X[i,],gm[i],gc[i],beta,nX)) + log(P_gc.ccmo(gm[i],gc[i],theta)) + log(P_gm.ccmo(gm[i],Xgm[i,],theta,eta))
    res <- res - log(n * (1 + lambda * (H.ccmo(X[i,],beta,theta,eta,nX) - f)))
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

P_gm.ccmo <- function(gm,Xgm,theta,eta){
  Pgm <- c((1-theta)^2,2*theta*(1-theta)*exp(sum(eta*Xgm)),theta*theta*exp(2*sum(eta*Xgm)))
  Pgm <- Pgm/sum(Pgm)
  return(Pgm[gm+1])
}

H.ccmo <- function(X,beta,theta,eta,nX){
  Xgm <- X[-(1:(nX[1]+nX[2]+nX[3]))]
  res <- P_Y.ccmo(1,X,0,0,beta,nX) * P_gc.ccmo(0,0,theta) * P_gm.ccmo(0,Xgm,theta,eta)
  res <- res + P_Y.ccmo(1,X,0,1,beta,nX) * P_gc.ccmo(0,1,theta) * P_gm.ccmo(0,Xgm,theta,eta)
  res <- res + P_Y.ccmo(1,X,1,0,beta,nX) * P_gc.ccmo(1,0,theta) * P_gm.ccmo(1,Xgm,theta,eta)
  res <- res + P_Y.ccmo(1,X,1,1,beta,nX) * P_gc.ccmo(1,1,theta) * P_gm.ccmo(1,Xgm,theta,eta)
  res <- res + P_Y.ccmo(1,X,1,2,beta,nX) * P_gc.ccmo(1,2,theta) * P_gm.ccmo(1,Xgm,theta,eta)
  res <- res + P_Y.ccmo(1,X,2,1,beta,nX) * P_gc.ccmo(2,1,theta) * P_gm.ccmo(2,Xgm,theta,eta)
  res <- res + P_Y.ccmo(1,X,2,2,beta,nX) * P_gc.ccmo(2,2,theta) * P_gm.ccmo(2,Xgm,theta,eta)
  return(res)
}

likeli.ccmo.na <- function(para,Y,X,gm,gc,f,lambda,n,nX){
  res <- 0
  nbeta <- 3 + nX[1] + nX[2] + nX[3]
  beta <- para[1:nbeta]
  theta <- exp(para[nbeta+1]) / (1 + exp(para[nbeta+1]))
  eta <- para[(nbeta+2):(nbeta+1+nX[4])]
  for(i in 1:length(gm)){
    res <- res + log(likeli0.ccmo.na(Y[i],X[i,],gm[i],gc[i],beta,theta,eta,nX))
    res <- res - log(n * (1 + lambda * (H.ccmo(X[i,],beta,theta,eta,nX) - f)))
  }
  return(res)
}

likeli0.ccmo.na <- function(Y,X,gm,gc,beta,theta,eta,nX){
  Xgm <- X[-(1:(nX[1]+nX[2]+nX[3]))]
  gg <- gg.na(gm,gc)
  res <- 0
  for(i in 1:nrow(gg))
    res <- res + P_Y.ccmo(Y,X,gg[i,1],gg[i,2],beta,nX) * P_gc.ccmo(gg[i,1],gg[i,2],theta) * P_gm.ccmo(gg[i,1],Xgm,theta,eta)
  return(res)
}

gg.na <- function(gm,gc){
  gg <- NULL
  if(all(is.na(c(gm,gc)) == c(1,1)))
    for(i in 0:2)
      for(j in 0:2)
        gg <- rbind(gg,c(i,j))
  if(all(is.na(c(gm,gc)) == c(0,1)))
    for(j in 0:2)
      gg <- rbind(gg,c(gm,j))
  if(all(is.na(c(gm,gc)) == c(1,0)))
    for(i in 0:2)
      gg <- rbind(gg,c(i,gc))
  return(gg)
}
