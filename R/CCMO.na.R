CCMO.na <- function(Y,gm,gc,X,f){
  n1 <- sum(Y == 1)
  n0 <- sum(Y == 0)
  n <- n1+n0
  lambda <- n1 / (n*f) - n0 / (n * (1-f))
  X <- as.matrix(X)
  nX <- max(1,ncol(X))
  
  group2 <- apply(is.na(cbind(gm,gc)) , 1, any)
  X1 <- as.matrix(X[!group2,])
  Y1 <- Y[!group2]
  gm1 <- gm[!group2]
  gc1 <- gc[!group2]
  X2 <- as.matrix(X[group2,])
  Y2 <- Y[group2]
  gm2 <- gm[group2]
  gc2 <- gc[group2]
  
  Z <- cbind(1,gm1,gc1,X1,gm1*X1,gc1*X1)
  fit = glm(Y1 ~ 0 + Z,family = binomial)
  res = summary(fit)$coef
  est.log = as.vector(res[,1])
  sd.log = as.vector(res[,2])
  
  beta0 <- est.log
  theta0 <- (2*sum(gm1 == 2)+sum(gm1 == 1))/(2*length(gm1))
  theta0 <- log(theta0/(1-theta0))
  eta0 <- 0
  para0 <- c(beta0,theta0,eta0)
  
  llik <- function(para)
    -likeli.ccmo(para,Y1,X1,gm1,gc1,f,lambda,n) - likeli.ccmo.na(para,Y2,X2,gm2,gc2,f,lambda,n)
  fit <- optim(par = para0,fn = llik,method = 'L-BFGS-B',hessian = TRUE)
  est <- fit$par #beta 6, theta 1, eta 1
  est <- c(est[1:(3+nX*3)],exp(est[4+nX*3])/(1 + exp(est[4+nX*3])),est[(5+nX*3):(4+nX*4)])
  sd <- sqrt(diag(solve(fit$hessian)))[1:(3+nX*3)]
  return(list(est = est,sd = sd,est.log = est.log,sd.log = sd.log))
}

likeli.ccmo <- function(para,Y,X,gm,gc,f,lambda,n){
  res <- 0
  nX <- max(1,ncol(X))
  beta <- para[1:(3+nX*3)]
  theta <- exp(para[4+nX*3]) / (1 + exp(para[4+nX*3]))
  eta <- para[(5+nX*3):(4+nX*4)]
  for(i in 1:length(gm)){
    res <- res + log(P_Y.ccmo(Y[i],X[i,],gm[i],gc[i],beta)) + log(P_gc.ccmo(gm[i],gc[i],theta)) + log(P_gm.ccmo(gm[i],X[i,],theta,eta))
    res <- res - log(n * (1 + lambda * (H.ccmo(X[i,],beta,theta,eta) - f)))
  }
  return(res)
}

P_Y.ccmo <- function(Y,X,gm,gc,beta){
  tmp <- exp(sum(c(1,gm,gc,X,gm*X,gc*X) * beta))
  return(1-Y+(2*Y-1)*tmp/(1+tmp))
}

P_gc.ccmo <- function(gm,gc,theta){
  Pg <- matrix(c(1-theta,theta,0,(1-theta)/2,0.5,theta/2,0,1-theta,theta),3,3)
  return(Pg[gc+1,gm+1])
}

P_gm.ccmo <- function(gm,X,theta,eta){
  Pgm <- c((1-theta)^2,2*theta*(1-theta)*exp(sum(eta*X)),theta*theta*exp(2*sum(eta*X)))
  Pgm <- Pgm/sum(Pgm)
  return(Pgm[gm+1])
}

H.ccmo <- function(X,beta,theta,eta){
  res <- P_Y.ccmo(1,X,0,0,beta) * P_gc.ccmo(0,0,theta) * P_gm.ccmo(0,X,theta,eta)
  res <- res + P_Y.ccmo(1,X,0,1,beta) * P_gc.ccmo(0,1,theta) * P_gm.ccmo(0,X,theta,eta)
  res <- res + P_Y.ccmo(1,X,1,0,beta) * P_gc.ccmo(1,0,theta) * P_gm.ccmo(1,X,theta,eta)
  res <- res + P_Y.ccmo(1,X,1,1,beta) * P_gc.ccmo(1,1,theta) * P_gm.ccmo(1,X,theta,eta)
  res <- res + P_Y.ccmo(1,X,1,2,beta) * P_gc.ccmo(1,2,theta) * P_gm.ccmo(1,X,theta,eta)
  res <- res + P_Y.ccmo(1,X,2,1,beta) * P_gc.ccmo(2,1,theta) * P_gm.ccmo(2,X,theta,eta)
  res <- res + P_Y.ccmo(1,X,2,2,beta) * P_gc.ccmo(2,2,theta) * P_gm.ccmo(2,X,theta,eta)
  return(res)
}

likeli.ccmo.na <- function(para,Y,X,gm,gc,f,lambda,n){
  res <- 0
  nX <- max(1,ncol(X))
  beta <- para[1:(3+nX*3)]
  theta <- exp(para[4+nX*3]) / (1 + exp(para[4+nX*3]))
  eta <- para[(5+nX*3):(4+nX*4)]
  for(i in 1:length(gm)){
    res <- res + log(likeli0.ccmo.na(Y[i],X[i,],gm[i],gc[i],beta,theta,eta))
    res <- res - log(n * (1 + lambda * (H.ccmo(X[i,],beta,theta,eta) - f)))
  }
  return(res)
}

likeli0.ccmo.na <- function(Y,X,gm,gc,beta,theta,eta){
  gg <- gg.na(gm,gc)
  res <- 0
  for(i in 1:nrow(gg))
    res <- res + P_Y.ccmo(Y,X,gg[i,1],gg[i,2],beta) * P_gc.ccmo(gg[i,1],gg[i,2],theta) * P_gm.ccmo(gg[i,1],X,theta,eta)
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
