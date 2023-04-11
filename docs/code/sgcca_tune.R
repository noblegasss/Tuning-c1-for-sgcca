# -------------------------------------------------------------------------------
# Tuning function for SGCCA 
# -------------------------------------------------------------------------------
library(RGCCA)
library(MASS)
library(tidyverse)
library(purrr)
# -------------------------------------------------------------------------------
sgcca_tune <- function(dat.list, nperms = 5, scheme = "centroid"){
  ### Inspired by tuning function in PMACCA
  len <- length(dat.list)
  min_dim <- dat.list %>% purrr::map(.,~ncol(.)) %>% unlist %>% min
  l1_grid <- c(seq(1/sqrt(min_dim) + 0.01 , 0.7, 0.05))

  
  doParallel::registerDoParallel(10)
  
  para <- plyr::ldply(
    l1_grid,
    .fun = function(l){
      s_org <- sgcca(dat.list, c1 = rep(l,2), scheme = scheme)
      cor_org <- s_org$AVE$AVE_inner
      
      cor_perm <- plyr::llply(
        1:nperms,
        .fun = function(j){
          dat.perm <- lapply(dat.list, function(d) {d[sample.int(nrow(d)),]})
          
          s_perm <- sgcca(dat.perm, c1 = rep(l,2), scheme = scheme)
          #cor_perm <- c(cor_perm, get_corr(s_perm$Y, scheme = scheme))
          c_perm <- s_perm$AVE$AVE_inner
          
          return(c_perm)
        }
      ) %>% unlist() 
      
      cc.norm <- ftrans(cor_org)
      ccperm.norm <- ftrans(cor_perm)
      t <- mean(ccperm.norm >= cc.norm)
      zs <- (cc.norm-mean(ccperm.norm))/(sd(ccperm.norm)+.05)
      
      data.frame(best_l1 = l, cor = cor_org, zs = zs, cor_perm = mean(cor_perm))
      
    },.parallel = T
  )
  
  ## Instead of selecting the penalty term using the maximum z score obtained by the permutation test, we calculated the descending index of z scores
  decreasing.ind <- sapply(1:(nrow(para)-1), function(i) (para[i,"zs"] - para[i+1,"zs"]))
  para$decent.ind <- c(decreasing.ind, NA)
  max.index <- which.max(para$decent.ind)
  tab <- para[max.index,]
  
  return(tab)
  
}

get_corr = function(dat, scheme = scheme, len = 2){
  corr = 0
  i = 1
  while(i <= len){
    j = i + 1
    while(j <= len){
      if(scheme == "horst"){
        corr = corr + cor(dat[[i]] , dat[[j]])
      }
      if(scheme == "factorial"){
        corr = corr + cor(dat[[i]] , dat[[j]])^2
      }else corr = corr + abs(cor(dat[[i]] , dat[[j]]))
      j = j+1
    }
    i = i + 1
  }
  return(corr)
}

sgcca_cv <- function(dat.list, scheme = "centroid", K = 3){
  
  len <- length(dat.list)
  min_dim <- dat.list %>% purrr::map(.,~ncol(.)) %>% unlist %>% min
  l1_grid <- c(seq(1/sqrt(min_dim) + 0.05 , 0.7, 0.05))
  
  doParallel::registerDoParallel(10)
  
  para <- plyr::ldply(
    l1_grid,
    .fun = function(l){
      fold <- caret::createFolds(1:nrow(dat.list[[1]]), k = K)
      
      cv.cor <- plyr::llply(
        1:K,
        .fun = function(cv){
          dat.list.train <- dat.list %>% purrr::map(.,~.[fold[-cv] %>% unlist,])
          dat.list.test <- dat.list %>% purrr::map(.,~.[fold[[cv]],])
          
          s1 <- sgcca(dat.list.train, c1 = rep(l,2), scheme = scheme)
          new.var <- list(s1$a[[1]][,1] %*% t(dat.list.test[[1]]) %>% t(), s1$a[[2]][,1] %*% t(dat.list.test[[2]]) %>% t())
          get_corr(new.var, scheme = scheme) %>% as.numeric()
        }
      ) %>% unlist() %>% mean()
      
      data.frame(l1 = l, cor = cv.cor)
      
    },.parallel = T
  )
 
  idx <- which.max(para$cor)
  
  return(para[idx,])
  
}

ftrans <- function (x) {
  return(atanh(x))
}
# -------------------------------------------------------------------------------
# Simulation function
# -------------------------------------------------------------------------------
sim.linear2 = function(n, p, j = 2, mu.sd = 2, sigma = 0.3, psel = 20){
  X = list()
  w = list()
  mu = rnorm(n, mean = 0, sd = mu.sd) 
  A = c("X", "Y", "Z", "V", "W")
  for(i in 1:j){
    w0 = runif(psel, -1, 1) 
    w[[i]] = c(w0/sum(w0), rep(0, p - psel))
    e1 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma, nrow = p))
    X[[i]] = data.frame(mu %*% t(w[[i]]) + e1)
    colnames(X[[i]]) = sapply(1:p, function(m) paste0(A[i],m))
  }
  names(X) = sapply(1:j, function(i) paste0("X",i))
  return(list(X = X, weight = w))
}
