# -------------------------------------------------------------------------------
# Tuning function for SGCCA 
# -------------------------------------------------------------------------------
sgcca_tune <- function(dat.list, nperms = 5, scheme = "centroid"){
  ### Inspired by tuning function in PMACCA
  len <- length(dat.list)
  min_dim <- dat.list %>% map(.,~ncol(.)) %>% unlist %>% min
  l1_grid <- c(seq(1/sqrt(min_dim) + 0.05 , 0.8, 0.05))

  t <- zs <- cor_org <- cv_cor <- c()
  
  for(i in c(1:length(l1_grid))){
    
    cor_perm <- c()
    s_org <- sgcca(dat.list, c1 = rep(l1_grid[i],2), scheme = scheme)
    cor_org <- c(cor_org, get_corr(s_org$Y, scheme = scheme))
    
    for(j in 1:nperms){
      
      dat.perm <- lapply(dat.list, function(d) {d[sample(nrow(d)),]})
      s_perm <- sgcca(dat.perm, c1 = rep(l1_grid[i],2), scheme = scheme)
      cor_perm <- c(cor_perm, get_corr(s_perm$Y, scheme = scheme))
      
    }
    
    cc.norm <- ftrans(cor_org)
    ccperm.norm <- ftrans(cor_perm)
    t <- c(t, mean(ccperm.norm >= cc.norm[i]))
    zs <- c(zs, (cc.norm[i]-mean(ccperm.norm))/(sd(ccperm.norm)+.05))
    
  }
  
  max.index <- which.max(zs)
  best_l1 <- l1_grid[max.index]
  
  return(list( best_l1 = best_l1, cor = cor_org[max.index], zs = zs, pval = t))
  
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
