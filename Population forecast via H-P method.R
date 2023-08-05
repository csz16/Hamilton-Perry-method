population.fore <- function(pop_f,pop_m,method,N,f.ratio,m.ratio,level){
  # intial value
  alpha <- 1-level
  
  pop_f.u <- pop_f
  pop_m.u <- pop_m
  pop_f.l <- pop_f
  pop_m.l <- pop_m
  orginal_f <- pop_f
  orginal_m <- pop_m
  
  for (k in 1:N) { 
    #one-step-ahead forecast
    n1 <- nrow(pop_f)-1
    m1 <- ncol(pop_f)-1
    
    # Compute the Cohort Component Ratio (CCR) for males
    CCR_m <- matrix(data = 0,nrow = n1,ncol = m1,dimnames = list(rownames(pop_f[-1,]),colnames(pop_m)[-1]))
    for (i in 1:n1) {
      for (j in 1:m1) {
        CCR_m[i,j] <- pop_m[i+1,j+1]/pop_m[i,j]
      }
    }
    
    # Compute the CCR for females
    CCR_f <- matrix(data = 0,nrow = n1,ncol = m1,dimnames = list(rownames(pop_m[-1,]),colnames(pop_f)[-1]))
    for (i in 1:n1) {
      for (j in 1:m1) {
        CCR_f[i,j] <- pop_f[i+1,j+1]/pop_f[i,j]
      }
    } 
    
    # Calculate the Child-Woman Ratio (CWR)
    child <- colSums(pop_f[1,])+colSums(pop_m[1,])
    women <- as.numeric(colSums(pop_f[c(16:50),]))
    
    # Forecast the CWR    
    CWR.m <- forecast(auto.arima(as.numeric(child/women)), h = 1, level = level)
    CWR <- as.numeric(CWR.m$mean)
    
    # Obtain upper and lower bound of the forecast
    CWR.u <- as.numeric(CWR.m$upper)
    CWR.l <- as.numeric(CWR.m$lower)
    
    # Forecast the age 0-1 female and male populations
    fore.pop.f <- women[length(women)]*CWR*f.ratio
    fore.pop.m <- women[length(women)]*CWR*m.ratio
    
    # Forecast the upper and lower bounds for age 0-1 female and male populations
    fore.pop.f.u <- women[length(women)]*CWR.u*f.ratio
    fore.pop.f.l <- women[length(women)]*CWR.l*f.ratio
    
    fore.pop.m.u <- women[length(women)]*CWR.u*m.ratio
    fore.pop.m.l <- women[length(women)]*CWR.l*m.ratio
    
    m1.1 <- fts(0:99,CCR_m)
    m1.2 <- ftsm(y = m1.1, order = 6, method = method, weight = FALSE)
    m1.3 <- forecast.ftsm(object = m1.2,h = 1,level = level, pimethod = "nonparametric", seed = 10)
    
    # Forecast the female and male populations
    fore.m <- c(fore.pop.m,m1.3$mean$y*pop_m[-nrow(pop_m),ncol(pop_m)])
    
    # Forecast the upper and lower bounds for female and male populations
    upper.m <- apply(m1.3$bootsamp, 1, function(x) quantile(x,1-alpha/2))
    lower.m <- apply(m1.3$bootsamp, 1, function(x) quantile(x,alpha/2))
    
    fore.m.u <- c(fore.pop.m.u,upper.m*pop_m.u[-nrow(pop_m.u),ncol(pop_m.u)])
    fore.m.l <- c(fore.pop.m.l,lower.m*pop_m.l[-nrow(pop_m.l),ncol(pop_m.l)])
    
    f1.1 <- fts(0:99,CCR_f)
    f1.2 <- ftsm(y = f1.1, order = 6, method = method, weight = FALSE)
    f1.3 <- forecast.ftsm(object = f1.2,h = 1,level = level, pimethod = "nonparametric", seed = 10)
    
    fore.f <- c(fore.pop.f,f1.3$mean$y*pop_f[-nrow(pop_f),ncol(pop_f)])
    upper.l <- apply(f1.3$bootsamp, 1, function(x) quantile(x,1-alpha/2))
    lower.l <- apply(f1.3$bootsamp, 1, function(x) quantile(x,alpha/2))
    
    fore.f.u <- c(fore.pop.f.u,upper.l*pop_f.u[-nrow(pop_f.u),ncol(pop_f.u)])
    
    fore.f.l <- c(fore.pop.f.l,lower.l*pop_f.l[-nrow(pop_f.l),ncol(pop_f.l)])
    
    #Combine the data set
    pop_f <- cbind(pop_f,fore.f)
    pop_m <- cbind(pop_m,fore.m)
    
    pop_f.u <- cbind(pop_f.u,fore.f.u)
    pop_m.u <- cbind(pop_m.u,fore.m.u)
    
    pop_f.l <- cbind(pop_f.l,fore.f.l)
    pop_m.l <- cbind(pop_m.l,fore.m.l)
    
    names(pop_f)[ncol(pop_f)] <- as.character(as.numeric(colnames(pop_f[ncol(pop_f)-1]))+1)
    names(pop_m)[ncol(pop_m)] <- as.character(as.numeric(colnames(pop_m[ncol(pop_m)-1]))+1)
    
    names(pop_f.u)[ncol(pop_f.u)] <- as.character(as.numeric(colnames(pop_f[ncol(pop_f)-1]))+1)
    names(pop_m.u)[ncol(pop_m.u)] <- as.character(as.numeric(colnames(pop_m.u[ncol(pop_m.u)-1]))+1)
    
    names(pop_f.l)[ncol(pop_f.l)] <- as.character(as.numeric(colnames(pop_f.l[ncol(pop_f.l)-1]))+1)
    names(pop_m.l)[ncol(pop_m.l)] <- as.character(as.numeric(colnames(pop_m.l[ncol(pop_m.l)-1]))+1)
  }
  
  pop.t <- pop_f+pop_m
  pop.t.1 <- pop.t[c(16:101),c(102:ncol(pop.t))]
  n2 <- nrow(pop.t.1)-1
  m2 <- ncol(pop.t.1)
  OADR <- matrix(data = 0,nrow = n2,ncol = m2,dimnames = list(rownames(pop.t.1[-1,]),colnames(pop.t.1)))
  
  for (j in 1:m2) {
    for (i in 1:n2) {
      OADR[i,j] <- sum(pop.t.1[c((i+1):(n2+1)),j])/sum(pop.t.1[c(1:i),j])
    }
  }
  
  #Calculate the sustainable pension age
  first_below <- apply(OADR, 2, function(x) min(which(x < 0.23)))
  pension.age <- substr(row.names(pop.t.1[first_below,]), start = 1, stop = 2)
  final.result <- cbind(colnames(pop.t.1), pension.age)
  colnames(final.result) <- c("year", "pension.age")
  ######################################## 
  #Calculate the sustainable lower bound pension age
  pop.t.l <- pop_f.l+pop_m.l
  pop.t.1l <- pop.t.l[c(16:101),c(102:ncol(pop.t.l))]
  n2 <- nrow(pop.t.1l)-1
  m2 <- ncol(pop.t.1l)
  OADR.l <- matrix(data = 0,nrow = n2,ncol = m2,dimnames = list(rownames(pop.t.1l[-1,]),colnames(pop.t.1l)))
  
  for (j in 1:m2) {
    for (i in 1:n2) {
      OADR.l[i,j] <- sum(pop.t.1l[c((i+1):(n2+1)),j])/sum(pop.t.1l[c(1:i),j])
    }
  }
  
  f.l <- apply(OADR.l, 2, function(x) min(which(x < 0.23)))
  p.l <- substr(row.names(pop.t.1l[f.l,]), start = 1, stop = 2)
  final.result.l <- cbind(colnames(pop.t.1l), p.l)
  colnames(final.result.l) <- c("year", "pension.age")
  ############################################### 
  #Calculate the sustainable upper bound pension age
  pop.t.u <- pop_f.u+pop_m.u
  pop.t.1u <- pop.t.u[c(16:101),c(102:ncol(pop.t.u))]
  n2 <- nrow(pop.t.1u)-1
  m2 <- ncol(pop.t.1u)
  OADR.u <- matrix(data = 0,nrow = n2,ncol = m2,dimnames = list(rownames(pop.t.1u[-1,]),colnames(pop.t.1u)))
  
  for (j in 1:m2) {
    for (i in 1:n2) {
      OADR.u[i,j] <- sum(pop.t.1u[c((i+1):(n2+1)),j])/sum(pop.t.1u[c(1:i),j])
    }
  }
  
  f.u <- apply(OADR.u, 2, function(x) min(which(x < 0.23)))
  p.u <- substr(row.names(pop.t.1u[f.u,]), start = 1, stop = 2)
  final.result.u <- cbind(colnames(pop.t.1u), p.u)
  colnames(final.result.u) <- c("year", "pension.age")
  ###########################################
  
  
  result_list <- list(pop_m = pop_m, pop_f = pop_f, orginal_f = orginal_f, orginal_m = orginal_m, CCR_m = CCR_m, CCR_f = CCR_f, pop_f.u = pop_f.u, pop_m.u = pop_m.u, pop_f.l = pop_f.l, pop_m.l = pop_m.l, final.result = final.result, final.result.l = final.result.l, final.result.u = final.result.u)
  
  return(result_list)
}

#Final results
pop.30 <- population.fore(pop_f = new_f,pop_m = new_m,method = method,N = 30,f.ratio = f.ratio,m.ratio = m.ratio,level = 0.95)
