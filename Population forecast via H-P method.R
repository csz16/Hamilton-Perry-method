require(demography)
require(ftsa)
require(TSA)

linear_interpolated_pop <- function(data, new_age, age, forecast_year) {
  # Initialize an empty list to store interpolated data for each year
  interpolated_data_list <- list()
  
  # Iterate over each year (column in the matrix)
  for (year_index in 1:ncol(data)) {
    year_data <- data[, year_index]
    interpolated_population <- approx(x = age, y = year_data, xout = new_age)$y
    interpolated_data_list[[year_index]] <- interpolated_population
  }
  
  # Convert the list to a matrix
  # Each column in this new matrix will correspond to a year, as in the original matrix
  interpolated_matrix <- do.call(cbind, interpolated_data_list)
  if (ncol(interpolated_matrix)==length(forecast_year)) {
    colnames(interpolated_matrix) <- forecast_year
  }
  rownames(interpolated_matrix) <- new_age
  return(interpolated_matrix)
}


pension_age_interpolated <- function(data)
{
  pop_total_1 <- data[c(which(rownames(data)==16):nrow(data)),]
  n2 <- nrow(pop_total_1)-1
  m2 <- ncol(pop_total_1)
  OADR <- matrix(data = 0, nrow = n2, ncol = m2, dimnames = list(rownames(pop_total_1[-1,]),colnames(pop_total_1)))
  
  for(j in 1:m2){
    for(i in 1:n2) 
    {
      OADR[i,j] <- sum(pop_total_1[c((i+1):(n2+1)),j])/sum(pop_total_1[c(1:i),j])
    }
  }
  
  result <- apply(OADR, 2, function(x) min(which(x < 0.23)))
  pension_age <- row.names(pop_total_1[result,])
  final_result <- cbind(colnames(pop_total_1), pension_age)
  colnames(final_result) <- c("year", "pension.age")
  final_result <- as.data.frame(final_result)
  final_result[] <- lapply(final_result, as.numeric)
  result_list <- list(final_result = final_result, OADR = OADR)
  return(result_list)
}



pension_age_bound_interpolated <- function(data, alpha)
{
  pop_total_1 <- data[c(which(rownames(data)==16):nrow(data)),]
  n2 <- nrow(pop_total_1)-1
  m2 <- ncol(pop_total_1)
  OADR <- matrix(data = 0,nrow = n2,ncol = m2)
  rownames(OADR) <- rownames(data[-1])
  
  for(j in 1:m2) 
  {
    for(i in 1:n2) 
    {
      OADR[i,j] <- sum(pop_total_1[c((i+1):(n2+1)),j])/sum(pop_total_1[c(1:i),j])
    }
  }
  
  result <- apply(OADR, 2, function(x) min(which(x < 0.23)))
  pension_age <- as.numeric(row.names(pop_total_1[result,]))
  final_result <- quantile(pension_age, alpha)
  return(final_result)
}

population_fore_interpolated <- function(pop_f, pop_m, method, N, f.ratio, m.ratio, level, year, age)
{
  # intial value
  alpha <- 1-level
  pop_f.u <- pop_f
  pop_m.u <- pop_m
  pop_f.l <- pop_f
  pop_m.l <- pop_m
  original_f <- pop_f
  original_m <- pop_m
  last_year <- year[length(year)]
  forecast_year <- c((last_year+1):(last_year+N))
  fore_boost_f <- list()
  fore_boost_m <- list()
  
  for(k in 1:N) 
  {   
    if(any(original_f == 0))
    {
      # Add 0.01 to all data points
      pop_f <- pop_f + 0.01
    } 
    
    if(any(original_m == 0))
    {
      # Add 0.01 to all data points
      pop_m <- pop_m + 0.01
    } 
    #one-step-ahead forecast
    n1 <- nrow(pop_f)-1
    m1 <- ncol(pop_f)-1
    
    # Compute the Cohort Component Ratio (CCR) for males
    CCR_m <- matrix(data = 0, nrow = n1, ncol = m1, dimnames = list(rownames(pop_f[-1,]),colnames(pop_m)[-1]))
    for(i in 1:n1) 
    {
      for(j in 1:m1) 
      {
        CCR_m[i,j] <- pop_m[i+1,j+1]/pop_m[i,j]
      }
    }
    
    # Compute the CCR for females
    CCR_f <- matrix(data = 0,nrow = n1,ncol = m1,dimnames = list(rownames(pop_m[-1,]),colnames(pop_f)[-1]))
    for(i in 1:n1) 
    {
      for(j in 1:m1) 
      {
        CCR_f[i,j] <- pop_f[i+1,j+1]/pop_f[i,j]
      }
    }
    
    # Calculate the Child-Woman Ratio (CWR)
    child <- pop_f[1,]+pop_m[1,]
    women <- as.numeric(colSums(pop_f[c(16:50),]))
    
    # Bootstrap forecast the CWR    
    CWR_bootstrap <- replicate(100, simulate(auto.arima(as.numeric(child/women)), nsim = 1))
    
    # Forecast the age 0-1 female and male populations
    fore.pop.f <- women[length(women)]*mean(CWR_bootstrap)*f.ratio
    fore.pop.m <- women[length(women)]*mean(CWR_bootstrap)*m.ratio
    
    m1.1 <- forecast.ftsm(ftsm(fts(0:(nrow(pop_m)-2),CCR_m), order = 6, method = method, weight = FALSE), h = 1,level = level, pimethod = "nonparametric")
    
    fore.m <- c(fore.pop.m, m1.1$mean$y*pop_m[-nrow(pop_m),ncol(pop_m)])
    fore_boost_m[[k]] <- rbind(women[length(women)]*CWR_bootstrap*m.ratio, m1.1$bootsamp[1:nrow(CCR_m), ,1]*pop_m[-nrow(pop_m),ncol(pop_m)])
    
    f1.1 <- forecast.ftsm(ftsm(fts(0:(nrow(pop_f)-2),CCR_f), order = 6, method = method, weight = FALSE), h = 1,level = level, pimethod = "nonparametric")
    
    # Forecast the upper and lower bounds for female and male populations
    fore.f <- c(fore.pop.f, f1.1$mean$y*pop_f[-nrow(pop_f),ncol(pop_f)])
    fore_boost_f[[k]] <- rbind(women[length(women)]*CWR_bootstrap*f.ratio, f1.1$bootsamp[1:nrow(CCR_f), ,1]*pop_f[-nrow(pop_f),ncol(pop_f)])
    
    #Combine the data set (Hanlin suggests below should be removed)
    
    pop_f <- cbind(pop_f, fore.f)
    pop_m <- cbind(pop_m, fore.m)
    
    pop_f.u <- cbind(pop_f.u, apply(fore_boost_f[[k]], 1, function(x) quantile(x,1-alpha/2)))
    pop_m.u <- cbind(pop_m.u, apply(fore_boost_m[[k]], 1, function(x) quantile(x,1-alpha/2)))
    
    pop_f.l <- cbind(pop_f.l, apply(fore_boost_f[[k]], 1, function(x) quantile(x,alpha/2)))
    pop_m.l <- cbind(pop_m.l, apply(fore_boost_m[[k]], 1, function(x) quantile(x,alpha/2)))
    
    colnames(pop_f)[ncol(pop_f)] <- as.character(as.numeric(colnames(pop_f)[ncol(pop_f)-1])+1)
    colnames(pop_m)[ncol(pop_m)] <- as.character(as.numeric(colnames(pop_m)[ncol(pop_m)-1])+1)
    
    colnames(pop_f.u)[ncol(pop_f.u)] <- as.character(as.numeric(colnames(pop_f)[ncol(pop_f)-1])+1)
    colnames(pop_m.u)[ncol(pop_m.u)] <- as.character(as.numeric(colnames(pop_m.u)[ncol(pop_m.u)-1])+1)
    
    colnames(pop_f.l)[ncol(pop_f.l)] <- as.character(as.numeric(colnames(pop_f.l)[ncol(pop_f.l)-1])+1)
    colnames(pop_m.l)[ncol(pop_m.l)] <- as.character(as.numeric(colnames(pop_m.l)[ncol(pop_m.l)-1])+1)
    
    if(any(original_f == 0))
    {
      # Subtract 0.01 to all data points
      pop_f <- pop_f - 0.01
      pop_f.u <- pop_f.u - 0.01
      pop_f.l <- pop_f.l - 0.01
    } 
    
    if(any(original_m == 0))
    {
      # Subtract 0.01 to all data points
      pop_m <- pop_m - 0.01
      pop_m.u <- pop_m.u - 0.01
      pop_m.l <- pop_m.l - 0.01
    } 
    print(k); rm(k)
  }
  
  new_age <- seq(0, 100, by = 0.1)
  
  fore_pop_total_mean <- pop_m[,which(colnames(pop_m)%in%forecast_year)] + pop_f[,which(colnames(pop_f)%in%forecast_year)]
  
  interpolated_fore_total_mean <- linear_interpolated_pop(data = fore_pop_total_mean, new_age = new_age, age = age, forecast_year = forecast_year)
  mean <- pension_age_interpolated(interpolated_fore_total_mean)
  
  
  fore_total_boost_list <- list()
  interpolated_fore_total_boost_list <- list()
  
  for (j in 1:N) {
    fore_total_boost_list[[j]] <- fore_boost_f[[j]]+fore_boost_m[[j]]
    interpolated_fore_total_boost_list[[j]] <- linear_interpolated_pop(data = fore_total_boost_list[[j]] , new_age = new_age, age = age, forecast_year = forecast_year)
  }

  
  
  
  pension_age_upper <- pension_age_lower <- pension_age_median <- vector("numeric", N)
  for(i in 1:N) 
  {
    pension_age_upper[i] <- pension_age_bound_interpolated(data = interpolated_fore_total_boost_list[[i]], alpha = 1-alpha/2)
    # below is the median
    pension_age_median[i] <- pension_age_bound_interpolated(data = interpolated_fore_total_boost_list[[i]], alpha = 0.5)
    pension_age_lower[i] <- pension_age_bound_interpolated(data = interpolated_fore_total_boost_list[[i]], alpha = alpha/2)
  }
  median <- cbind(forecast_year, pension_age_median)
  lower <- cbind(forecast_year, pension_age_lower)
  upper <- cbind(forecast_year, pension_age_upper)
  
  ###########################################
  
  result_list <- list(fore_pop_m = pop_m, 
                      fore_pop_f = pop_f, 
                      original_f = original_f, 
                      original_m = original_m, 
                      CCR_m = CCR_m, 
                      CCR_f = CCR_f, 
                      pop_f.u = pop_f.u, 
                      pop_m.u = pop_m.u, 
                      pop_f.l = pop_f.l, 
                      pop_m.l = pop_m.l, 
                      mean = mean$final_result, 
                      median = median,
                      upper = upper, 
                      lower = lower)
  
  return(result_list)
}


pop_30_au_interpolated <- list()
for (i in 1:8) {
  pop_30_au_interpolated[[i]] <- population_fore_interpolated(pop_f = female_matrices[[i]], pop_m = male_matrices[[i]], method = "classical", N = 30, f.ratio = f.ratio, m.ratio = m.ratio, level = 0.95, year = year, age = age)
  print(i); rm(i)
}


names(pop_30_au_interpolated) <- paste(valid_States, 30, sep = "_")


pop_30_total_au_interpolated <- population_fore_interpolated(pop_f = au_pop_f, pop_m = au_pop_m, method = "classical", N = 30, f.ratio = f.ratio, m.ratio = m.ratio, level = 0.95, year = year, age = age)
