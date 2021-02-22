##################################################################################
### Quality Control functions
##################################################################################

runQC <- function(filepath){
  txt_folder <-  file.path(filepath)
  r = createReport(txt_folder)
  cat(paste0("\nReport generated as '", r$report_file, "'\n\n"))
}


#####################################################################################
### Impute
#####################################################################################

count_zeros <- function(x){
  return(sum(is.na(x)))
}

min_col <- function(x){
  return(names(x[x==min(x)]))
}
### function: generate uniform distribution for given number, mean and std
na_zeros_impute <- function(na_zeros , mu,sd){
  return(matrix(runif(sum(na_zeros), mu-3*sd, mu-2*sd),ncol=1))
}

calculate_stats_nonzeros <- function(data){
  ## Predicting all missing replicates values
  colnames <- colnames(data)
  sample <- data
  sample[sample==0] <- NA
  if(length(colnames)==2){
    print(paste0("Lenght of the columns is ",length(colnames)))
    
    mu1 <- mean(sample[,1], na.rm=T)
    sd1 <- sd(sample[,1], na.rm=T)
    mu2 <- mean(sample[,2], na.rm=T)
    sd2 <- sd(sample[,2], na.rm=T)
    stats <- rbind(cbind(mu1,sd1),
                   cbind(mu2,sd2))
  }
  if(length(colnames)==3){
    print(paste0("Lenght of the columns is ",length(colnames)))
    
    mu1 <- mean(sample[,1], na.rm=T)
    sd1 <- sd(sample[,1], na.rm=T)
    mu2 <- mean(sample[,2], na.rm=T)
    sd2 <- sd(sample[,2], na.rm=T)
    mu3 <- mean(sample[,3], na.rm=T)
    sd3 <- sd(sample[,3], na.rm=T)
    stats <- rbind(cbind(mu1,sd1),
                   cbind(mu2,sd2),
                   cbind(mu3,sd3))
  }
  if(length(colnames)==4){
    print(paste0("Lenght of the columns is ",length(colnames)))
    
    mu1 <- mean(sample[,1], na.rm=T)
    sd1 <- sd(sample[,1], na.rm=T)
    mu2 <- mean(sample[,2], na.rm=T)
    sd2 <- sd(sample[,2], na.rm=T)
    mu3 <- mean(sample[,3], na.rm=T)
    sd3 <- sd(sample[,3], na.rm=T)
    mu4 <- mean(sample[,4], na.rm=T)
    sd4 <- sd(sample[,4], na.rm=T)
    stats <- rbind(cbind(mu1,sd1),
                   cbind(mu2,sd2),
                   cbind(mu3,sd3),
                   cbind(mu4,sd4))
  }
  
  return(stats)
}

impute_all_zeros <- function(data,amm = "2", pmm ="6"){
  ## Predicting all missing replicates values
  colnames <- colnames(data)
  na_zeros <- rowSums(data,na.rm=T)==0
  print(paste("There are", length(na_zeros), "total zero records in either caseor control replicates to be imputed."))
  if(sum(na_zeros)>0){
    min_zero_col <- min_col(apply(data,2,count_zeros))
    if(amm == "2"){
      sample_min_zeros <- unname(unlist(data[,min_zero_col]))
      mu <- mean(sample_min_zeros, na.rm=T)
      sd <- sd(sample_min_zeros, na.rm=T)
      data[na_zeros,colnames] <- matrix(na_zeros_impute(na_zeros* length(colnames) , mu,sd),ncol= length(colnames) , nrow = sum(na_zeros))
    } 
    if(amm == "1"){
      sample <- data
      stats <- calculate_stats_nonzeros(data)
      
      dat <- matrix(NA,nrow= sum(na_zeros),ncol=0)
      for(i in 1:dim(stats)[1]){
        dat <- cbind(dat,na_zeros_impute(na_zeros , stats[1,1],stats[1,2]))
      }
      data[na_zeros,colnames] <- matrix(dat,
                                        ncol= length(colnames),
                                        nrow = sum(na_zeros))
    }
  }
  return(data)
}


impute_partial_zeros <- function(data, count_na, colnames, nb=0){
  if(length(colnames)==2){
    y <- data[count_na>nb,colnames]
  } else{
    y <- data[count_na>nb,colnames]
  }
  
  if(length(colnames)>2){
    for(i in 1:dim(y)[1]){
      #print(i)
      col_NAs <- colnames(y)[is.na(unlist(y[i,]))]
      for(j in col_NAs){
        col_NA <- j
        col_select <- names(y[i,colnames(y)!=col_NA])
        sample <- data[complete.cases(data[,col_select]),col_select]
        delta <- (sample[,1]-sample[,2])/mean(unlist(unname(sample)))
        mu <- mean(unlist(unname(delta)), na.rm=T)
        std <- sd(unlist(unname(delta)),na.rm=T)
        cor <- cor(data[count_na==0,colnames])
        mean_cor <- mean(cor[col_NA,col_select])
        deltanew <- rnorm(1,mu, std/(sqrt(2)*mean_cor))
        y[i,col_NA] <- mean(unlist(unname(y[i,col_select])),na.rm=T)*abs(1+deltanew)
      }
    }
  }else{
    for(i in 1:dim(y)[1]){
      col_NAs <- colnames(y)[is.na(unlist(y[i,]))]
      for(j in col_NAs){
        col_NA <- j
        col_select <- names(y)
        sample <- data[complete.cases(data[,col_select]),col_select]
        delta <- (sample[,1]-sample[,2])/mean(unlist(unname(sample)))
        mu <- mean(unlist(unname(delta)), na.rm=T)
        std <- sd(unlist(unname(delta)),na.rm=T)
        cor <- cor(data[count_na==0,colnames])
        mean_cor <- mean(cor[col_NA,col_select])
        deltanew <- rnorm(1,mu, std/(sqrt(2)*mean_cor))
        y[i,col_NA] <- mean(unlist(unname(y[i,col_select])),na.rm=T)*abs(1+deltanew)
      }
    }
  }
  
  if(length(colnames)==2){
    data[count_na>nb,colnames] <- y
  }else{
    data[count_na>nb,colnames] <- y
  }
  
  return(data)
}


impute <- function(data, amm = "2", pmm = "6"){
  data[data==0] <- NA
  stats <- calculate_stats_nonzeros(data)
  data <- impute_all_zeros(data,amm )
  colnames <- colnames(data)
  count_na <- apply(data,1, count_zeros)
  
  if(pmm == "7"){
    data <- impute_partial_zeros(data, count_na, colnames, nb = 0)
  }
  
  if(pmm == "6"){
    data <- impute_partial_zeros(data, count_na, colnames, nb = 1)
    
    for( i in 1:length(colnames)){
      data_missing <- unlist(lapply(data[,colnames[i]],  function(x) is.na(x)))
      data[data_missing,colnames[i]] <- na_zeros_impute(sum(data_missing) , stats[i,1], stats[i,2])
    }
  }
  
  return(data)
}


############################