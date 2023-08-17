
# 

library(R.matlab)
library(glmnet)
# 
lambdas <- 10^seq(-10, 0, length=200)


N_nodes = 1000  #
N_folds = 100  # 
mode_seq = c(434:1000) #


# nested cross-validation outer loop
for (i in c(1:N_folds)) {
  # 
  train_data <- readMat(paste0('data/train_test',  i, '.mat'))
  
  # 初始化需要保存的变量
  sele_lambda <- rep(NA, N_nodes)  # 
  sele_mode <- matrix(0, nrow = N_nodes, ncol = length(mode_seq)) #
  acc <- rep(NA, N_nodes)   #
  # 
  x_train_out <- NULL
  y_train_out <- NULL
  for (j in c(1:20)) {
    x_train_out <- rbind(x_train_out, train_data$BEC.trainK[, , j], train_data$BEC.testK[, , j])
    y_train_out <- rbind(y_train_out, train_data$FC.trainK[, , j], train_data$FC.testK[, , j])
  }
  x_train_out <- x_train_out[, mode_seq]
  
  x_test_out <- train_data$BEC.test
  x_test_out <- x_test_out[, mode_seq]
  y_test_out <- train_data$FC.test
  
  # inner-loop cross validation
  corr <- array(dim = c(20, length(lambdas), N_nodes))
  for (j in c(1:20)) {
    # 
    x_train <- train_data$BEC.trainK[, , j]
    x_train <- x_train[, mode_seq]
    y_train <- train_data$FC.trainK[, , j]
    x_test <- train_data$BEC.testK[, , j]
    x_test <- x_test[, mode_seq]
    y_test <- train_data$FC.testK[, , j]
    
    for (k in c(1:N_nodes)) {
      # 
      #x_train_k <- x_train[-k, ]
      x_train_k <- x_train
      # 
      y_train_k <- y_train[, k]
      # y_train_k <- y_train_k[-k]
      # 
      #x_test_k <- x_test[-k, ]
      x_test_k <- x_test
      # 
      y_test_k <- y_test[, k]
      # y_test_k <- y_test_k[-k]
      #      
      lasso_models <- glmnet(x_train_k, y_train_k, alpha = 1, lambda = lambdas, family = 'gaussian', intercept = TRUE)
      y_predict <- predict(lasso_models, x_test_k)
      corr[j, , k] <- cor(y_predict, y_test_k)
    }
  }
  # 
  for (j in c(1:N_nodes)) {
    corr_temp <- corr[, , j]
    # 
    d_mean <- colMeans(corr_temp)
    # 
    ind <- grep('TRUE', d_mean==max(d_mean, na.rm = TRUE))
    max_lambda <- lasso_models$lambda[ind[1]]
    sele_lambda[j] <- max_lambda
    # 
    # rmows <- seq(from=0, by=1000, length=20) + rep(j, 20)
    # x_train_out_j <- x_train_out[-rm_rows, ]
    x_train_out_j <- x_train_out
    y_train_out_j <- y_train_out[, j]
    # y_train_out_j <- y_train_out_j[-rm_rows]
    
    # 
    lasso_models1 <- glmnet(x_train_out_j, y_train_out_j, alpha = 1, lambda = max_lambda, family = 'gaussian', intercept = TRUE)
    sele_ind <- grep('TRUE', lasso_models1$beta != 0)
    sele_mode[j, sele_ind] <- 1
    #
    # x_test_out_j <- x_test_out[-j, ]
    x_test_out_j <- x_test_out
    y_test_out_j <- y_test_out[, j]
    # y_test_out_j <- y_test_out_j[-j]
    
    y_predict <- predict(lasso_models1, x_test_out_j)
    acc[j] <- cor(y_predict, y_test_out_j)
  }
  # 
  save_path <- 'res/'
  # 
  write.table(sele_lambda, file = paste0(save_path, 'lambda_', i, '.txt'), quote = FALSE,row.names = FALSE,col.names = FALSE,sep = ' ')
  write.table(sele_mode, file = paste0(save_path, 'selected_mode_', i, '.txt'), quote = FALSE,row.names = FALSE,col.names = FALSE,sep = ' ')
  write.table(acc, file = paste0(save_path, 'acc_', i, '.txt'), quote = FALSE,row.names = FALSE,col.names = FALSE,sep = ' ')
}



















