# This function is a modification of the SLOG algorithm to compute Lasso coefficients (described in Rajaratnam, 2016). 
# The modification adds adaptive weights as well as a cross-validation component to the original algorithm sequence.

adalasso.SLOG.withCV <- function(x,y,times,thresh,start=NULL){

  # x: covariate data
  # y: response data
  # lambda.seq: range of values of lasso regularization parameter
  # times: convergence criteria - difference between successive coefficient vectors
  # thresh: the value below which estimates are set to 0 (runs rSLOG)
  # start: allows starting values other than sign(xty)*lambda.seq/p to be specified.
  # wts: vector of weights used for adaptive lasso stage, obtained from previous stage.

  
  
  # Define vectors to save values of B and k over all values of lambda chosen
  
  b.lambda <- NULL
  k.lambda <- NULL
  mse.lambda <- NULL
  # Prepare sequence of lambda to try over
  
  lambda.seq <- seq (0.01,10,0.01)
  l.ada <- NULL
  # Loop over lambda
      
  for (j in 1:length(lambda.seq)) {
  
    # Randomize data
    full.data <- cbind(x, y)
    full.data<-full.data[sample(nrow(x)), ]
  
    # Create 10 equally size folds
    folds <- cut(seq(1,nrow(full.data)),breaks=10,labels=FALSE)
  
    # Prepare vectors for storing values from cross-validation
    B.CUR.CV <- NULL
    K.CUR.CV <- NULL
    xty <- NULL
    coeffs.stage1 <- NULL
    mse.cv <- NULL
    # Perform 10 fold cross validation for each value of lambda
    for (i in 1:10) {
      testIndexes <- which(folds==i,arr.ind=TRUE)
      testData <- full.data[testIndexes, ]
      trainData <- full.data[-testIndexes, ]
      x.values <- data.matrix(full.data[-testIndexes, 1:(ncol(full.data)-1)])
      y.values <- data.matrix(full.data[-testIndexes, ncol(full.data)])
      xtx<-crossprod(x.values)
      xty<-crossprod(x.values,y.values)
      p<-length(xty)
      n<-length(y.values)
  
      # First, we find coefficients using OLS and create a weights vector inversely proportional to OLS coefficients.
  
      model.stage1 <- lm(y.values~x.values)
      coeffs.stage1 <- model.stage1$coefficients[-1]
      wts <- 1/abs(as.vector(coeffs.stage1))
      
      # Following is the SLOG algorithm with the modified lambda value
  
      l.ada <- lambda.seq[j]*wts
      b.cur <- sign(xty)*l.ada/p
      b.old <- b.cur
      vin<-1:p
      temp<-rep(0,p)
      conv=FALSE
      k<-1
      while(conv==FALSE){
      p <- length(b.cur[vin])
      
      # Following is the deterministic sequence that defines the coefficients.
      
      B.inv  <- diag(l.ada/abs(b.cur[vin] + 0001), nrow = p, ncol = p)
      b.new<-tcrossprod(chol2inv(chol(B.inv+xtx[vin,vin])),t(xty[vin]))
      temp[vin]<-as.vector(b.new)
      temp[abs(temp)<=thresh]<-0
      b.new<-temp
      vin<-which(b.new!=0)
      b.cur <- as.vector(b.new)
      conv<-(sqrt(sum((b.cur-b.old)^2))/sqrt(sum(b.old^2)))<times
      b.old<-b.cur
      k<-k+1
      }
      # save values from cross-validation into these objects
  
      B.CUR.CV <- cbind(B.CUR.CV, b.cur)
      K.CUR.CV <- c(K.CUR.CV, k)
  

      response.pred <- NULL
      # Find mean square error of each loop of cross-validation
      for (m in 1:nrow(testData)){
      response.pred[m] <- sum(b.cur*testData[m,1:(ncol(full.data)-1)])
      }
      resid <- full.data[(testIndexes), ncol(full.data)] - response.pred
      mse.cv <- c(mse.cv, mean(resid^2))

    } # END OF CROSS VALIDATION LOOP
  
    # average over cv objects to find values for each particular lambda
  
    mse.lambda[j] <- mean(mse.cv)
    b.lambda <- rbind(b.lambda, rowMeans(B.CUR.CV))
    k.lambda[j] <- mean(K.CUR.CV)
      
  } # END OF LOOP OVER ALL VALUES OF LAMBDA
  
  # Finally, pick the value of b for which MSE is minimum
  
  b.final <- b.lambda[which.min(mse.lambda), ]
  k.final <- k.lambda[which.min(mse.lambda)]
  lambda.final <- lambda.seq[which.min(mse.lambda)]
  return(c(lambda.final, k.final, b.final))
}
  

