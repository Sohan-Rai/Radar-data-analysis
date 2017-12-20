trainfilenames=c()
testfilenames=c()
for (i in 1:5)
{
  trainfilename=paste("ionosphere_train_", i, ".csv",sep="")
  trainfilenames=c(trainfilenames,trainfilename)
  testfilename=paste("ionosphere_test_", i, ".csv",sep = "")
  testfilenames=c(testfilenames,testfilename)
}

level1='g'
level2='b'
levels=c(level1,level2)

euclidean_distance = function (a,b)
{
  dist = sqrt(apply(t(apply(b,1,'-',a) )^2,1,sum))
  return (dist)
}

cosine_distance = function (a,b)
{
  dist = abs((t(apply(b,1,'%*%',a)))-1)/(sqrt(a%*%a %*% t(rowSums(b^2))))
  return (dist)
}



knn = function(xtrain,ytrain,xtest,k,level)
{
  return_min_indices = function(a)
  {
    min_indices=which(a <= sort(a)[k], arr.ind=TRUE)
    return(min_indices)
  }
  euc_dist=matrix(NA,nrow=nrow(xtest),ncol=nrow(xtrain))
  cb_dist=matrix(NA,nrow=nrow(xtest),ncol=nrow(xtrain))
  xtest=as.matrix(xtest)
  for (i in 1:nrow(xtest))
  {
    euc_dist[i,]=euclidean_distance(xtest[i,],as.matrix(xtrain))
    cb_dist[i,]=cosine_distance(xtest[i,],as.matrix(xtrain))
  }
  xtest=as.data.frame(xtest)
  ##########################################################
  closest_euc=matrix(data=NA, nrow = nrow(euc_dist), ncol = nrow(xtrain))
  for (i in 1:nrow(euc_dist))
  {
    temp=return_min_indices(euc_dist[i,])
    for (j in 1:length(temp))
    {
      closest_euc[i,j]=temp[j]
    }
  }
  #a=euc_dist[1,]
  euc_count=matrix(data=0,nrow(xtest),ncol=length(levels))
  pred_euc=rep(0,nrow(xtest))
  for (i in 1:nrow(closest_euc))
  {
    for (x in 1:length(levels))
    {
      for (j in 1:k)
      {
        if (ytrain[closest_euc[i,j]]==levels[x])
        {
          euc_count[i,x]=euc_count[i,x] + (1/euc_dist[i,closest_euc[i,j]])
        } 
      }
    }
    y=as.list(which(euc_count[i,]==max(euc_count[i,])))
    pred_euc[i]=levels[as.numeric(y[1])]
    
  }
  
  ###############################################################
  closest_cb=apply(cb_dist,1,return_min_indices)
  closest_cb=matrix(data=NA, nrow = nrow(euc_dist), ncol = nrow(xtrain))
  for (i in 1:nrow(euc_dist))
  {
    temp=return_min_indices(euc_dist[i,])
    for (j in 1:length(temp))
    {
      closest_cb[i,j]=temp[j]
    }
  }
  cb_count=matrix(data=0,nrow(xtest),ncol=length(levels))
  pred_cb=rep(0,nrow(xtest))
  for (i in 1:nrow(closest_cb))
  {
    for (x in 1:length(levels))
    {
      for (j in 1:k)
      {
        if (ytrain[closest_cb[i,j]]==levels[x])
        {
          cb_count[i,x]=cb_count[i,x] + (1/cb_dist[i,closest_cb[i,j]])
        } 
      }
    }
    y=as.list(which(cb_count[i,]==max(cb_count[i,])))
    pred_cb[i]=levels[as.numeric(y[1])]
    
  }
  
  return(rbind(pred_euc,pred_cb))
}  
###################################################################

euc_errors=c()
cb_errors=c()
for (K in 2:6)
{
  errors_euc=c()
  errors_cb=c()
  for (i in 1:5)
  {
    train=read.csv(trainfilenames[i],header=TRUE)
    train=train[,-2]
    test=read.csv(testfilenames[i],header=TRUE)
    test=test[,-2]
    Ytrain=train[,ncol(train)]
    Xtrain=train[,-(ncol(train))]
    Ytest=test[,ncol(test)]
    Xtest=test[,-(ncol(test))]
    Xtrain=as.data.frame(lapply(Xtrain, as.numeric))
    Xtrain=scale(Xtrain)
    Xtest=as.data.frame(lapply(Xtest, as.numeric))
    Xtest=scale(Xtest)
    output=knn(Xtrain,Ytrain,Xtest,K,levels)
    misclassified_euc=length(which(output[1,]!=Ytest))
    misclassified_cb=length(which(output[2,]!=Ytest))
    error_euc=misclassified_euc/length(Ytest)
    errors_euc=c(errors_euc,error_euc)
    error_cb=misclassified_cb/length(Ytest)
    errors_cb=c(errors_cb,error_cb)
  }
  avg_error_euc=mean(errors_euc)
  euc_errors=c(euc_errors,avg_error_euc)
  avg_error_cb=mean(errors_cb)
  cb_errors=c(cb_errors,avg_error_cb)
}
par(mfrow=c(1,2))
plot(euc_errors~c(2:6),main="K vs Errors(Ionosphere) - Euclidean dist",ylab ="Classification error" ,xlab=" Value of K")
lines(euc_errors~c(2:6))
plot(cb_errors~c(2:6),main="K vs Errors(Ionosphere) - Cosine dist",ylab ="Classification error" ,xlab=" Value of K")
lines(cb_errors~c(2:6))