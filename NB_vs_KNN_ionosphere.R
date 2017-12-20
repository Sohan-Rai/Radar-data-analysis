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


#city_block_distance = function (a,b)
#{
#  dist = apply(abs(a-b),1,sum)
#  return (dist)
#}

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
  errors_euc=c()
  errors_cb=c()
  K=2
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
    error_euc=misclassified_euc/length(Ytest)
    errors_euc=c(errors_euc,error_euc)
  }
  avg_error_knn=mean(errors_euc)

############################################################################################

  library(mvtnorm)
  naiveBayes = function(train,test)
  {
    train=train[,-2]
    test=test[,-2]
    cat_cols=1
    y_col=ncol(train)
    level1='g'
    level2='b'
    levels=c(level1,level2)
    ##########################
    ytrain=train[,y_col]
    xtrain_cat=train[,cat_cols]
    xtrain=train[,-c(y_col,cat_cols)]
    ytest=test[,ncol(test)]
    xtest_cat=test[,cat_cols]
    xtest=test[,-c(y_col,cat_cols)]
    xtrain=as.data.frame(lapply(xtrain, as.numeric))
    #xtrain_cat=as.data.frame(lapply(xtrain_cat, as.factor))
    xtrain_cat=t(t(xtrain_cat))
    xtest_cat=t(t(xtest_cat))
    xtest=as.data.frame(lapply(xtest, as.numeric))
    m=length(ytrain)
    ########## Continuous ##########################
    gtrain=xtrain[which(train$V35==level1),]
    btrain=xtrain[which(train$V35==level2),]
    #################################################
    m1=nrow(gtrain)
    pc1=m1/m
    mu1=apply(gtrain,2,mean)
    gcentered=as.matrix(sweep(gtrain,2,mu1))
    var1=(1/m1)*t(gcentered)%*%gcentered
    ##################################################
    m2=nrow(btrain)
    pc2=m2/m
    mu2=apply(btrain,2,mean)
    bcentered=as.matrix(sweep(btrain,2,mu1))
    var2=(1/m2)*t(bcentered)%*%bcentered
    #################################################
    op_prob=matrix(data=NA,nrow=length(ytest),ncol=2)
    op_prob[,1]=dmvnorm(xtest, mean=mu1 ,sigma = var1)*pc1
    op_prob[,2]=dmvnorm(xtest, mean=mu2 ,sigma = var2)*pc2
    ###########  Categorical #########################
    match_count = function(xtrainc,ytrainc,case,class)
    {
      #match=c()
      count=0
      for (i in 1:length(xtrainc))
      {
        if(xtrainc[i]==case)
        {
          #match=c(match,i)
          if (ytrainc[i] == class)
          {
            count=count+1
          }
        }
      }
      return (count)
    }
    
    
    op_prob_cat=matrix(data=NA,nrow=length(ytest),ncol=2)
    p=2
    mc=3
    for (i in 1:nrow(xtest_cat))
    {
      probs1=1
      probs2=1
      for (j in 1:ncol(xtest_cat))
      { 
        nc1=match_count(xtrain_cat[,j],ytrain,xtest_cat[i,j],level1)
        n1=length(which(ytrain==level1))
        p1=1/p
        mc1=mc
        prob1=(nc1+mc1*p1)/(n1+mc1)
        probs1=probs1*prob1
        #####################
        nc2=match_count(xtrain_cat[,j],ytrain,xtest_cat[i,j],level2)
        n2=length(which(ytrain==level2))
        p2=1/p
        mc2=mc
        prob2=(nc2+mc2*p2)/(n2+mc2)
        probs2=probs2*prob2
      }
      op_prob_cat[i,1]=probs1
      op_prob_cat[i,2]=probs2
    }
    ##################################################
    op=matrix(data=NA,nrow=length(ytest),ncol=2)
    op[,1]=op_prob_cat[,1]*op_prob[,1]
    op[,2]=op_prob_cat[,2]*op_prob[,2]
    output=rep(NA,nrow(xtest))
    for (i in 1:nrow(xtest))
    {
      x=which(op[i,]==max(op[i,]))
      output[i]=levels[x]
    }
    output
    misclassified=length(which(output!=ytest))
    error=misclassified/length(ytest)
    return (error)
  }
  errs = rep(NA,5)
  for (k in 1:5)
  { 
    ip_train=read.csv(trainfilenames[k],header=TRUE)
    ip_test=read.csv(testfilenames[k],header=TRUE)
    errs[k]=naiveBayes(ip_train,ip_test)
  }
  avg_error_bayes=mean(errs)
  
############################################################################################  
combined_errors=c(errors_euc,errs)
par(mfrow=c(1,1))
plot(c(combined_errors,avg_error_knn,avg_error_bayes)~c(1:5,1:5,3,3),
     main="KNN vs Naive bayes (Ionosphere) ",
      ylab ="Classification error",
      xlab=" Validation set number",
     pch = c(rep(1,5),rep(4,5),15,15), 
     col = c("red","red","red","red","red", "blue",
             "blue", "blue", "blue", "blue","red","blue")
     )
#abline(avg_error_knn, 0,col="red",pch=15)
#abline(avg_error_bayes, 0,col="blue",pch=15)
legend("topright", pch = c(1,4,15,15), 
       col = c("red", "blue","red", "blue"), 
       legend = c("KNN", "Naive bayes","KNN avg error", "Naive bayes avg error"))
