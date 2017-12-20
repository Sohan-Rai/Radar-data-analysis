library(mvtnorm)
trainfilenames=c()
testfilenames=c()
for (i in 1:5)
{
  trainfilename=paste("ionosphere_train_", i, ".csv",sep="")
  trainfilenames=c(trainfilenames,trainfilename)
  testfilename=paste("ionosphere_test_", i, ".csv",sep = "")
  testfilenames=c(testfilenames,testfilename)
}
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
avg_err=mean(errs)
plot(errs~c(1:5),main="Error for each validation set - Ionosphere",ylab ="Classification error" ,xlab="Validation set number")
lines(errs~c(1:5))
avg_err
