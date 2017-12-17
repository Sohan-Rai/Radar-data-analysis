#inputs section
library(MASS)
library(mvtnorm)
library(ggplot2)
K=2  ### Set the number of clusters needed
T=0.2
del=10
inner_itr=0
itr=1
df <- read.csv("ionosphere_data.csv", header=FALSE)#Input the data fromm the CSV file and store in Dataframe
input_data=df[,1:ncol(df)-1]                  #Retaining only the Numeric columns of the dataframe
input_data=data.matrix(input_data, rownames.force = NA)
input_original=input_data
zero_cols=apply(input_data,2,sum)
zero_cols=which(zero_cols==0)
if (length(zero_cols)>0)
{
  input_data=input_data[,-zero_cols]
}
input_data=apply(input_data, 2, function(x) {(x - mean(x))/sd(x)})  ### Mean normalize the data.
input_quality= df[,ncol(df)]
input_quality=t(t(input_quality))
#End of inputs section
n=nrow(input_data)
cols=ncol(input_data)

#######  Random initialization of cluster centroids  ##############
initialize_centroids <- function()
{
  rand_centroids=input_data[1:K, ]
  return(rand_centroids)
}

####### Initialization of mean and variance #################
  mus=initialize_centroids()
  sig=lapply(1:K, function(x) diag(cols))
  prob_mat=matrix(0,nrow=n,ncol=K)
  pc=c()
  for (i in 1:K)
  {
    pc=c(pc,1/K)
  }
  denom=c()
  dellogs=c()
  itr_errors=c()
  cluster_div=c()
  inner_itr=0
while(inner_itr<100)
  {
    #############  Expectatioin step  ##########################
    prob_mat_last=prob_mat
      for (j in 1:K)
      {
        prob_mat[,j]=dmvnorm(input_data, mean=mus[j,],sigma = sig[[j]])*pc[j]
      }
    prob_mat=prob_mat/rowSums(prob_mat)
    collapse_check=colSums(prob_mat)
    clusters=c()
    for (i in 1:n)
    {
      clusters[i]=which.max(prob_mat[i,])
    }
    
    if(0 %in% collapse_check | NaN %in% collapse_check) 
    { 
      break;
    }
    errors=c()
    wcolsum=c()
    sig_nums=c()
    old_mus=mus
    #############  Maximization step  #########################
    for (j in 1:K)
    {
      wcolsum[j]=sum(prob_mat[,j])
      pc[j]=wcolsum[j]/n
      mus_num=(prob_mat[,j]*input_data)
      mus[j,]=colSums(mus_num)/wcolsum[j]
      sig_new=matrix(0,ncol=cols,nrow=cols)
      for (c in 1:cols)
      {
        sig_nums=prob_mat[,j]*(input_data[,c]-mus[j,c])*(input_data[,c]-mus[j,c])
        sigcc=sum(sig_nums)/wcolsum[j]
        sig_new[c,c]=ifelse(sigcc==0,0.001,sigcc)   #### To handle singularity
      }
      sig[[j]]=sig_new
    }
    errors=c()
    #### Cluster centroid shift based stopping criteria #########
    for (i in 1:K)
    {
      errors[i] <- sqrt(sum((old_mus[i,]-mus[i,])**2)) # calculates the difference b for each K
    }
    del <- sum(errors)/K
    itr_errors=c(itr_errors,del)
    inner_itr=inner_itr+1
    if (del<0.00001){break}
    #### Log likelyhood based stopping criteria #########
    lkhs=c()
    log_lkhs=c()
    for (i in 1:n)
    {
      for (j in 1:K)
      {
        lkhs[j]=dmvnorm(input_data[i,], mean=mus[j,],sigma = sig[[j]])*pc[j]
      }
      log_lkhs[i]=log(sum(lkhs))
    }
    dellog=sum(log_lkhs)
    dellogs=c(dellogs,dellog)
    if (inner_itr>1)
    {
      if (dellogs[length(dellogs)] < dellogs[length(dellogs)-1]) {break;}
    }
}
  ########### Interpret the hard clusters from the soft clustering performed ###########
  clusters=ifelse(prob_mat[,1]>prob_mat[,2],1,2)
  cluster_div=c(cluster_div,length(which(clusters==1)),length(which(clusters==2)))
  cluster_div
  true_clusters=ifelse(input_quality=="g",1,2)
  true_cluster_div=c(length(which(true_clusters==1)),length(which(true_clusters==2)))
  true_cluster_div
  Accuracy_percent = 100*max((length(which(clusters == true_clusters))/length(clusters)),(length(which(clusters != true_clusters))/length(clusters)))
  Accuracy_percent
