#inputs section
ionosphere_df <- read.csv("ionosphere_data.csv", header=FALSE)#Input the data fromm the CSV file and store in Dataframe
input=ionosphere_df[,1:ncol(ionosphere_df)-1] # exclude the last column of the data as it is the true class atribute
K=2                                           # Change number of clusters K to the value needed
T=0.1                                         # Set the value of tolerance of change in centroid position
del=10                                            
itr=0                                         #iteration count variable
#End of inputs section
n=nrow(input)                                 #Find total number of observation
true_class= ifelse(ionosphere_df[,35]=='g',1,2)
#Function definition section begins here
#Function to assign random observations as the intial centroids
initialize_centroids <- function(c)                        
{
  rand_centroids=input[sample(nrow(input), c), ]
  return(rand_centroids)
}

euclidean_distance = function (a,b) # Distance of list a from Matrix b
{
  dist = sqrt(apply(t(apply(b,1,'-',a) )^2,1,sum))
  return (dist)
}


#Function to assign the closest cluster number to each observation
assign_cluster <- function(mat, centroids)
  
{
  
  cluster = c()
  dist=t(apply(input,1,euclidean_distance,as.matrix(centroids)))
  mindist=apply(dist,1,min)
  for (i in 1:n)
  {
    cluster[i]=which(dist[i,]==mindist[i])
  }
  return (cluster)
  
}

# Function to calculate the new position of the centroids
update_centroids <- function(mat, clust)
  
{
  
  cent <- c()
  
  for(c in 1:K)
  {
    group <- which(clust == c)
    # compute the mean point of all points in cluster c
    means <- colMeans(mat[group,])
    cent <- c(cent, means)         #Store the mean of each attribute for all centroids in a list       
    
  }
  
  cent <- matrix(cent, ncol=ncol(input),byrow=TRUE)  #Convert the list into a Matrix with same number of columns as the input
  
  return(cent)
  
}

#Function definition section ends here
#Main Program begins here
#Initalize the random centroids
centroids=initialize_centroids(K)

while(itr<25)                    #This is one stopping criteria : Upper bound on the number of iterations permitted
  
{
  
  clusters <- assign_cluster(input,centroids)
  
  if (length(unique(clusters))!=K)                    #To handle cluster collapse
  {
    missing=setdiff(c(1:K),clusters)
    for (m in missing)
    {
      centroids[m,]=input_data[sample(nrow(input_data),1), ]
    }
    clusters <- assign_cluster(input_data,centroids)
  }
  
  old_centroids=centroids
  
  centroids <- update_centroids(input, clusters)    #Update centroids
  
  error=t(apply(old_centroids,1,euclidean_distance,as.matrix(centroids))) 
  errors=diag(error) #To calculate the error - Shift in the position of the centroids between iterations
  del <- sum(errors)/K                        
   
  itr <- itr+1
  
  if(del<T) {break}            #Another stopping criteria : if the shift in centroids is less than tolerance value
  
}
Accuracy_percent = 100*max((length(which(clusters == true_class))/length(clusters)),(length(which(clusters != true_class))/length(clusters)))
centroids
clusters
itr
Accuracy_percent
