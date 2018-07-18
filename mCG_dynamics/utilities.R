library(fastcluster)

smooth_matrix <- function(x,smooth_factor=20,ncluster=50)
{
  x = na.omit(x)
  k = kmeans(x,centers=ncluster,iter.max=100)
  count = table(k$cluster)
  x_sm = NULL;
  for(c in 1:ncluster)
  {
    num = count[c]
    if(num == 1) #throw away small clusters
    {next}
    n = floor( num / smooth_factor ) 
    
    x_c = x[k$cluster==c,]
    h = hclust(dist(x_c),method="ward")
    #h = hclust(as.dist(1-cor(t(x_c))),method="ward")
    x_c = x_c[h$order,]
    
    if(n>0)
    {
      for(i in 1:n)
      {
        x_sm = rbind(x_sm, colMeans(x_c[((i-1)*smooth_factor):(i*smooth_factor-1),]))
      }      
    }
    
    if((num-(n)*smooth_factor) > 1)
    {
      x_sm = rbind(x_sm, colMeans(x_c[((n)*smooth_factor):(num),]))
    }
  }
  return(x_sm)
}


univ_entropy <- function(x) {
    p = x / sum(x,na.rm=T)
    p = p[p!=0]
    return( sum(- p * log(p)))
    
}

tissue_specificity_index <- function(x) {
    if( min(x) == 0  ) {
        return(0)
    }
    y = x / max(x)   
    return( sum(1-x) / (length(x) - 1))
}
