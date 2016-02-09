cl   <- function(dat,fm, cluster){
              attach(dat, warn.conflicts = F)
              library(sandwich)
              M <- length(unique(cluster))   
              N <- length(cluster)  	        
              K <- fm$rank		             
              dfc <- (M/(M-1))*((N-1)/(N-K))  
              uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
              vcovCL <- sandwich(fm, meat=crossprod(uj)/N)
              #coeftest(fm, vcovCL) 
}
