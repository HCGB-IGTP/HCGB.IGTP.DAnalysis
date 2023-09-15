##--------------------------
## Parallel code in R
##--------------------------
# https://www.blasbenito.com/post/02_parallelizing_loops_with_r/
# https://cran.r-project.org/web/packages/foreach/vignettes/foreach.html

parallel::detectCores()
n.cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


library(foreach)
library(doParallel)

x <- foreach(
  i = 1:10, 
  .combine = 'c'
) %dopar% {
  sqrt(i)
}
##--------------------------

