library(igraph)
x=100
temp=barabasi.game(x,power = 3)
network=get.edgelist(temp)
weights=c(runif(n=nrow(network),min = .4,max = 1))
network=cbind(network,weights)

L=network
Niter=10000
N=nrow(network)
M=x
Lindices=sample(N)
C=10
alpha=10
beta=.01
q=matrix(0,C,M)
ncores=4


z <- sample(1:C, N, replace=TRUE)
n <- as.vector(table(factor(z, levels=1:C)))
for (l in 1:N) {
  q[z[l], L[l,1]] <- q[z[l], L[l,1]]+L[l,3]
  q[z[l], L[l,2]] <- q[z[l], L[l,2]]+L[l,3]
}
conv <- vector("numeric",M)
q_orig=q


OrigMICA<-  function(L, Niter, N, M, Lindices,
           C, z, q, n, alpha, beta, conv) {
    ## Main iteration loop
    date()
    for (s in 1:Niter) {
      cat(".")
      ## Sample new component for each link
      for (li in 1:N) {

        l = Lindices[li]
        i = L[l,1]
        j = L[l,2]
        weight=L[l,3]

        ## Subtract the contribution of the link from the counts
        q[z[l],i] <- q[z[l],i]-weight
        q[z[l],j] <- q[z[l],j]-weight
        n[z[l]] <- n[z[l]]-1

        ## Loop for computing probabilities for the components to be sampled
        uz <- vector("numeric", C)
        for (p in 1:C) {
          A <- (n[p] +alpha)
          B <- (q[p,i] +beta)*(q[p,j] +beta) / ((2*n[p] +1 +M*beta)*(2*n[p] +M*beta))
          uz[p] <- A *B
        }

        ## Draw a new component for the links and update the counts */
        newz <- ICMg.multinom.single(uz)
        conv[l] <- uz[newz]/sum(uz)
        n[newz] <- n[newz]+1
        q[newz,i] <- q[newz,i]+weight
        q[newz,j] <- q[newz,j]+weight
        z[l] <- newz
      }
    }
  date()
    cat("\n")
    return(list(z=z, q=q, n=n, conv=conv))
  }


####MICA APPLY WORKHORSE FUNCTION####
  MAWorkhorse <-function(l){
    i = L[l,1]
    j = L[l,2]
    weight=L[l,3]
    q[z[l],i] <<- q[z[l],i]-weight
    q[z[l],j] <<- q[z[l],j]-weight
    n[z[l]] <<- n[z[l]]-1

    ## Loop for computing probabilities for the components to be sampled
    uz <- vector("numeric", C)
    for (p in 1:C) {
      A <- (n[p] +alpha)
      B <- (q[p,i] +beta)*(q[p,j] +beta) / ((2*n[p] +1 +M*beta)*(2*n[p] +M*beta))
      uz[p] <- A *B
    }
    newz <- ICMg.multinom.single(uz)
    conv[l] <<- uz[newz]/sum(uz)
    n[newz] <<- n[newz]+1
    q[newz,i] <<- q[newz,i]+weight
    q[newz,j] <<- q[newz,j]+weight
    z[l] <<- newz
  }

MICA_Apply<-  function(L, Niter, N, M, Lindices,
                     C, z, q, n, alpha, beta, conv) {
  ## Main iteration loop
  date()
  for (s in 1:Niter) {
    cat(".")
    ## Sample new component for each link
    temp=sapply(Lindices,MAWorkhorse)
  }
  date()

  cat("\n")
  return(list(z=z, q=q, n=n, conv=conv))
}

















MAPWorkhorse <-function(l){
  i = L[l,1]
  j = L[l,2]
  weight=L[l,3]
  q[z[l],i] <- q[z[l],i]-weight
  q[z[l],j] <- q[z[l],j]-weight
  n[z[l]] <- n[z[l]]-1
  q_new=array(0,dim=dim(q))
  n_new=rep(0,length(n))
  conv_new=rep(0,length(conv))
  z_new=rep(0,length(z))

  ## Loop for computing probabilities for the components to be sampled
  uz <- vector("numeric", C)
  for (p in 1:C) {
    A <- (n[p] +alpha)
    B <- (q[p,i] +beta)*(q[p,j] +beta) / ((2*n[p] +1 +M*beta)*(2*n[p] +M*beta))
    uz[p] <- A *B
  }
  newz <- ICMg.multinom.single(uz)
  conv_new[l] <- uz[newz]/sum(uz)
  n_new[newz] <- n[newz]+1
  q_new[newz,i] <- q[newz,i]+weight
  q_new[newz,j] <- q[newz,j]+weight
  z_new[l] <- newz
  return(list(n=n_new,q=q_new,z=z_new,conv=conv_new))
}




library(foreach)
library(doParallel)
MICA_Par<-  function(L, Niter, N, M, Lindices,
                       C, z, q, n, alpha, beta, conv,ncores) {
  ## Main iteration loop
  date()
  for (s in 1:Niter) {
    cat(".")
    ## Sample new component for each link
    c1<-makeCluster(ncores-1)
    clusterExport(c1,c("L","Niter","N","M","C","z","q","n","alpha","beta",
                       "conv","ICMg.multinom.single"))
    this_out=parLapply(c1,Lindices,MAPWorkhorse)
    stopCluster(c1)

    q=array(0,dim=dim(q))
    n=rep(0,length(n))
    conv=rep(0,length(conv))
    z=rep(0,length(z))

    for(i in 1:length(Lindices)){
      q=q_new+this_out[[i]]$q
      n=n_new+this_out[[i]]$n
      conv=conv_new+this_out[[i]]$conv
      z=z_new+this_out[[i]]$z
    }
  }




  cat("\n")
  return(list(z=z, q=q, n=n, conv=conv))
}














