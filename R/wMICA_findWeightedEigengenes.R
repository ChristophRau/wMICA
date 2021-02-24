#Inputs
#expr: MxN matrix of expression, with M = samples and N = genes  
#MM: PxN matrix of module memberships, with N = genes and P = modules
#Outputs:
#WeightedEigengenes:  MxP matrix of Weighted first PCs, with M= samples and P = modules

findWeightedEigengenes <- function(expr,MM){
  library(ade4)
  tokeep=!is.na(MM[1,]) #Prune data such that only genes that are in the modules are used to calculate eigengenes
  MM=MM[,tokeep]
  expr=expr[,tokeep]
  
  datExpr=t(expr)
  MM=apply(MM,2,as.numeric)
  
  nSamples=ncol(datExpr)
  WMEs2=rep(0,nSamples)
  for(j in 1:nrow(MM)){   #Calculate weighted eigengenes
    print(j)
    outdata=dudi.pca(t(datExpr),col.w=as.numeric(MM[j,]),scannf=FALSE,nf=1)
    value=outdata$li
    WMEs2=cbind(WMEs2,value)
  }
  WMEs2=WMEs2[,-1]
  
  WeightedEigengenes=WMEs2
  return(WeightedEigengenes)
  
}