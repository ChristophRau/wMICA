#Inputs
#input:  Rows are Samples, Columns are Genes
#output_prefix: prefix to affix to outputted files.
#n.cores:  number of cores.  Defaults to total number minus 1  (detectCores())
#threshold:  lowest MINE value to send to ICMg.  Defaults to .3.  Lower equals 'better' but takes a LOT longer
#C:  Number of Modules in ICMg
#Seeds:  Which genes to seed.  Names need to be the same as in the input file
#Weight:  Weight to give each Seed.
#alpha:  for ICMg
#beta:  for ICMg
#B.num:  Number of Burnin rounds
#B.size:  Size of one burnin round
#S.num:  Number of Sample Rounds
#S.size:  Size of one sample round

'wMICA.Run'<-function(input,output_prefix="MICA",n.cores=detectCores()-1,threshold=.4,C=20,Seeds=c(),Weight=0,
                      alpha=10,beta=.01,B.num=10,B.size=10,S.num=10,S.size=10,links_file=c()){

  if(length(links_file)==0){
  wMICA_GetMIC(input,output_prefix,n.cores,threshold)
    links_file=paste0(output_prefix,"_Links.csv")}

  gene_names=colnames(input)
  if(length(Seeds)>0){
    temp=match(Seeds[,1],gene_names)
    Seeds[,1]=temp
    Seeds=apply(Seeds,2,as.numeric)
  }

LINKS=read.csv(links_file)
LINKS=apply(LINKS,2,as.numeric)
nodes=c(LINKS[,1],LINKS[,2])
nodes=length(table(nodes))
edges=nrow(LINKS)

print(paste0("There are ",nodes," genes and ", edges," edges.  Starting Weighted ICMg."))
wICMg=ICMg.links.sampler.SeedsAndWeights(L = LINKS,C = C,Seeds = Seeds,Weight = Weight,
                                         alpha = alpha,beta = beta,B.num = B.num,
                                         B.size = B.size,S.num = S.num,S.size = S.size,
                                         C.boost = 1)
print("Calculating Memberships...")
wICMg$comp.memb=ICMg.get.comp.memberships(LINKS,wICMg)
outdata=wICMg$comp.memb
rownames(outdata)=c(1:C)
colnames(outdata)=gene_names
outname=paste0(output_prefix,"_MODULEMEMBERSHIPS.csv")
print(paste0("Writing Output File ",outname))
write.csv(outdata,file=outname)
print("Done!")
}

