wMICA_RunMIC <- function(input,output_prefix="MICA",n.cores=detectCores()-1,threshold=.4){
  input=apply(input,2,as.numeric)
  print(paste0("Starting MINERVA with ",n.cores, " cores.  Please Wait."))
  MINE=mine(x=input,n.cores=n.cores)$MIC
  outname=paste0(output_prefix,"_MICARRAY.csv")
  print(paste0("saving MIC Array as ",outname))
  write.csv(MINE,file=outname)
  outname=paste0(output_prefix,"_Links.csv")
  print(paste0("Creating Links File ",outname))
  outfile=file(outname,"w")
  firstRow="Gene1,Gene2,Strength"
  cat(firstRow,file=outfile,sep="\n")
  MINE_THRESH=MINE>threshold
  edges=0
  for(i in 1:(nrow(MINE)-1)){
    for(j in (i+1):nrow(MINE)){
      if(MINE_THRESH[i,j]){
        edges=edges+1
        outrow=paste(i,j,MINE[i,j],sep=",")
        cat(outrow,file=outfile,sep="\n")
      }
    }
  }
  close(outfile)

}
