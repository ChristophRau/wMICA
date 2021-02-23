wMICA_makeLinksFromMICFile <-function(MICname="MICA_MICARRAY.csv",output_prefix="MICA",threshold=.4){
  MINE=read.csv(MICname,row.names=1)
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
