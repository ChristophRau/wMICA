
findRecommendedModuleNumber<-function(input,explained_threshold=.95,buffer=10){

if(buffer>1) buffer=buffer/100
total_explained=summary(prcomp(input))$importance["Cumulative Proportion",]

x_cross=length(total_explained[(total_explained<explained_threshold)])+1
if(x_cross>length(total_explained)) stop("Cumulative Variance Explained Never Crosses Threshold.  Select Lower Threshold.")
x_max=x_cross+ceiling(x_cross*buffer)+20
if(x_max>length(total_explained)) x_max=length(total_explained)

toplot=total_explained[0:x_max]

x_recommended=x_cross+ceiling(x_cross*buffer)


plot(toplot,ylab="Cumulative Variance Explained",xlab="Principle Component",
     main="Module Number Estimation",cex=1)
abline(v=x_cross,col="red",lty="dotted",lwd=3)
text(x=x_cross,y=min(toplot)+.01*min(toplot),"Minimum",pos=2,col="red")

abline(v=x_recommended,col="blue",lty="dotted",lwd=3)
text(x=x_recommended,y=min(toplot)+.01*min(toplot),"Recommended",pos=4,col="blue")

abline(h=explained_threshold,col="black",lty="dotted",lwd=1)
text(x=2,y=explained_threshold,labels = "Threshold",pos=3)

message(paste0(x_cross," modules crosses threshold"))
message(paste0(x_recommended," modules are recommended"))
return(x_recommended)
}
