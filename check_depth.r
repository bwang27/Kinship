
depth_of_one_sample<-function(fn) {
  dt<-read.table(fn,header=TRUE)
  med<-as.integer(10*median(dt[,"nT"]))/10
  avg<-as.integer(10*mean(dt[,"nT"]))/10
  n20<-sum(dt[,"nT"]>20)
  return(c(med,avg,n20))
}

args<-commandArgs(trailingOnly=TRUE)

fns<-readLines(args[1])
project_id<-args[2]

#fns<-grep("ai_check_",dir(path=".",pattern=".txt$"),value=TRUE)
smps<-gsub("ai_check_|.txt","",basename(fns))
info<-cbind(Sample=smps,File=fns)

rlt<-NULL
for(i in 1:length(fns)) {
  rlt<-rbind(rlt, depth_of_one_sample(fns[i]))
  print(fns[i])
}

rlt<-cbind(smps,rlt)
colnames(rlt)<-c("Sample","Median","Mean","nMkr_20X")
write.table(rlt, paste(project_id,"_depth_stats.txt",sep=""),row.names=F,col.names=TRUE,quote=F,sep="\t")

kk<-NULL
for(i in 1:(length(fns)-1)) {
  for(j in (i+1):length(fns)) {
    kk<-rbind(kk,c(fns[i],fns[j]))
  }
}
write.table(kk, paste(project_id,"_all_pairs.txt",sep=""),row.names=F,col.names=F,sep="\t",quote=F)

