
check_pair_1_v_all<-function(one,all,ofn) {
  IDdata<-NULL
  dt_A<-read.table(one,header=TRUE,sep="\t")
  row.names(dt_A)<-apply(cbind(as.character(dt_A[,1]),as.character(dt_A[,2]),as.character(dt_A[,3]),as.character(dt_A[,4])),1,paste,collapse=".")
  dt_A<-dt_A[dt_A[,"nT"]>5,]
  for(i in 1:length(all)) {
    fnB<-all[i]
    dt_B<-read.table(fnB,header=TRUE,sep="\t")
    row.names(dt_B)<-apply(cbind(as.character(dt_B[,1]),as.character(dt_B[,2]),as.character(dt_B[,3]),as.character(dt_B[,4])),1,paste,collapse=".")
    dt_B<-dt_B[dt_B[,"nT"]>5,]
    commonIdx<-intersect(row.names(dt_A),row.names(dt_B))
    theA<-dt_A[commonIdx,]
    theB<-dt_B[commonIdx,]
    A_carrier<-theA[,"nM"]>1 & theA[,"nM"]>0.1*theA[,"nT"]
    B_carrier<-theB[,"nM"]>1 & theB[,"nM"]>0.1*theB[,"nT"]
    ApBp<-sum(A_carrier & B_carrier)
    AnBp<-sum(!A_carrier & B_carrier)
    ApBn<-sum(A_carrier & !B_carrier)
    AnBn<-sum(!A_carrier & !B_carrier)
    tmp<-c(basename(one),basename(fnB),ApBp,ApBn,AnBp,AnBn)
    IDdata<-rbind(IDdata,tmp)
  }
  colnames(IDdata)<-c("SampleA","SampleB","Ap_Bp","Ap_Bn","An_Bp","An_Bn")
  write.table(IDdata,ofn,row.names=F,col.names=TRUE,sep="\t",quote=F)
}

#make_lists() 

args<-commandArgs(trailingOnly=TRUE)

if(length(args)==3) {
  #fn_pair_matrix<-read.table(args[1],header=FALSE,sep="\t",colClass="character")
  all<-readLines(args[2])
  check_pair_1_v_all(args[1],all,args[3])
} else {
  print("Rscript ~ <wgs.one.samp.ifn> <all.list> <ofilename>")
  print("each of fn1/fn2 have 6 columns: Chr Pos Ref Alt nM nT")
}

