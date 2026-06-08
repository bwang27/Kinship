
check_pair<-function(fn_pair_matrix,ofn) {
  IDdata<-NULL
  for(i in 1:nrow(fn_pair_matrix)) {
    fnA<-fn_pair_matrix[i,1]
	fnB<-fn_pair_matrix[i,2]
    dt_A<-read.table(fnA,header=TRUE,sep="\t")
	dt_B<-read.table(fnB,header=TRUE,sep="\t")
    row.names(dt_A)<-apply(cbind(as.character(dt_A[,1]),as.character(dt_A[,2]),as.character(dt_A[,3]),as.character(dt_A[,4])),1,paste,collapse=".")
	row.names(dt_B)<-apply(cbind(as.character(dt_B[,1]),as.character(dt_B[,2]),as.character(dt_B[,3]),as.character(dt_B[,4])),1,paste,collapse=".")
    dt_A<-dt_A[dt_A[,"nT"]>20,]
    dt_B<-dt_B[dt_B[,"nT"]>20,]
    commonIdx<-intersect(row.names(dt_A),row.names(dt_B))
    theA<-dt_A[commonIdx,]
    theB<-dt_B[commonIdx,]
    A_carrier<-theA[,"nM"]>5 & theA[,"nM"]>0.1*theA[,"nT"]
    B_carrier<-theB[,"nM"]>5 & theB[,"nM"]>0.1*theB[,"nT"]
    ApBp<-sum(A_carrier & B_carrier)
    AnBp<-sum(!A_carrier & B_carrier)
    ApBn<-sum(A_carrier & !B_carrier)
    AnBn<-sum(!A_carrier & !B_carrier)
    tmp<-c(basename(fnA),basename(fnB),ApBp,ApBn,AnBp,AnBn)
	IDdata<-rbind(IDdata,tmp)
	print(i)
  }
  colnames(IDdata)<-c("SampleA","SampleB","Ap_Bp","Ap_Bn","An_Bp","An_Bn")
  write.table(IDdata,ofn,row.names=F,col.names=TRUE,sep="\t",quote=F)
}

args<-commandArgs(trailingOnly=TRUE)

if(length(args)==2) {
  #fn_pair_matrix<-read.table(args[1],header=FALSE,sep="\t",colClass="character")
  fn_pair_matrix<-read.table(args[1],header=TRUE,sep="\t",colClass="character")
  check_pair(fn_pair_matrix,args[2])
} else {
  print("Rscript ~ <filename_pair_matrix> <ofilename>")
  print("filename_pair_matrix has two columns: fn1 fn2")
  print("each of fn1/fn2 have 6 columns: Chr Pos Ref Alt nM nT")
}

