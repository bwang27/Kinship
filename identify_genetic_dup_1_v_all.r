
check_pair_1_v_all<-function(one,all,ofn) {
  IDdata<-NULL
  dt_A<-read.table(one,header=TRUE,sep="\t")
  row.names(dt_A)<-apply(cbind(as.character(dt_A[,"Chr"]),as.character(dt_A[,"Pos"]),as.character(dt_A[,"Ref"]),as.character(dt_A[,"Alt"])),1,paste,collapse=".")
  dt_A<-dt_A[dt_A[,"nT"]>20,]
  for(i in 1:length(all)) {
    print(paste(c(i, length(all)),collapse=" / "))
    fnB<-all[i]
    dt_B<-read.table(fnB,header=TRUE,sep="\t")
    row.names(dt_B)<-apply(cbind(as.character(dt_B[,"Chr"]),as.character(dt_B[,"Pos"]),as.character(dt_B[,"Ref"]),as.character(dt_B[,"Alt"])),1,paste,collapse=".")
    dt_B<-dt_B[dt_B[,"nT"]>20,]
    commonIdx<-intersect(row.names(dt_A),row.names(dt_B))
    theA<-dt_A[commonIdx,]
    theB<-dt_B[commonIdx,]
    A_carrier<-theA[,"nM"]>=0.2*theA[,"nT"]
    B_carrier<-theB[,"nM"]>=0.2*theB[,"nT"]

    ApBp<-sum(A_carrier & B_carrier)
    AnBp<-sum(!A_carrier & B_carrier)
    ApBn<-sum(A_carrier & !B_carrier)
    AnBn<-sum(!A_carrier & !B_carrier)
    pcntID<-(ApBn+AnBp)/(ApBp+ApBn+AnBp+AnBn)
    pcntID<-as.integer((1000*pcntID))/10
    tmp<-c(basename(one),basename(fnB),ApBp,ApBn,AnBp,AnBn,pcntID)
    IDdata<-rbind(IDdata,tmp)
  }
  colnames(IDdata)<-c("SampleA","SampleB","Ap_Bp","Ap_Bn","An_Bp","An_Bn","pcntIdentity")
  write.table(IDdata,ofn,row.names=F,col.names=TRUE,sep="\t",quote=F)
}

#make_lists() 

args<-commandArgs(trailingOnly=TRUE)

if(length(args)==3) {
  all<-readLines(args[2])
  check_pair_1_v_all(args[1],all,args[3])
} else {
  print("Rscript ~ <wgs.one.samp.ifn> <all.list> <ofilename>")
  print("each of fn1/fn2 have 6 columns: Chr Pos Ref Alt nM nT")
}

