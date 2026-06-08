

prep_marker_file<-function() {
  print("reading data")
  mkrs<-read.table("gnomad.genomes.v4.1.sites.filtered_maf_1_to_99.tsv",header=TRUE,sep="\t")  
  print("splitting SNV")
  snvs<-matrix(unlist(strsplit(mkrs[,1],split=".",fixed=TRUE)),byrow=TRUE,ncol=4)
  print("making chr.pos")
  chr.pos<-apply(snvs[,1:2],1,paste,collapse=".")
  ref<-snvs[,3]
  alt<-snvs[,4]
  nmkrs<-cbind(mkrs,chr.pos,ref,alt)
  print("writting")
  write.table(nmkrs,"gnomad.genomes.v4.1.sites.filtered_maf_1_to_99_mod.tsv",row.names=F,sep="\t",quote=F)
}

do_calling<-function(ifn,ofn,depth_cut=15) {
  mkrs<-read.table("gnomad.genomes.v4.1.sites.filtered_maf_1_to_99_mod.tsv",header=TRUE,sep="\t",colClass="character")
  dt<-read.table(ifn,header=F,sep="\t")
  row.names(dt)<-apply(cbind(dt[,1],as.character(dt[,2])),1,paste,collapse=".")
  colnames(dt)<-c("chr","pos","A","C","G","T","N")

  ndt<-NULL
  lets<-c("A","C","G","T")

  for(ref in 1:4) {
    for(alt in setdiff(1:4,ref)) {
      tmpmkr<-mkrs[mkrs[,"ref"]==lets[ref] & mkrs[,"alt"]==lets[alt],]
      tmpdt<-dt[row.names(dt) %in% tmpmkr[,"chr.pos"],]
      SNP<-apply(cbind(row.names(tmpdt),lets[ref],lets[alt]),1,paste,collapse=".")
      geno<-rep(-1,nrow(tmpdt))
      nM<-tmpdt[,lets[alt]]
      AF<-tmpdt[,lets[alt]]/(0.01+tmpdt[,"N"])
      AF<-as.integer(100*AF)
      nT<-tmpdt[,"N"]
      geno[AF<=20 & nT>=depth_cut]<-0
      geno[AF>=80 & nT>=depth_cut]<-2
      geno[AF>20 & AF<80 & nT>=depth_cut]<-1
      tmpdt<-cbind(SNP,nM,nT,geno)
      ndt<-rbind(tmpdt,ndt)
    }
  }
  write.table(ndt,ofn,row.names=F,col.names=TRUE,sep="\t",quote=F)
}

#prep_marker_file()
args<-commandArgs(trailingOnly=TRUE)
if(length(args)==2) {
  do_calling(args[1],args[2],depth_cut=15)
} else {
  print("Rscript ~ <ifn> <ofn>")
  print("ifn: count file; 7 columns: chr pos A C G T N")
}


