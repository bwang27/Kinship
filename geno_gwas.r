

#rm jobs.sh; touch jobs.sh
#while read line
#do
#  cat sample.lst | awk -v c="${line}" '{print "grep -e Ref.Alt -e " c " " $0 " > " c "." $0}'  >> jobs.sh
#done < mychrs
#split -l100 jobs.sh; rm jobs.sh
#ls x* | awk '{print "bsub -q compbio -M8000 -eo eo.txt -oo oo.txt sh " $1}' | sh
##rm x??




split_markers<-function(n_Size) {
  dt<-read.table("gnomad.genomes.v4.1.sites.filtered_maf_1_to_99_mod.tsv",header=TRUE,sep="\t",colClass="character")
  chrs<-matrix(unlist(strsplit(dt[,"Chr.Pos"],split=".",fixed=TRUE)),ncol=2,byrow=TRUE)[,1]
  mychrs<-unique(chrs)
  ofns<-NULL
  NN<-1
  for(i in 1:length(mychrs)) {
    print(mychrs[i])
    ndt<-dt[chrs==mychrs[i],]
    starts<-seq(from=1,to=nrow(ndt),by=n_Size)
	if(length(starts)>1) {
      ends<-c(starts[2:length(starts)]-1,nrow(ndt))
	} else {
	  ends<-nrow(ndt)
	}
    for(j in 1:length(starts)) {
      ofn<-paste(c("part_",NN,".mkr"),collapse="")
	  NN<-NN+1
	  write.table(ndt[starts[j]:ends[j],],ofn,row.names=F,col.names=TRUE,sep="\t",quote=F)  
	  ofns<-c(ofns,ofn)
	  print(ofn)
    }
  }
  return(ofns)
}

merge_one_part<-function(mkr_part) {
  smplst_ifn<-"sample.lst"
  mkrs<-read.table(mkr_part,header=TRUE,sep="\t",colClass="character")
  row.names(mkrs)<-mkrs[,1]
  chr<-unique(matrix(unlist(strsplit(mkrs[,"Chr.Pos"],split=".",fixed=TRUE)),ncol=2,byrow=TRUE)[,1])
  chr<-paste(chr,".S",sep="")
  smps<-grep(chr,readLines(smplst_ifn),value=TRUE)
  for(i in 1:length(smps)) {
    SJID<-unlist(strsplit(smps[i],split=".",fixed=TRUE))[2]
	print(paste(c(i, length(smps)), collapse=" / "))
	print(SJID)	
	print("reading data")
	system("date")
    tmp<-read.table(smps[i],header=TRUE,sep="\t",colClass="character")
	print("data read finished")
	system("date")
	row.names(tmp)<-apply(cbind(tmp[,"Chr"],as.character(tmp[,"Pos"]),tmp[,"Ref.Alt"]),1,paste,collapse=".")
	system("date")
	tmp<-tmp[row.names(tmp) %in% row.names(mkrs),c("nM","nT","geno")]
	other_mkr<-setdiff(row.names(mkrs),row.names(tmp))
	other<-matrix(-1,nrow=length(other_mkr),ncol=3)
	row.names(other)<-other_mkr
	colnames(other)<-colnames(tmp)
	tmp<-rbind(tmp,other)[row.names(mkrs),]
	colnames(tmp)<-apply(cbind(colnames(tmp),SJID),1,paste,collapse=".")
	mkrs<-cbind(mkrs,tmp)
  }
  ofn<-paste(c(mkr_part,"_w_geno"),collapse="")
  write.table(mkrs,ofn,row.names=F,col.names=TRUE,sep="\t",quote=F)
}


#rm jobs.sh
#touch jobs.sh
#aiplot="/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/script/ai_check/ai_plot_cnvsv.r"
#chrsize="/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/script/ai_check/chr_sizes_hg19.txt"
#ls SJ* | grep _S > s.lst
#ls SJ* | grep _E > e.lst
#while read line
#do
#  sid=`echo $line | cut -f1 -d'_'`
#  smp=`echo $line | cut -f1 -d'.'`
#  theE=`grep $sid e.lst`
#  etag=`echo $theE | cut -f1 -d'.' | cut -f2 -d'_' `
#  echo "bsub -q compbio -M10000 -eo eo.txt -oo oo.txt Rscript $aiplot $chrsize ${line}+${theE} D1N1 NULL ${smp}_${etag}.jpg ${smp}" >> jobs.sh
#done < s.lst
#rm e.lst s.lst


p_fisher<-function(nAA,nAa,naa) {
  PP<-0-10*log10(0.9)
  nCase<-nAA+nAa+naa
  if(nCase>300) {
    nA<-2*nAA+nAa
    nN<-2*(nAA+nAa+naa)
    pA<-nA/(nN+0.00000001)
    exp_nAA<-as.integer(sprintf("%.0f",nCase*pA*pA))
    exp_nAa<-as.integer(sprintf("%.0f",nCase*pA*(1-pA)*2))
    exp_naa<-as.integer(sprintf("%.0f",nCase*(1-pA)*(1-pA)))
	mat<-matrix(c(nAA,nAa,naa,exp_nAA,exp_nAa,exp_naa),ncol=3,byrow=TRUE)
	#print(mat)
	uu<-abs(mat[1,]-mat[2,])
	MM<-apply(mat,2,max)
	vv<-max(uu/(MM+0.1))
	#print(uu)
	#print(vv)
	if(max(uu)>=0.03*nCase & vv>0.08) {
	  PP<-fisher.test(mat,simulate.p.value=FALSE)$p.value
	  PP<-0-10*log10(PP)
	}  
	PP<-as.integer(PP*10)/10
  } else {
    PP<-0-9
  }  
  return(PP)
}

consolidate_by_case<-function(part_geno) {
  dt<-read.table(part_geno,header=TRUE,sep="\t")
  info<-dt[,1:7]
  samples<-grep("geno.S",colnames(dt),value=TRUE)
  nMM<-grep("nM.S",colnames(dt),value=TRUE)
  nTT<-grep("nT.S",colnames(dt),value=TRUE)
  USI<-unique(matrix(unlist(strsplit(gsub("geno.S","S",samples),split="_")),ncol=2,byrow=TRUE)[,1])
  rlt<-NULL  
  for(i in 1:length(USI)) {
    if(i%%100==0) print(i)
    theG<-grep(paste(USI[i],"_G",sep=""),samples,value=TRUE)
	theD<-setdiff(grep(USI[i],samples,value=TRUE),theG)
	nMs<-gsub("geno.S","nM.S",c(theG,theD))
	nTs<-gsub("geno.S","nT.S",c(theG,theD))
	if(length(nMs)>1) {
	  sumM<-apply(dt[,nMs],1,sum)
	  sumT<-apply(dt[,nTs],1,sum)
	  sumGeno<-rep(-1,nrow(dt))
	  sumGeno[sumT>=15 & sumM>=0.8*sumT]<-2
	  sumGeno[sumT>=15 & sumM<0.8*sumT & sumM>0.2*sumT]<-1
	  sumGeno[sumT>=15 & sumM<=0.2*sumT]<-0
	  sumGeno[sumT>=10 & sumM>=0.95*sumT]<-2
	  sumGeno[sumT>=10 & sumM<=0.05*sumT]<-0
	}
	genocall<-rep(-1,nrow(dt))
	if(length(theG)==1) {
	  genocall<-dt[,theG]
	}	
	if(length(theD)>0) {
	  for(j in 1:length(theD)) {
	    call_na<-which(genocall<0)
		genocall[call_na]<-dt[call_na,theD[j]]		
	  }
	}
	if(length(nMs)>0) {
	  call_na<-which(genocall<0)
	  genocall[call_na]<-sumGeno[call_na]	
	}
	rlt<-cbind(rlt,genocall)  
  }
  colnames(rlt)<-USI
  rownames(rlt)<-rownames(info)
  nCalled<-apply(rlt>0-0.5,1,sum)
  nTwo<-apply(rlt==2,1,sum)
  nOne<-apply(rlt==1,1,sum)
  nZero<-apply(rlt==0,1,sum)
  nPos<-nTwo*2+nOne*1
  nTot<-(nTwo+nOne+nZero)*2
  PhredP_HWE<-rep(-1,nrow(dt))
  dAF<-sprintf("%.4f",nPos/(0.001+nTot))  
  for(i in 1:nrow(dt)) {
    PhredP_HWE[i]<-p_fisher(nTwo[i],nOne[i],nZero[i])
  }
  ret<-cbind(info,nCalled,nTwo,nOne,nZero,nPos,nTot,dAF,PhredP_HWE,rlt)
  ofn<-paste(part_geno,"_consolidated",sep="")
  write.table(ret,ofn,row.names=F,col.names=TRUE,sep="\t",quote=F)
}

#parts<-split_markers(n_Size=50000)

#args<-commandArgs(trailingOnly=TRUE)
#if(length(args)==1) {
#  merge_one_part(args[1])
#} else {
#  print("Rscript ~ <marker.part>")
#}

args<-commandArgs(trailingOnly=TRUE)
if(length(args)==1) {
  consolidate_by_case(args[1])
} else {
  print("Rscript ~ <part_geno>")
}
