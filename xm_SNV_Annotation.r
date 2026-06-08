
prepare_gene_granges <- function(refseq_df, promoter_flank) {
  gr_tx <- GenomicRanges::GRanges(seqnames = refseq_df$chr,
                                  ranges = IRanges::IRanges(refseq_df$txS, refseq_df$txE),
                                  strand = refseq_df$strand,gene = refseq_df$gName,nmID = refseq_df$nmID,
                                  cdsS = refseq_df$cdsS,cdsE = refseq_df$cdsE)
  ps <- ifelse(refseq_df$strand == "+",refseq_df$txS - promoter_flank - 1,refseq_df$txE + 1)
  pe <- ifelse(refseq_df$strand == "+",refseq_df$txS - 1,refseq_df$txE + promoter_flank + 1)
  ps <- pmax(ps, 1)
  pe <- pmax(pe, ps)
  gr_prom <- GenomicRanges::GRanges(seqnames = refseq_df$chr,ranges = IRanges::IRanges(ps, pe),gene = refseq_df$gName,nmID = refseq_df$nmID)
  exonS <- strsplit(refseq_df$exonS, ",", fixed=TRUE)
  exonE <- strsplit(refseq_df$exonE, ",", fixed=TRUE)
  exonS <- lapply(exonS, function(x) as.integer(x[nzchar(x)]))
  exonE <- lapply(exonE, function(x) as.integer(x[nzchar(x)]))
  n <- lengths(exonS)
  gr_exon <- GenomicRanges::GRanges(seqnames = rep(refseq_df$chr, n),
                                    ranges = IRanges::IRanges(start = pmax(1, unlist(exonS) - 20), end = pmax(unlist(exonS) - 20,unlist(exonE) + 20)),
                                    gene = rep(refseq_df$gName, n), nmID = rep(refseq_df$nmID, n), cdsS = rep(refseq_df$cdsS, n), cdsE = rep(refseq_df$cdsE, n))
  ret<-list(tx = gr_tx, prom = gr_prom, exon = gr_exon)
  return(ret)
}


obtain_centro_telo_one_build <- function(build, ofn) {
  chr_keep <- paste0("chr", c(1:22, "X", "Y","M","MT"))
  rlt <- data.frame(Chr=character(), Start=character(),End=character(), repName=character(),repClass=character(),repFamily=character(), stringsAsFactors = FALSE)
  if (build == "hg38") {
    cmd <- "wget -q -O centromeres.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz && gunzip -f centromeres.txt.gz"
    if (system(cmd) != 0) stop("Failed to download centromeres.txt.gz")
    dt <- read.table("centromeres.txt", header = FALSE, sep = "\t",colClasses = "character")
    cen <- data.frame(Chr = dt[,2],Start = dt[,3],End = dt[,4],repName = "centromere", repClass="centromere", repFamily="centromere", stringsAsFactors = FALSE)
    cen <- cen[cen$Chr %in% chr_keep, ]
    rlt <- rbind(rlt, cen)
    unlink("centromeres.txt")
  }
  myurl <- paste0("https://hgdownload.soe.ucsc.edu/goldenPath/",build, "/database/gap.txt.gz")
  cmd <- paste0("wget -q -O gap.txt.gz ", shQuote(myurl)," && gunzip -f gap.txt.gz")
  if (system(cmd) != 0) stop("Failed to download gap.txt.gz")
  dt <- read.table("gap.txt", header = FALSE, sep = "\t",colClasses = "character",quote = "", comment.char = "")
  keep_types <- c("telomere", "short_arm")
  if(build != "hg38") keep_types<-c(keep_types,"centromere")
  ndt <- dt[dt[,8] %in% keep_types & dt[,2] %in% chr_keep, ]
  gaps <- data.frame(Chr = ndt[,2],Start = ndt[,3],End = ndt[,4],repName = ndt[,8],repClass=ndt[,8], repFamily=ndt[,8], stringsAsFactors = FALSE)  
  rlt <- rbind(rlt, gaps)
  write.table(rlt, ofn, row.names = FALSE,col.names = TRUE, sep = "\t", quote = FALSE)
  unlink("gap.txt")
}
obtain_masked_repeats_one_build<-function(build,ofn) {
  obtain_centro_telo_one_build(build, ofn)
  chr_keep <- paste0("chr", c(1:22, "X", "Y","M","MT"))
  rlt<-read.table(ofn,header=TRUE,sep="\t",colClasses="character")  
  myurl <- paste0("https://hgdownload.gi.ucsc.edu/goldenPath/",build,"/database/rmsk.txt.gz")
  cmd <- paste0("wget -q -O rmsk.txt.gz ", shQuote(myurl)," && gunzip -f rmsk.txt.gz")
  if (system(cmd) != 0) stop("Failed to download or gunzip rmsk.txt.gz")
  dt <- read.table("rmsk.txt",header = FALSE,sep = "\t",colClasses = "character",quote = "",comment.char = "")
  ndt<-dt[dt[,6] %in% chr_keep,]  
  repeats <- data.frame(Chr = ndt[,6],Start = ndt[,7],End = ndt[,8],repName = ndt[,11],repClass=ndt[,12], repFamily=ndt[,13], stringsAsFactors = FALSE)
  rlt<-rbind(rlt,repeats)
  rlt<-unique(rlt)  
  write.table(rlt,ofn,row.names=F,col.names=TRUE,sep="\t",quote=F)
  unlink("rmsk.txt")
}
prepare_repeat_masks_include_telo_centro<-function() {
  print("doing hg19")
  obtain_masked_repeats_one_build("hg19","repeat_telo_centro.hg19.bed")
  print("doing hg19")
  obtain_masked_repeats_one_build("hg38","repeat_telo_centro.hg38.bed")
}
#prepare_repeat_masks_include_telo_centro()

get_repeats_SNPs_gene_definition<-function(build,file_loc,promoter_flank) {
  repeat_hg19<-file_loc$repeat_hg19; repeat_hg38<-file_loc$repeat_hg38; refseq_hg19<-file_loc$refseq_hg19; refseq_hg38<-file_loc$refseq_hg38
  chr_keep <- paste0("chr", c(1:22, "X", "Y","M","MT"))
  if(!(build %in% c("hg19","hg38"))) stop("build must be in hg19/hg38")
  repeat_df_fn<- if(build=="hg19") repeat_hg19 else repeat_hg38
  print("reading repeat/telo/centro file")
  repeat_df<-read.table(repeat_df_fn,header=TRUE,sep="\t",stringsAsFactors=FALSE) #header is Chr, Start, End, repName, repClass, repFamily  
  print("done")
  required <- c("Chr","Start","End")
  if (!all(required %in% colnames(repeat_df))) stop("repeat_df missing required columns: Chr, Start, End")  
  repeat_df$Start<-as.integer(repeat_df$Start)
  repeat_df$End<-as.integer(repeat_df$End)
  if(sum(repeat_df$End<repeat_df$Start)) stop("repeat file has entries where End < Start")
  repeat_df$Chr <- ifelse(grepl("^chr", repeat_df$Chr),repeat_df$Chr,paste0("chr", repeat_df$Chr))
  repeat_df<-repeat_df[!is.na(repeat_df$Start) & !is.na(repeat_df$End),]
  repeat_df<-repeat_df[repeat_df$Chr %in% chr_keep,]  
  repeat_gr <- GenomicRanges::GRanges(seqnames = repeat_df$Chr,ranges   = IRanges::IRanges(start=repeat_df$Start,end=repeat_df$End))
  refseq_df_fn<- if(build=="hg19") refseq_hg19 else refseq_hg38
  refseq_df<-read.table(refseq_df_fn,header=TRUE,sep="\t",stringsAsFactors=FALSE)
  refseq_df<-refseq_df[refseq_df$isDefault=="Yes",]
  refseq_df$chr<-as.character(refseq_df$chr)  
  refseq_df$chr <- ifelse(grepl("^chr", refseq_df$chr),refseq_df$chr,paste0("chr", refseq_df$chr))  
  refseq_df$txS<-as.integer(refseq_df$txS)+1
  refseq_df$txE<-as.integer(refseq_df$txE)  
  refseq_df$txS_m <- refseq_df$txS - 10000
  refseq_df$txE_m <- refseq_df$txE + 10000
  refseq_df$txS_m[refseq_df$txS_m<=1]<-1
  refseq_df$cdsS<-as.integer(refseq_df$cdsS)+1
  refseq_df$cdsE<-as.integer(refseq_df$cdsE)
  refseq_df<-refseq_df[refseq_df$chr %in% chr_keep,]
  refseq_df$exonS <- vapply(strsplit(refseq_df$exonS, ",", fixed=TRUE),
                            function(x) { x <- x[nzchar(x)]; xi <- suppressWarnings(as.integer(x)); 
                                          if (anyNA(xi)) stop("Malformed exonS entry detected");
                                          paste0(xi + 1, collapse=",")},character(1))
  gene_gr <- prepare_gene_granges(refseq_df, promoter_flank)
  rlt<-list(repeat_df = repeat_df,refseq_df=refseq_df, gene_gr=gene_gr, repeat_gr=repeat_gr)
  return(rlt)
}


further_annotate_gene_region_df <- function(dt, chr_column, pos_s_column, pos_e_column, gene_gr) {
  chr <- dt[,chr_column]; pos_s <- as.integer(dt[,pos_s_column]); pos_e <- as.integer(dt[,pos_e_column]);
  if(any(is.na(chr)) | any(is.na(pos_s)) | any(is.na(pos_e))) stop("dt has na")
  gr_var <- GenomicRanges::GRanges(chr, IRanges(pos_s, pos_e))
  n <- nrow(dt); RefSeq <- GeneName <- Region <- rep("NoNoNoNoNoData",n);
  add_val<-function(x,i,v) {  x[i]<-paste(x[i],v,sep=";"); x } 
  tx_mc   <- mcols(gene_gr$tx); exon_mc <- mcols(gene_gr$exon); prom_mc <- mcols(gene_gr$prom)
  ov <- findOverlaps(gr_var, gene_gr$exon)  #want all transcript-level exon hits
  if (length(ov)) {
    q <- queryHits(ov); s <- subjectHits(ov);
    cdsS <- exon_mc$cdsS[s]-10; cdsE <- exon_mc$cdsE[s]+5;
    coding <- pos_s[q] >= pmin(cdsS, cdsE) & pos_e[q] <= pmax(cdsS, cdsE)
    reg <- ifelse(coding, "coding", "UTR")
    Region <- add_val(Region,q,reg); GeneName <- add_val(GeneName,q,exon_mc$gene[s]); RefSeq <- add_val(RefSeq,q,exon_mc$nmID[s]);
  }
  ov <- findOverlaps(gr_var, gene_gr$prom)  #want all transcript-level exon hits
  if (length(ov)) {
    q <- queryHits(ov); s <- subjectHits(ov);
    Region<-add_val(Region,q,rep("promoter",length(q))); GeneName<-add_val(GeneName,q,prom_mc$gene[s]); RefSeq<-add_val(RefSeq,q,prom_mc$nmID[s]); 
  }
  is_not_annotated <- Region == "NoNoNoNoNoData"
  if (any(is_not_annotated)) {
    gr_sub <- gr_var[is_not_annotated]
    ov <- findOverlaps(gr_sub, gene_gr$tx)
    if (length(ov)) {
      q_sub <- queryHits(ov); q <- which(is_not_annotated)[q_sub]  # map back to original indices
      s <- subjectHits(ov)
      Region<-add_val(Region, q, "intron"); GeneName<-add_val(GeneName,q,tx_mc$gene[s]); RefSeq<-add_val(RefSeq,q,tx_mc$nmID[s]);
    }
  }
  Region <- gsub("NoNoNoNoNoData", "T3",gsub("NoNoNoNoNoData;","",Region))
  GeneName <- gsub("NoNoNoNoNoData", "T3",gsub("NoNoNoNoNoData;","",GeneName))
  RefSeq <- gsub("NoNoNoNoNoData", "T3",gsub("NoNoNoNoNoData;","",RefSeq))
  dt$Region <- Region; dt$GeneName <- GeneName; dt$RefSeq <- RefSeq;
  return(dt)
}

tier_snvs_indels_per_sample <- function(markerfile, chr_column, pos_s_column,pos_e_column, marker_column, build, snp_repeat_gene_df,grep_out_file) {
  if (!(build %in% c("hg19","hg38"))) stop("build must be hg19/hg38")
  if (!requireNamespace("GenomicRanges", quietly=TRUE) || !requireNamespace("IRanges", quietly=TRUE)) stop("Packages GenomicRanges and IRanges are required")
  if (!requireNamespace("S4Vectors", quietly=TRUE)) stop("Package S4Vectors is required")

  snp_df<-read.table(grep_out_file,header=FALSE,sep="\t",colClasses="character")
  colnames(snp_df)<-c("Marker_hg38","Marker_hg19","nM_gnomAD","nT_gnomAD")
  snp_column <- if(build=="hg38") 1 else 2
  row.names(snp_df)<-snp_df[,snp_column]
  dt <- read.table(markerfile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  dt[,chr_column]<- ifelse(grepl("^chr", dt[,chr_column]),dt[,chr_column],paste0("chr", dt[,chr_column]))  
  
  chr_keep <- paste0("chr", c(1:22, "X", "Y", "M", "MT"))
  dt <- dt[dt[, chr_column] %in% chr_keep, ]
  dt[, pos_s_column] <- as.integer(dt[, pos_s_column])
  dt[, pos_e_column] <- as.integer(dt[, pos_e_column])
  dt <- dt[!is.na(dt[, pos_s_column]) & !is.na(dt[, pos_e_column]), ]
  ofn_SNP  <- paste0(markerfile, "_in_SNP.txt")
  ofn_rep <- paste0(markerfile, "_in_Repeat.txt")
  ofn_to_Process<-paste0(markerfile, "_to_Process.txt")
  ofn_T3<-paste0(markerfile, "_T3.txt")
  if (nrow(dt) == 0) {
    write.table(data.frame(Result="Nothing"),ofn_SNP,row.names=FALSE, sep="\t", quote=FALSE)
    write.table(data.frame(Result="Nothing"),ofn_rep,row.names=FALSE, sep="\t", quote=FALSE)
    write.table(data.frame(Result="Nothing"),ofn_to_Process,row.names=FALSE, sep="\t", quote=FALSE)
    return(invisible(dt))
  }
  repeat_gr <- snp_repeat_gene_df$repeat_gr
  gene_gr <- snp_repeat_gene_df$gene_gr
  gr_vars <- GenomicRanges::GRanges(seqnames = dt[, chr_column],ranges   = IRanges::IRanges(start=dt[, pos_s_column],end=dt[, pos_e_column]))
  ov   <- GenomicRanges::findOverlaps(gr_vars, repeat_gr)
  hits <- S4Vectors::queryHits(ov)  
  is_inRepeat <- logical(nrow(dt))
  is_inRepeat[unique(hits)] <- TRUE
  is_inSNP <- dt[,marker_column] %in% row.names(snp_df)
  is_droppedSNP <- is_inSNP 
  if (any(is_droppedSNP)) {
    inSNP <- dt[is_droppedSNP, ]
    row.names(inSNP) <- inSNP[, marker_column]
    inSNP$nM_gnomAD <- snp_df[row.names(inSNP), "nM_gnomAD"]
    inSNP$nT_gnomAD <- snp_df[row.names(inSNP), "nT_gnomAD"]
    write.table(inSNP, ofn_SNP, row.names=FALSE, sep="\t", quote=FALSE)
  } else {
    write.table(data.frame(Result="Nothing"),ofn_SNP, row.names=FALSE, sep="\t", quote=FALSE)
  }
  is_droppedRepeat <- is_inRepeat & !is_inSNP 
  if (any(is_droppedRepeat)) {
    write.table(dt[is_droppedRepeat, ],ofn_rep, row.names=FALSE, sep="\t", quote=FALSE)
  } else {
    write.table(data.frame(Result="Nothing"),ofn_rep, row.names=FALSE, sep="\t", quote=FALSE)
  }
  is_remaining <- !is_droppedSNP & !is_droppedRepeat
  if (any(is_remaining) == 0) {
    write.table(data.frame(Result="Nothing"),ofn_to_Process,row.names=FALSE, sep="\t", quote=FALSE)
  } else {    
    print("obtaining intron/UTR/coding/promoter regions")
    remaining<-dt[is_remaining,]
    processed<-further_annotate_gene_region_df(remaining, chr_column, pos_s_column, pos_e_column, gene_gr)
    write.table(processed[!(processed$Region %in% c("T3","intron")),],ofn_to_Process,row.names=FALSE, sep="\t", quote=FALSE)
    write.table(processed[processed$Region %in% c("T3","intron"),],ofn_T3,row.names=FALSE, sep="\t", quote=FALSE)
  }
  return(invisible(dt))
}

extract_SNVP_from_one_high20<-function(h20, ofn,ofn_4_grep) {
  chr_keep <- paste0("chr", c(1:22, "X", "Y","M","MT"))
  print("reading high20")
  dt<-read.table(h20,header=TRUE,sep="\t",colClasses="character")
  print("done")
  ndt<-dt[dt$Type=="SNP" & dt$Chr %in% chr_keep,]
  markers <- paste(ndt$Chr,ndt$Pos,ndt$Chr_Allele,ndt$Alternative_Allele,sep=".")
  Pos_E <- as.integer(ndt$Pos) + nchar(ndt$Chr_Allele)
  rlt<-cbind(ndt[,c("Chr","Pos")],Pos_E=Pos_E,Marker=markers,Ref=ndt$Chr_Allele,Mut=ndt$Alternative_Allele)

  colnames(rlt)<-c("Chr","Pos_S","Pos_E","Marker","Ref","Alt")
  write.table(rlt, ofn, row.names=F, col.names=TRUE,sep="\t",quote=F)
  write.table(markers, ofn_4_grep, row.names=F, col.names=FALSE, sep="\t",quote=F)
}

process_SNPV_from_many_high20<-function(fns_hg20,build,promoter_flank,file_loc) {
  fns<-read.table(fns_hg20,header=TRUE,sep="\t",colClasses="character")
  if(!(nrow(fns)>0 & ncol(fns)==2)) stop("no data to process")  
  snp_repeat_gene_df<-get_repeats_SNPs_gene_definition(build,file_loc,promoter_flank)
  polymorphism_file <- file_loc$polymorphism_file
  chr_column<-1
  pos_s_column<-2
  pos_e_column<-3
  marker_column<-4
  for(i in 1:nrow(fns)) {
    h20<-fns[i,2]
    prefix<-fns[i,1]
    markerfile<-paste(c(prefix,"_step0"),collapse="")
    grep_file <- paste(c(prefix,"_step0_4_grep"),collapse="")
    grep_out_file <- paste(c(prefix,"_step0_grepout"),collapse="")
    extract_SNVP_from_one_high20(h20, markerfile, grep_file)
    cmd <- sprintf("fgrep -f %s %s > %s",shQuote(grep_file),shQuote(polymorphism_file),shQuote(grep_out_file))
    rc <- system(cmd)
    tier_snvs_indels_per_sample(markerfile, chr_column, pos_s_column, pos_e_column, marker_column, build, snp_repeat_gene_df,grep_out_file)
    unlink(grep_file)
    unlink(grep_out_file)
    unlink(markerfile)
  }
}

suppressPackageStartupMessages(library(GenomicRanges))

polymorphism_file<-"/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/refs/gnomAD_v3.1.2_hg38_genome/pcnt0.01_hg38hg19.snv"
repeat_hg19<-"/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/refs/repeat_telo_centro.hg19.bed"
repeat_hg38<-"/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/refs/repeat_telo_centro.hg38.bed"
refseq_hg19<-"/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/refs/ncbiRefSeqCurated_hg19_default.txt"
refseq_hg38<-"/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/refs/ncbiRefSeqCurated_hg38_default.txt"

file_loc<-list(polymorphism_file=polymorphism_file,repeat_hg19=repeat_hg19,repeat_hg38=repeat_hg38,refseq_hg19=refseq_hg19,refseq_hg38=refseq_hg38)


do_SNV_anno_by_args_from_high20<-function(args) {
  if(length(args)==5) {
    fns_hg20<-args[3]
    build<-args[4]
    promoter_flank<-as.integer(args[5])
    process_SNPV_from_many_high20(fns_hg20,build,promoter_flank,file_loc)
  } else {
    print("Rscript ~ <ann> <snv> <fns_hg20> <build> <promoter_flank>")  
  }
}

args<-commandArgs(trailingOnly=TRUE)
do_SNV_anno_by_args_from_high20(args)

