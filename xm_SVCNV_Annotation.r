is_chr <- function(x) {
  valid <- c(as.character(1:22), "X", "Y", "M", "MT")
  as.character(x) %in% c(paste0("chr", valid),"?")
}

check_sv_type <- function(sv) {
  is_sv  <- FALSE
  type   <- dists <- "NoNoNoNoData"
  if (is.na(sv)) return(list(is_sv=is_sv, SV_Type=type, SV_Dist=dists))
  vals <- unlist(strsplit(sv, ".", fixed=TRUE))
  if (length(vals) < 6 || length(vals) > 7) return(list(is_sv=is_sv, SV_Type=type, SV_Dist=dists))
  cA <- vals[1]; cB <- vals[4]
  pA <- suppressWarnings(as.integer(vals[2]))
  pB <- suppressWarnings(as.integer(vals[5]))
  oA <- vals[3]; oB <- vals[6]
  #if (is.na(pA) || is.na(pB)) return(list(is_sv=is_sv, SV_Type=type, SV_Dist=dists))
  if (xor(cA=="?", cB=="?"))  return(list(is_sv=TRUE, SV_Type="side_1", SV_Dist=999999)) #cA and cB will never be coded as chr?
  if (is.na(pA) || is.na(pB)) return(list(is_sv=is_sv, SV_Type=type, SV_Dist=dists))
  if (!is_chr(cA) || !is_chr(cB)) return(list(is_sv=is_sv, SV_Type=type, SV_Dist=dists))
  if (!(oA %in% c("+","-") && oB %in% c("+","-"))) return(list(is_sv=is_sv, SV_Type=type, SV_Dist=dists))
  is_sv <- TRUE
  if (cA == cB) {
    dists <- abs(pB - pA)
    if (oA == oB) {
      if (oA == "+") {
        type <- if (pB > pA) "DEL" else "INS"
      } else {
        type <- if (pB > pA) "INS" else "DEL"
      }
    } else {
      type <- "INV"
    }
  } else {
    type <- "CTX"
    dists <- 999999
  }
  list(is_sv=is_sv, SV_Type=type, SV_Dist=dists)
}

check_cnv_type <- function(cnv) {
  is_cnv  <- FALSE
  dists <- "NoNoNoNoData"
  if (is.na(cnv)) return(list(is_cnv=is_cnv, CNV_Dist=dists))
  vals <- unlist(strsplit(cnv, ".", fixed=TRUE))
  if (length(vals) != 3) return(list(is_cnv=is_cnv, CNV_Dist=dists))
  
  cA <- vals[1]
  pA <- suppressWarnings(as.integer(vals[2]))
  pB <- suppressWarnings(as.integer(vals[3]))
  if (!is_chr(cA) || is.na(pA) || is.na(pB) || cA=="?") return(list(is_cnv=is_cnv, CNV_Dist=dists))
  dists<-abs(1+pB-pA)
  is_cnv <- TRUE  
  list(is_cnv=is_cnv, CNV_Dist=dists, pA=min(pA,pB), pB=max(pA,pB),chr=cA)
}


get_sv_anno_one_refseq <- function(gGeneName, gRefSeq, gStrand, gnExon,gcdsS, gcdsE, gexonS, gexonE, gframeS,pos, Promoter_Flank) {
  GeneName <- as.character(gGeneName)
  RefSeq   <- as.character(gRefSeq)
  Strand   <- as.character(gStrand)
  ExonIntron <- ExonFrame <- ExonIntronPos <- WithinUTR <- "NoNoNoNoData"
  if (length(grep("NM", RefSeq, value=TRUE)) == 0) {
    return(list(GeneName = GeneName,RefSeq = RefSeq, Strand = Strand, ExonIntron = " ", ExonFrame = " ", ExonIntronPos = " ", WithinUTR = " "))
  }
  exonS  <- as.integer(strsplit(gexonS, ",")[[1]])
  exonE  <- as.integer(strsplit(gexonE, ",")[[1]])
  frameS <- as.integer(strsplit(gframeS, ",")[[1]])
  if (length(exonS) != gnExon) stop("Exon length mismatch")
  o <- order(exonS, exonE)
  exonS <- exonS[o]; exonE <- exonE[o]; frameS <- frameS[o]
  pos <- as.integer(pos)
  if (!is.na(gcdsS) && !is.na(gcdsE) && gcdsS > 0 && gcdsE > 0) {
    if (pos >= gcdsS && pos <= gcdsE) {
      WithinUTR <- "coding"
    } else if (pos < gcdsS) {
      WithinUTR <- if (Strand=="+") "UTR5" else "UTR3"
    } else {
      WithinUTR <- if (Strand=="+") "UTR3" else "UTR5"
    }
  }
  firstExonStart <- exonS[1]
  lastExonEnd <- exonE[length(exonE)]
  if ((Strand=="+" && pos >= firstExonStart-Promoter_Flank && pos < firstExonStart) ||
      (Strand=="-" && pos > lastExonEnd && pos <= lastExonEnd+Promoter_Flank)) {
    ExonIntron <- "promoter"
  }
  hit_exon <- which(exonS <= pos & pos <= exonE)
  if (length(hit_exon)) {
    i <- hit_exon[1]
    dist5 <- if (Strand=="+") pos-exonS[i] else exonE[i]-pos
    dist3 <- if (Strand=="+") exonE[i]-pos else pos-exonS[i]
    ExonIntronPos <- paste(dist5, dist3, sep="_")
    exonNum <- if (Strand=="-") gnExon+1-i else i
    ExonIntron <- paste0("exon_", exonNum)
    ExonFrame <- as.character(frameS[i])
  } else if (gnExon > 1) {
    intronS <- exonE[-length(exonE)] + 1
    intronE <- exonS[-1] - 1
    hit_intron <- which(intronS <= pos & pos <= intronE)
    if (length(hit_intron)) {
      j <- hit_intron[1]
      dist5 <- if (Strand=="+") pos-intronS[j] else intronE[j]-pos
      dist3 <- if (Strand=="+") intronE[j]-pos else pos-intronS[j]
      ExonIntronPos <- paste(dist5, dist3, sep="_")
      intronNum <- if (Strand=="-") (gnExon-1)-(j-1) else j
      ExonIntron <- paste0("intron_", intronNum)
    }
  }
  list(GeneName=GeneName, RefSeq=RefSeq, Strand=Strand,
       ExonIntron=ExonIntron, ExonFrame=ExonFrame,
       ExonIntronPos=ExonIntronPos, WithinUTR=WithinUTR)
}

get_sv_anno_one_side <- function(def_by_chr, chr, pos, Promoter_Flank) {
  ret<-list(GeneName="NoNoNoNoData", RefSeq="NoNoNoNoData",Strand="NoNoNoNoData", 
            ExonIntron="NoNoNoNoData",ExonFrame="NoNoNoNoData", ExonIntronPos="NoNoNoNoData",WithinUTR="NoNoNoNoData")
  if (is.null(def_by_chr)) return(ret)
  chr_def <- def_by_chr[[chr]]
  if (is.null(chr_def)) return(ret)
  hit <- which(chr_def$txS_i - Promoter_Flank < pos & pos < chr_def$txE_i + Promoter_Flank)
  if (!length(hit)) return(ret)
  res <- lapply(hit, function(i) {
    get_sv_anno_one_refseq(chr_def$gName[i], chr_def$nmID[i], chr_def$strand[i],
      chr_def$exonN_i[i], chr_def$cdsS_i[i], chr_def$cdsE_i[i],
      chr_def$exonS[i], chr_def$exonE[i], chr_def$exonFrames[i],
      pos, Promoter_Flank)})

  Reduce(function(a,b) {
    lapply(names(a), function(n) { gsub("NoNoNoNoData;","",paste(a[[n]],b[[n]],sep=";")) }) |> setNames(names(a))
  }, res)
}

find_neighbor_genes_one_side <- function(def_by_chr, chr, pos, hasGenes, nNeighborGeneToInclude) {
  if (is.null(def_by_chr)) return("")
  hasGene <- if (hasGenes=="" || hasGenes=="NoNoNoNoData") character(0)
             else strsplit(hasGenes,";")[[1]]
  chr_def <- def_by_chr[[chr]]
  if (is.null(chr_def)) return("")

  df <- chr_def[!(chr_def$gName %in% hasGene),]
  left  <- df[df$txE_i < pos,]
  right <- df[df$txS_i > pos,]
  out <- character(0)
  if (nrow(left)) {
    o <- order(pos - left$txE_i)
    for (i in o[seq_len(min(nNeighborGeneToInclude,length(o)))]) {
      out <- c(out, paste0(left$gName[i],"[",ceiling((pos-left$txE_i[i])/1000),"k]"))
    }
  }
  if (nrow(right)) {
    o <- order(right$txS_i - pos)
    for (i in o[seq_len(min(nNeighborGeneToInclude,length(o)))]) {
      out <- c(out, paste0(right$gName[i],"[",ceiling((right$txS_i[i]-pos)/1000),"k]"))
    }
  }
  paste(out, collapse=";")
}

get_sv_anno_two_sides <- function(def_by_chr, sv, Promoter_Flank, nNeighborGeneToInclude) {
  qc <- check_sv_type(sv)
  if (!qc$is_sv) return(paste(rep(" ",18), collapse="\t"))
  vals <- strsplit(sv,"\\.")[[1]]
  chrA <- vals[1]; chrB <- vals[4]
  resA <- resB <- rep(" ",7)
  moreA <- moreB <- " "
  if (chrA!="?") {
    posA <- as.integer(vals[2])
    rA <- get_sv_anno_one_side(def_by_chr, chrA, posA, Promoter_Flank)
    resA <- unlist(rA)
    moreA <- find_neighbor_genes_one_side(def_by_chr, chrA, posA, rA$GeneName, nNeighborGeneToInclude)
  }
  if (chrB!="?") {
    posB <- as.integer(vals[5])
    rB <- get_sv_anno_one_side(def_by_chr, chrB, posB, Promoter_Flank)
    resB <- unlist(rB)
    moreB <- find_neighbor_genes_one_side(def_by_chr, chrB, posB, rB$GeneName, nNeighborGeneToInclude)
  }
  paste(c(qc$SV_Type, qc$SV_Dist, resA, resB, moreA, moreB), collapse="\t")
}

annotate_many_svs <- function(genedef_ifn, sv_ifn, sv_column, ofn, Promoter_Flank, nNeighborGeneToInclude) {
  genedef <- read.table(genedef_ifn, header=TRUE, sep="\t", colClasses="character")
  chrs <- unique(genedef$chr)
  bad  <- chrs[!grepl("^chr(\\d+|X|Y|M|MT)$", chrs)]
  if (length(bad) > 0) {
    stop(paste0("Check gene definition ", genedef_ifn," for chromosome names; Bad values: ", paste(bad, collapse=", ")))
  }
  genedef$txS_i <- as.integer(genedef$txS)
  genedef$txE_i <- as.integer(genedef$txE)
  genedef$cdsS_i <- as.integer(genedef$cdsS)
  genedef$cdsE_i <- as.integer(genedef$cdsE)
  genedef$exonN_i <- as.integer(genedef$exonN)
  default<-genedef[genedef$isDefault=="Yes",]
  def_by_chr <- split(default, default$chr)
  sv_data<-read.table(sv_ifn,header=TRUE,sep="\t",colClasses="character")
  uniq_sv <- unique(sv_data[,sv_column])
  good_svs<- NULL
  for(i in 1:length(uniq_sv)) {
    qcSV<-check_sv_type(uniq_sv[i])
    if(qcSV$is_sv) good_svs<-c(good_svs,uniq_sv[i])
  }
  rlt<-rep("",length(good_svs))
  default_anno<-paste(rep(" ",18),collapse="\t")
  names(rlt)<-good_svs
  for (i in 1:length(good_svs)) {
    rlt[i] <- gsub("NoNoNoNoData"," ",get_sv_anno_two_sides(def_by_chr, good_svs[i],Promoter_Flank, nNeighborGeneToInclude))
    if (i %% 500 == 0) message(i," / ",length(uniq_sv))
  }
  names_newcol <- paste(c("SV_Type","SV_Dist",
            "GeneName_A","RefSeq_A","Strand_A","ExonIntron_A","ExonFrame_A","ExonIntronPos_A","WithinUTR_A",
            "GeneName_B","RefSeq_B","Strand_B","ExonIntron_B","ExonFrame_B","ExonIntronPos_B","WithinUTR_B",
            "MoreGene_A","MoreGene_B"), collapse="\t")
  new_col <- rep(default_anno, nrow(sv_data))
  m <- match(sv_data[, sv_column], names(rlt))
  ok <- !is.na(m)
  new_col[ok] <- rlt[m[ok]]
  final<-cbind(sv_data,new_col) 
  colnames(final)[ncol(final)]<-names_newcol
  write.table(final, ofn, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
}


get_cnv_anno_two_sides <- function(def_by_chr, cnv, CNV_GeneFlank) {
  GeneName_Inside <- GeneName_B <- GeneName_A <- dists <- " "
  qc <- check_cnv_type(cnv)
  if (!qc$is_cnv) return(paste(rep(" ", 4), collapse="\t"))
  chrA <- qc$chr;   pA   <- qc$pA;   pB   <- qc$pB;   dists <- qc$CNV_Dist
  if (!is_chr(chrA) || chrA == "?") return(paste(rep(" ", 4), collapse="\t"))
  def_dt <- def_by_chr[[chrA]]
  if (is.null(def_dt) || nrow(def_dt) == 0) return(paste(rep(" ", 4), collapse="\t"))
  span_A_gene <- def_dt$txS_i - CNV_GeneFlank <= pA & def_dt$txE_i + CNV_GeneFlank >= pA
  span_B_gene <- def_dt$txS_i - CNV_GeneFlank <= pB & def_dt$txE_i + CNV_GeneFlank >= pB
  inside_Gene <- def_dt$txE_i >= pA & def_dt$txS_i <= pB
  if (any(span_A_gene, na.rm=TRUE))   GeneName_A <- paste(unique(def_dt$gName[span_A_gene]), collapse=";")
  if (any(span_B_gene, na.rm=TRUE))   GeneName_B <- paste(unique(def_dt$gName[span_B_gene]), collapse=";")
  if (any(inside_Gene, na.rm=TRUE))   {
    already_has<-union(def_dt$gName[span_A_gene], def_dt$gName[span_B_gene])
    myGenes<-unique(setdiff(def_dt$gName[inside_Gene],already_has))
	if(length(myGenes)>4) {
	  GeneName_Inside<-paste(c(paste(myGenes[1:4],collapse=";")," + ", length(myGenes)-4, " more"),collapse= "")	  
	} else {
	  GeneName_Inside<-paste(myGenes,collapse=";")
	}    
  }
  paste(c(dists, GeneName_A, GeneName_Inside, GeneName_B), collapse="\t")
}

annotate_many_cnvs <- function(genedef_ifn, cnv_ifn, cnv_column, ofn,CNV_GeneFlank) {
  genedef <- read.table(genedef_ifn, header=TRUE, sep="\t", colClasses="character")
  chrs <- unique(genedef$chr)
  bad  <- chrs[!grepl("^chr(\\d+|X|Y|M|MT)$", chrs)]
  if (length(bad) > 0) {
    stop(paste0("Check gene definition ", genedef_ifn," for chromosome names; Bad values: ", paste(bad, collapse=", ")))
  }
  genedef$txS_i <- as.integer(genedef$txS)
  genedef$txE_i <- as.integer(genedef$txE)
  genedef$cdsS_i <- as.integer(genedef$cdsS)
  genedef$cdsE_i <- as.integer(genedef$cdsE)
  genedef$exonN_i <- as.integer(genedef$exonN)
  default<-genedef[genedef$isDefault=="Yes",]
  def_by_chr <- split(default, default$chr)
  cnv_data<-read.table(cnv_ifn,header=TRUE,sep="\t",colClasses="character")
  uniq_cnv <- unique(cnv_data[,cnv_column])
  good_cnvs<- NULL
  for(i in 1:length(uniq_cnv)) {
    qccnv<-check_cnv_type(uniq_cnv[i])
    if(qccnv$is_cnv) good_cnvs<-c(good_cnvs,uniq_cnv[i])
  }
  rlt<-rep("",length(good_cnvs))
  default_anno<-paste(rep(" ",4),collapse="\t")
  names(rlt)<-good_cnvs
  for (i in 1:length(good_cnvs)) {
    rlt[i] <- get_cnv_anno_two_sides(def_by_chr, good_cnvs[i],CNV_GeneFlank)
    if (i %% 500 == 0) message(i," / ",length(uniq_cnv))
  }
  names_newcol <- paste(c("cnv_Dist","GeneName_A","GeneName_Inside","GeneName_B"), collapse="\t")
  new_col <- rep(default_anno, nrow(cnv_data))
  m <- match(cnv_data[, cnv_column], names(rlt))
  ok <- !is.na(m)
  new_col[ok] <- rlt[m[ok]]
  final<-cbind(cnv_data,new_col) 
  colnames(final)[ncol(final)]<-names_newcol
  write.table(final, ofn, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
}


genedef_ifn_hg38 <- "/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/refs/ncbiRefSeqCurated_hg38_default.txt"
genedef_ifn_hg19 <- "/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/refs/ncbiRefSeqCurated_hg19_default.txt"
Promoter_Flank <- 2000
CNV_GeneFlank  <- 5000
nNeighborGeneToInclude <- 2

do_sv_anno_by_args<-function(args) {
  if(length(args)>=6) {
    sv_ifn <- args[3]; sv_column <- as.integer(args[4]); build <- args[5]; ofn <- args[6]
	Promoter_Flank_local <- if(length(args)==8) as.integer(args[7]) else Promoter_Flank
	nNeighborGeneToInclude_local <- if(length(args)==8) as.integer(args[8]) else nNeighborGeneToInclude
	genedef_ifn <- if (build=="hg38") genedef_ifn_hg38 else genedef_ifn_hg19
	annotate_many_svs(genedef_ifn, sv_ifn, sv_column, ofn, Promoter_Flank_local, nNeighborGeneToInclude_local)
  } else {
	print("Rscript ~ <anno> <sv> <sv_ifn> <sv_column> <build|hg19|hg38> <ofn> | <promoter_flank> <n_neighbors> ")
	print("          sv_ifn: must have header; if not, first row will be skipped")
	print("          promoter_flank: 2000 default, n_neighbors: 2 nearest default; sv_column: 1-based")
  }
}

do_cnv_anno_by_args<-function(args) {
  if(length(args)>=6) {
    cnv_ifn <- args[3]; cnv_column <- as.integer(args[4]); build <- args[5]; ofn <- args[6]
	genedef_ifn <- if (build=="hg38") genedef_ifn_hg38 else genedef_ifn_hg19
	CNV_GeneFlank_local <- if(length(args)==7) as.integer(args[7]) else CNV_GeneFlank
	annotate_many_cnvs(genedef_ifn, cnv_ifn, cnv_column, ofn, CNV_GeneFlank_local)
  } else {
	print("Rscript ~ <anno> <cnv> <cnv_ifn> <cnv_column> <build|hg19|hg38> <ofn> | <cnv_geneflank>")
	print("          cnv_ifn: must have header; if not, first row will be skipped")
	print("          cnv_column: 1-based")
  }
}

parse_args<-function(args) {
  if(toupper(args[1])=="ANNO") {    
    if(length(args)>=2) {	  
      if(toupper(args[2])=="SV") do_sv_anno_by_args(args)
      if(toupper(args[2])=="CNV") do_cnv_anno_by_args(args)
	} else {
	  print("Rscript ~ <anno> <cnv|sv>")
	}
  }
}

args <- commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  print("Rscript ~ <anno>")
} else {
  parse_args(args)
}


