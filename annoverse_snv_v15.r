###############################################################################
## initialization 
###############################################################################

initialize_constants <- function(build) {
  stopifnot(build %in% c("hg19","hg38"))
  ## ---------------- Tool paths ----------------
  Tools <- list(samtools = "/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/apps/samtools")
  stopifnot(file.exists(Tools$samtools))
  stopifnot(file.access(Tools$samtools, 1) == 0)  # executable

  ## ---------------- Motif builder ----------------
  build_motif <- function(counts,flank5_len,core_len,flank3_len,score_cutoff,delta_cutoff,key) {
    possible_BASES <- c("A","C","G","T")
    genome_freq <- c(A=0.3, C=0.2, G=0.2, T=0.3)
    genome_freq <- genome_freq / sum(genome_freq)
    if (!all(possible_BASES %in% names(genome_freq))) stop("genome_freq missing entries for: ", paste(setdiff(possible_BASES, names(genome_freq)), collapse = ","))
    if (is.null(rownames(counts)) || !all(possible_BASES %in% rownames(counts))) stop("counts must have rownames: ", paste(possible_BASES, collapse = ","))
    counts <- counts[possible_BASES, , drop = FALSE]
    stopifnot(ncol(counts) == flank5_len + core_len + flank3_len)
    pseudo <- 0.001
    prob <- (counts + pseudo) / (colSums(counts) + length(possible_BASES) * pseudo)
    score_mat <- prob
    for (b in possible_BASES) { score_mat[b, ] <- log(score_mat[b, ] / genome_freq[b]) }
    if (any(!is.finite(score_mat))) stop("Non-finite values in motif score matrix")
    N<-rep(min(score_mat),ncol(score_mat)) #to punish N bases in reference genome
    score_mat<-rbind(score_mat,N=N)  # to account for possible N bases in the scoring window    
    list(
      Score_Mat     = score_mat,
      Flank5_Length = flank5_len,
      Core_Length   = core_len,
      Flank3_Length = flank3_len,
      Motif_Length  = flank5_len + core_len + flank3_len,
      Score_Cutoff  = score_cutoff,
      Delta_Score   = delta_cutoff,
      Core_Key      = key)
  }

  ## ---------------- Donor motif ----------------
  Donor_counts <- rbind(
    A=c(69078,133203,20306,344,68,124042,144269,18185,37190,61333,46865,45907,47024,46191,45237),
    C=c(76100,22636,5909,90,2306,6204,16214,12413,31747,39876,52585,55648,50072,51271,50437),
    G=c(38921,24475,169126,208677,140,72594,24918,162563,40437,62928,50038,51226,54130,55092,54561),
    T=c(25093,28878,13851,81,206678,6352,23791,16031,99818,45055,59704,56411,57966,56638,58957)
  )
  colnames(Donor_counts) <- c(
    "N_3","N_2","N_1","G","T",
    "P_1","P_2","P_3","P_4","P_5",
    "P_6","P_7","P_8","P_9","P_10"
  )

  ## ---------------- Acceptor motif ----------------
  Acceptor_counts <- rbind(
    A=c(40734,36613,33102,29628,27245,24987,23178,21133,19501,17699,18037,20296,22224,23471,18254,19145,48523,12085,204984,83,52969,50620,52803),
    C=c(53377,54812,55903,56448,57231,58971,57997,57779,58066,54220,59175,61115,66801,68831,71330,61779,56893,133511,129,345,29914,39336,48198),
    G=c(30363,30345,29272,28350,27850,26412,25050,24009,21881,20869,22035,23326,21519,18864,13311,13277,42431,1020,125,204804,98960,40403,48740),
    T=c(80855,83559,87052,90903,93003,94959,99104,102408,105881,112541,106082,100592,94785,94163,102434,111128,57482,58713,91,97,23486,74970,55588)
  )
  colnames(Acceptor_counts) <- c(
    "N_18","N_17","N_16","N_15","N_14","N_13","N_12","N_11","N_10","N_9",
    "N_8","N_7","N_6","N_5","N_4","N_3","N_2","N_1",
    "A","G","P_1","P_2","P_3"
  )

  ## ---------------- Codon table ----------------
  CodonTable <- c(
    TTT="F",TTC="F",TTA="L",TTG="L",TCT="S",TCC="S",TCA="S",TCG="S",
    TAT="Y",TAC="Y",TAA="*",TAG="*",TGT="C",TGC="C",TGA="*",TGG="W",
    CTT="L",CTC="L",CTA="L",CTG="L",CCT="P",CCC="P",CCA="P",CCG="P",
    CAT="H",CAC="H",CAA="Q",CAG="Q",CGT="R",CGC="R",CGA="R",CGG="R",
    ATT="I",ATC="I",ATA="I",ATG="M",ACT="T",ACC="T",ACA="T",ACG="T",
    AAT="N",AAC="N",AAA="K",AAG="K",AGT="S",AGC="S",AGA="R",AGG="R",
    GTT="V",GTC="V",GTA="V",GTG="V",GCT="A",GCC="A",GCA="A",GCG="A",
    GAT="D",GAC="D",GAA="E",GAG="E",GGT="G",GGC="G",GGA="G",GGG="G"
  )

  ## ---------------- Reference files ----------------
  RefSeq <- list(
    hg19 = "/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/refs/ncbiRefSeqCurated_hg19_default.txt",
    hg38 = "/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/refs/ncbiRefSeqCurated_hg38_default.txt"
  )

  GenomeFa <- list(
    hg19 = "/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/refs/genome_fa/hg19.no_alt.fa",
    hg38 = "/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/refs/genome_fa/hg38.no_alt.fa"
  )

  stopifnot(file.exists(RefSeq[[build]]))
  stopifnot(file.exists(GenomeFa[[build]]))

  ## ---------------- Ensure FASTA index ----------------
  fai_file <- paste0(GenomeFa[[build]], ".fai")
  if (!file.exists(fai_file)) {
    cmd <- sprintf("%s faidx %s",shQuote(Tools$samtools),shQuote(GenomeFa[[build]]))
    rc <- system(cmd)
    if (rc != 0) stop("samtools faidx failed with exit code ", rc)
    stopifnot(file.exists(fai_file))
  }

  ## ---------------- Final constants object ----------------
  list(
    Donor_Motif    = build_motif(Donor_counts, 3, 2, 10, 3, 1, "GT"),
    Acceptor_Motif = build_motif(Acceptor_counts, 18, 2, 3, 3, 1, "AG"),
    CodonTable     = CodonTable,
    RefSeq         = RefSeq[[build]],
    GenomeFa       = GenomeFa[[build]],
    Tools          = Tools
  )
}

###############################################################################
## FASTA extraction via samtools faidx — FIXED
###############################################################################
extract_DNA_segs_using_samtools <- function(genomefa,samtools,chr,start,end,prefix) {
  ## ----------------------------
  ## Input validation
  if (length(chr) != length(start) || length(start) != length(end)) stop("chr, start, and end must have the same length")
  if (anyNA(chr) || anyNA(start) || anyNA(end)) stop("chr, start, and end must not contain NA")
  if (any(start < 1L) || any(end < start)) stop("Invalid genomic coordinates")
  stopifnot(file.exists(genomefa),file.exists(paste0(genomefa, ".fai")),file.exists(samtools),file.access(samtools,1) == 0)
  n <- length(chr)
  if (n == 0L) {
    return(data.frame(Chr = character(),Start = integer(),End = integer(),Seqs = character(),stringsAsFactors = FALSE))
  }
  ## ----------------------------
  ## Prepare region list
  regions <- paste0(chr, ":", start, "-", end)
  bed_file <- tempfile(paste0(prefix, ".bed"))
  fa_file  <- tempfile(paste0(prefix, ".fa"))
  writeLines(regions, bed_file)
  ## Ensure cleanup even on error
  on.exit(unlink(c(bed_file, fa_file)), add = TRUE)
  ## ----------------------------
  ## Call samtools faidx (CHECK STATUS)
  cmd <- paste(c(shQuote(samtools), "faidx", shQuote(genomefa), "-r", shQuote(bed_file), ">",  shQuote(fa_file)),collapse=" ") #keep this command final
  rc <- system(cmd) #keep this command final
  if (rc != 0 || !file.exists(fa_file) || file.info(fa_file)$size == 0) {
    stop("samtools faidx failed.")
  }
  ## ----------------------------
  ## Parse FASTA
  lines <- readLines(fa_file, warn = FALSE)
  hdr <- which(startsWith(lines, ">"))
  if (length(hdr) != n) {
    stop("FASTA entry count (", length(hdr),") does not match requested regions (", n, ")")
  }
  seqs <- vapply(seq_along(hdr), function(i) {
    s1 <- hdr[i] + 1L
    s2 <- if (i < length(hdr)) hdr[i + 1L] - 1L else length(lines)
    if (s1 <= s2) {
      paste(lines[s1:s2], collapse = "")
    } else {
      ""
    }
  }, character(1))
  if (any(nchar(seqs) == 0L)) stop("One or more regions returned empty sequences (check contig names and coordinates)")
  ## ----------------------------
  ## Return result
  data.frame(Chr = chr,Start = start,End = end,Seqs = toupper(seqs),stringsAsFactors = FALSE)
}

################################
### Extract one refseq sequence
################################
extract_sequences_for_one_transcript_in_genomic_order <- function(genomefa,samtools,gene_def,nmID) {
  idx <- which(gene_def$nmID == nmID)
  if (length(idx) != 1L) stop("nmID must match exactly one transcript: ", nmID)
  chr    <- gene_def$chr[idx]
  stopifnot(length(unique(chr)) == 1L)
  strand <- gene_def$strand[idx]
  cdsS <- as.integer(gene_def$cdsS[idx]) + 1L
  cdsE <- as.integer(gene_def$cdsE[idx])
  if (cdsS >= cdsE) stop("Transcript has no CDS: ", nmID)
  ## Parse exon coordinates
  exonS <- as.integer(strsplit(gene_def$exonS[idx], ",", fixed = TRUE)[[1]]) + 1L
  exonE <- as.integer(strsplit(gene_def$exonE[idx], ",", fixed = TRUE)[[1]])
  if (length(exonS) != length(exonE) || !length(exonS)) stop("Invalid exon definition for ", nmID)
  if (any(exonS > exonE)) stop("Exon start > exon end for ", nmID)
  ord <- order(exonS) #keep the increasing order so that we only extract the genomic position, regardless of gene_strand
  exonS <- exonS[ord]
  exonE <- exonE[ord]
  exonS[1] <- exonS[1] - 30L #extending a bit distance to ensure splice analysis; we have verified that all genes have txS and txE at least 30 bp away from chromosome ends
  exonE[length(exonE)]<-exonE[length(exonE)]+30 #extending a bit distance to ensure splice analysis; we have verified that all genes have txS and txE at least 30 bp away from chromosome ends

  ## Extract exon sequences
  exon_df <- extract_DNA_segs_using_samtools(genomefa,samtools,rep(chr, length(exonS)),exonS,exonE,paste0("tx_", nmID))
  if (nrow(exon_df) == 0L) stop("no exons extracted")
  ## Build transcript bases ordered by genomic positions, agnostic to the strand of the transcript
  if (any(exon_df$Start > exon_df$End, na.rm = TRUE)) stop("extraction region has wrong order")
  gpos_list  <- vector("list", nrow(exon_df))
  gbase_list <- vector("list", nrow(exon_df))
  for (i in seq_len(nrow(exon_df))) {
    tmp_gpos  <- exon_df$Start[i] : exon_df$End[i]
    tmp_gbase <- strsplit(exon_df$Seqs[i], "", fixed = TRUE)[[1]]
    if (length(tmp_gpos) != length(tmp_gbase)) stop("extraction region has problem at exon row ", i)
    gpos_list[[i]]  <- tmp_gpos
    gbase_list[[i]] <- tmp_gbase
  }
  gpos  <- unlist(gpos_list,  use.names = FALSE)
  gbase <- unlist(gbase_list, use.names = FALSE)
  ord   <- order(gpos)
  gpos  <- gpos[ord]
  gbase <- gbase[ord]
  ## Annotate regions in genomic space
  region <- rep("CDS", length(gpos)) #only coding genes in refseq will be analyzed in this code; this is fine
  region[gpos < cdsS] <- "UTR5"
  region[gpos > cdsE] <- "UTR3"
  tx_base<-gbase
  tx_pos<-seq_along(gbase)
  if(strand=="-") {
    region[gpos < cdsS] <- "UTR3"
    region[gpos > cdsE] <- "UTR5"    
    tx_base[gbase=="A"]<-"T"
    tx_base[gbase=="T"]<-"A"
    tx_base[gbase=="C"]<-"G"
    tx_base[gbase=="G"]<-"C"
    tx_base[!(gbase %in% c("A","C","G","T"))]<-"N"
    tx_pos<-rev(tx_pos)
  }
  ord<-order(tx_pos)
  rlt<-data.frame(Genomic_Pos=as.integer(gpos), Genomic_Base=gbase, Region=region,Transcript_Base=tx_base,Transcript_Pos=tx_pos)[ord,]
  rlt
}

#################################################
## extract splice sequence, plus strand only
#################################################
extract_splice_seqs_plus <- function(refseq,consts,splice_type,flank5,flank3) {
  CORE_LEN <- 2L
  WIN_LEN  <- flank5 + CORE_LEN + flank3
  all_regions <- character(0)
  for (i in seq_len(nrow(refseq))) {    
    if (refseq$strand[i] != "+") next
    chr <- refseq$chr[i]
    exonS <- as.integer(strsplit(refseq$exonS[i], ",")[[1]]) + 1L
    exonE <- as.integer(strsplit(refseq$exonE[i], ",")[[1]])
    if (length(exonS) < 2L) next  # no splice sites
    ord <- order(exonS)
    exonS <- exonS[ord]
    exonE <- exonE[ord]
    if (splice_type == "donor") {
      ## exon → intron junction
      for (j in seq_len(length(exonS) - 1L)) {
        core_start <- exonE[j] + 1L   # GT
        s <- core_start - flank5
        e <- core_start + CORE_LEN + flank3 - 1L
        all_regions <- c(all_regions, sprintf("%s:%d-%d", chr, s, e))
      }
    }
    if (splice_type == "acceptor") {
      ## intron → exon junction
      for (j in 2:length(exonS)) {
        core_start <- exonS[j] - CORE_LEN  # AG
        s <- core_start - flank5
        e <- core_start + CORE_LEN + flank3 - 1L
        all_regions <- c(all_regions, sprintf("%s:%d-%d", chr, s, e))
      }
    }
  }
  if (!length(all_regions)) return(character(0))
  ## -------- batch fetch --------
  bed <- tempfile()
  fa  <- tempfile()
  writeLines(all_regions, bed)
  cmd <- sprintf("%s faidx %s -r %s > %s",shQuote(consts$Tools$samtools),shQuote(consts$GenomeFa),shQuote(bed),shQuote(fa))
  system(cmd)

  fa_lines <- readLines(fa)
  hdr <- which(startsWith(fa_lines, ">"))

  seqs <- vapply(seq_along(hdr), function(i) {
    s1 <- hdr[i] + 1L
    s2 <- if (i < length(hdr)) hdr[i + 1L] - 1L else length(fa_lines)
    paste(fa_lines[s1:s2], collapse = "")
  }, character(1))
  unlink(c(bed, fa))
  toupper(seqs[nchar(seqs) == WIN_LEN])
}

#################################################
## extract splice sequence, minus strand only
#################################################
extract_splice_seqs_minus <- function(refseq, consts, splice_type, flank5, flank3) {
  CORE_LEN <- 2L
  WIN_LEN  <- flank5 + CORE_LEN + flank3
  all_regions <- character(0)
  comp <- c(A="T", T="A", C="G", G="C", N="N")
  for (i in seq_len(nrow(refseq))) {
    if (refseq$strand[i] != "-") next
    chr <- refseq$chr[i]
    exonS <- as.integer(strsplit(refseq$exonS[i], ",")[[1]]) + 1L
    exonE <- as.integer(strsplit(refseq$exonE[i], ",")[[1]])
    if (length(exonS) < 2L) next
    ord <- order(exonS)
    exonS <- exonS[ord]
    exonE <- exonE[ord]
    if (splice_type == "donor") {
      ## minus strand donor = genomic acceptor boundary
      for (j in 2:length(exonS)) {
        core_start <- exonS[j] - CORE_LEN   # AC on genome (→ GT in transcript)
        s <- core_start - flank3
        e <- core_start + CORE_LEN + flank5 - 1L
        all_regions <- c(all_regions, sprintf("%s:%d-%d", chr, s, e))
      }
    }
    if (splice_type == "acceptor") {
      ## minus strand acceptor = genomic donor boundary
      for (j in seq_len(length(exonS) - 1L)) {
        core_start <- exonE[j] + 1L          # CT on genome (→ AG in transcript)
        s <- core_start - flank3
        e <- core_start + CORE_LEN + flank5 - 1L
        all_regions <- c(all_regions, sprintf("%s:%d-%d", chr, s, e))
      }
    }
  }
  if (!length(all_regions)) return(character(0))

  ## -------- batch fetch --------
  bed <- tempfile()
  fa  <- tempfile()
  writeLines(all_regions, bed)

  cmd <- sprintf("%s faidx %s -r %s > %s",shQuote(consts$Tools$samtools),shQuote(consts$GenomeFa),shQuote(bed),shQuote(fa))
  system(cmd)

  fa_lines <- readLines(fa)
  hdr <- which(startsWith(fa_lines, ">"))

  seqs <- vapply(seq_along(hdr), function(i) {
    s1 <- hdr[i] + 1L
    s2 <- if (i < length(hdr)) hdr[i + 1L] - 1L else length(fa_lines)
    paste(fa_lines[s1:s2], collapse = "")
  }, character(1))

  unlink(c(bed, fa))

  ## reverse-complement to transcript orientation
  seqs <- toupper(seqs)
  seqs <- vapply(seqs, function(s) {
    b <- strsplit(s, "")[[1]]
    paste(rev(comp[b]), collapse = "")
  }, character(1))

  seqs[nchar(seqs) == WIN_LEN]
}

#################################################
## pileup motifs
#################################################

pileup_motif_seqs <- function(seqs) {
  bases <- c("A","C","G","T","N")
  L <- unique(nchar(seqs))
  if (length(L) != 1L) stop("Unequal sequence lengths")

  mat <- matrix(
    0L,
    nrow = 5,
    ncol = L,
    dimnames = list(bases, paste0("Pos", seq_len(L)))
  )

  for (s in seqs) {
    b <- strsplit(s, "")[[1]]
    for (i in seq_along(b)) {
      base <- ifelse(b[i] %in% bases, b[i], "N")
      mat[base, i] <- mat[base, i] + 1L
    }
  }
  mat
}

obtain_donor_acceptor_count_file<-function(build,Donor_Flank5,Donor_Flank3,Acceptor_Flank5,Acceptor_Flank3,ofn_donor,ofn_acceptor) {
  consts<-initialize_constants(build)
  refseq<-read.table(consts$RefSeq,header=TRUE,sep="\t")
  refseq<-refseq[refseq$isDefault=="Yes",]
  donor_seqs_P<-extract_splice_seqs_plus(refseq,consts,"donor",Donor_Flank5,Donor_Flank3)
  donor_pileup_P<-pileup_motif_seqs(donor_seqs_P)
  acceptor_seqs_P<-extract_splice_seqs_plus(refseq,consts,"acceptor",Acceptor_Flank5,Acceptor_Flank3)
  acceptor_pileup_P<-pileup_motif_seqs(acceptor_seqs_P)
  donor_seqs_M<-extract_splice_seqs_minus(refseq,consts,"donor",Donor_Flank5,Donor_Flank3)
  donor_pileup_M<-pileup_motif_seqs(donor_seqs_M)
  acceptor_seqs_M<-extract_splice_seqs_minus(refseq,consts,"acceptor",Acceptor_Flank5,Acceptor_Flank3)
  acceptor_pileup_M<-pileup_motif_seqs(acceptor_seqs_M)  
  donor_pileup<-donor_pileup_P + donor_pileup_M
  acceptor_pileup<-acceptor_pileup_P + acceptor_pileup_M
  colnames(donor_pileup)<-c(paste0("N_",Donor_Flank5:1),"G","T",paste0("P_",1:Donor_Flank3))
  colnames(acceptor_pileup)<-c(paste0("N_",Acceptor_Flank5:1),"G","T",paste0("P_",1:Acceptor_Flank3))
  donor_pileup<-cbind(Base=row.names(donor_pileup),donor_pileup)
  acceptor_pileup<-cbind(Base=row.names(acceptor_pileup),acceptor_pileup)
  write.table(donor_pileup,ofn_donor,row.names=F,col.names=TRUE,sep="\t",quote=F)
  write.table(acceptor_pileup,ofn_donor,row.names=F,col.names=TRUE,sep="\t",quote=F)
}

#########################################
## analyze coding changes for both strand, by storing genomic positions but analyzing transcript postions and transcript bases
#########################################
analyze_coding_aachanges<- function(cds_info, chr, gene_strand, codon_table) {
  required_cols <- c("Genomic_Pos", "Genomic_Base", "Region","Transcript_Base","Transcript_Pos")
  if (!all(required_cols %in% names(cds_info))) stop("cds_info must contain columns: ",paste(required_cols, collapse = ", "))
  gene_strand <- as.character(gene_strand)
  if (!gene_strand %in% c("+","-")) stop("gene_strand must be '+' or '-'") 

  the_cds<-cds_info; ord<-order(the_cds$Transcript_Pos); the_cds<-the_cds[ord,]
  keep <- the_cds$Region == "CDS"
  if (!any(keep)) {
    return(data.frame(Mutation = character(),Region = character(),Class = character(),AAChange = character(),stringsAsFactors = FALSE))
  }
  comp <- c(A="T", T="A", C="G", G="C", N="N")
  tbases <- toupper(the_cds$Transcript_Base[keep])
  gpos  <- as.integer(the_cds$Genomic_Pos[keep])
  gbases<-toupper(the_cds$Genomic_Base[keep])
  if (any(!tbases %in% c("A","C","G","T","N"))) stop("Non-ACGTN bases found in CDS")

  CDS <- paste0(tbases, collapse = ""); L <- nchar(CDS) #allow CDS to not start from ATG; such as genes MYC and PTEN
  if (L %% 3L != 0L) stop("CDS length is not a multiple of 3")

  codon_start <- seq.int(1L, L, by = 3L) 
  codon_end   <- codon_start + 2L
  codons <- substring(CDS, codon_start, codon_end)

  ref_aa <- codon_table[codons]
  ref_aa[is.na(ref_aa)] <- "X"
  theM <- which(ref_aa == "M")
  nextM <- theM[theM > 1L]
  dist_to_nextM <- if (length(nextM)) { min(nextM) - 1L } else { length(ref_aa) - 1L }

  ## Base index codon index & codon position
  codon_idx <- ((seq_len(L) - 1L) %/% 3L) + 1L
  codon_pos <- ((seq_len(L) - 1L) %% 3L) + 1L

  ## Preallocate result storage
  valid_BASES <- c("A","C","G","T"); k <- 0L
  n_rows <- sum(tbases %in% valid_BASES) * 3L
  Mutation <- character(n_rows); Region <- rep("CDS", n_rows); Class <- character(n_rows); AAChange <- character(n_rows)
  for (i in seq_len(L)) {
    t_ref_base  <- tbases[i]
    genomic_ref_base<-gbases[i]
    if(t_ref_base %in% valid_BASES) {
      ref_codon <- codons[codon_idx[i]]
      ref_aa_i  <- ref_aa[codon_idx[i]]
      for (t_alt in valid_BASES) {
        if (t_alt == t_ref_base) next
        k <- k + 1L
        mut_codon <- ref_codon
        substr(mut_codon, codon_pos[i], codon_pos[i]) <- t_alt
        alt_aa <- codon_table[mut_codon]
        if (is.na(alt_aa)) alt_aa <- "X"
        if (codon_idx[i] == 1L && ref_aa_i == "M" && alt_aa != "M") {
          Class[k] <- "startloss"
          AAChange[k]<-paste0("M1",alt_aa,"_",dist_to_nextM,"AA_to_nextM")
        } else if (ref_aa_i == alt_aa) {
          Class[k] <- "synonymous"
          AAChange[k] <- paste0(ref_aa_i, codon_idx[i], alt_aa)
        } else if (ref_aa_i != "*" && alt_aa == "*") {
          Class[k] <- "nonsense"
          AAChange[k] <- paste0(ref_aa_i, codon_idx[i], alt_aa)
        } else if (ref_aa_i == "*" && alt_aa != "*") {
          Class[k] <- "stoploss"
          AAChange[k] <- paste0(ref_aa_i, codon_idx[i], alt_aa)
        } else {
          Class[k] <- "missense"
          AAChange[k] <- paste0(ref_aa_i, codon_idx[i], alt_aa)
        }
        genomic_alt_base <- if (gene_strand == "-") comp[t_alt] else t_alt
        Mutation[k] <- paste(chr, gpos[i], genomic_ref_base, genomic_alt_base, sep = ".")
      }
    }
  }
  unique(data.frame(Mutation = Mutation[seq_len(k)],Region = Region[seq_len(k)],Class = Class[seq_len(k)],AAChange = AAChange[seq_len(k)],stringsAsFactors = FALSE))
}

########################
## analyze neo start codon 
########################
detect_neo_start_codons_UTR5 <- function(cds_info, chr, gene_strand, codon_table) {
  gene_strand <- as.character(gene_strand)
  if (!gene_strand %in% c("+","-")) stop("gene_strand must be '+' or '-'")

  comp <- c(A="T", T="A", C="G", G="C", N="N")
  the_cds <- cds_info
  ord <- order(the_cds$Transcript_Pos)
  the_cds <- the_cds[ord, ]

  UTR5_idx <- which(the_cds$Region == "UTR5")
  CDS_idx  <- which(the_cds$Region == "CDS")
  if (length(UTR5_idx) <= 2) {
    return(data.frame(Mutation = character(),Region = character(),Class = character(),AAChange = character(),stringsAsFactors = FALSE))
  }

  valid_BASES <- c("A","C","G","T")
  tx_seq <- toupper(the_cds$Transcript_Base)
  cds_start_tx_pos <- min(CDS_idx)

  ## list accumulation
  res_list <- vector("list", length(UTR5_idx) * 3L * 3L)  # safe upper bound
  k <- 0L

  ## helper: translate until stop or end
  translate_until_stop <- function(seq, start_pos) {
    aa <- character(0)
    for (j in seq(from = start_pos, to = length(seq) - 2, by = 3)) {
      codon <- paste0(seq[j:(j + 2)], collapse = "")
      aa_sym <- codon_table[[codon]]
      if (is.null(aa_sym)) aa_sym <- "X"
      if (aa_sym == "*") return(list(stop = TRUE, len = length(aa)))
      aa <- c(aa, aa_sym)
    }
    list(stop = FALSE, len = length(aa))
  }
  ## enumerate all single-base mutations in UTR5
  for (i in UTR5_idx) {
    t_ref_base <- tx_seq[i]
    if (!t_ref_base %in% valid_BASES) next
    for (t_alt_base in setdiff(valid_BASES, t_ref_base)) {
      mut_seq <- tx_seq
      mut_seq[i] <- t_alt_base
      ## ATG could start at i-2, i-1, or i
      for (A_pos_of_ATG in (i - 2):i) {
        if (A_pos_of_ATG < 1 || A_pos_of_ATG + 2 > length(mut_seq)) next
        if (!(i %in% A_pos_of_ATG:(A_pos_of_ATG + 2))) next
        if (paste0(mut_seq[A_pos_of_ATG:(A_pos_of_ATG + 2)], collapse = "") != "ATG") next
        if (paste0(tx_seq[A_pos_of_ATG:(A_pos_of_ATG + 2)], collapse = "") == "ATG") next
        ## ATG must be fully contained in UTR5
        if (!all((A_pos_of_ATG:(A_pos_of_ATG + 2)) %in% UTR5_idx)) next
        frame_offset <- (cds_start_tx_pos - A_pos_of_ATG) %% 3
        tr <- translate_until_stop(mut_seq, A_pos_of_ATG)
        if (tr$stop) {
          orf_last_coding_base <- A_pos_of_ATG + tr$len * 3 - 1
          if (orf_last_coding_base < cds_start_tx_pos + 2) {
            class <- "ruORF"
            aachg <- paste0(tr$len, "_AA")
          } else if (frame_offset == 0 && A_pos_of_ATG < cds_start_tx_pos && orf_last_coding_base >= cds_start_tx_pos + 2) {
            class <- "euORF"
            ext_aa <- ((cds_start_tx_pos - A_pos_of_ATG) %/% 3)
            aachg <- paste0(ext_aa, "_AA_extended_at_N_terminus")
          } else {
            class <- "iuORF"
            cds_aa_run <- floor((orf_last_coding_base - cds_start_tx_pos + 1) / 3)
            aachg <- paste0(cds_aa_run, "_AA_runaway")
          }
        } else {
          class <- "iuORF"
          aachg <- paste0(tr$len, "_AA_total_no_stop")
        }
        genomic_alt_base <- if (gene_strand == "-") comp[t_alt_base] else t_alt_base
        mut_str <- paste0(chr, ".", the_cds$Genomic_Pos[i], ".",the_cds$Genomic_Base[i], ".", genomic_alt_base)
        k <- k + 1L
        res_list[[k]] <- data.frame(Mutation = mut_str,Region = "UTR5",Class = class,AAChange = aachg,stringsAsFactors = FALSE)
      }
    }
  }
  if (k == 0L) {
    return(data.frame(Mutation = character(),Region = character(),Class = character(),AAChange = character(),stringsAsFactors = FALSE))
  }
  unique(do.call(rbind, res_list[seq_len(k)]))
}

###############################################################################
## score abolish (donor/acceptor)
###############################################################################
score_affinity_with_mutagenesis <- function(the_dna,start_rown,end_rown,gene_strand,motif,gene_chr, splice_type, region) {
  gene_strand <- as.character(gene_strand)
  if (!gene_strand %in% c("+","-")) stop("gene_strand must be '+' or '-'")
  valid_BASES <- c("A","C","G","T")
  comp <- c(A="T", T="A", C="G", G="C", N="N")
  .empty_df <- function() {
    data.frame(Mutation = character(),Region = character(),Class = character(),AAChange = character(),stringsAsFactors = FALSE)
  }
  ## ----------------------------
  ## Column normalization (accept both naming styles)
  if (!(sum(c("Transcript_Pos","Genomic_Pos","Genomic_Base","Transcript_Base") %in% names(the_dna))==4)) stop("the_dna must contain Transcript_Pos Genomic_Pos Genomic_Base Transcript_Base")
  tx_pos<- as.integer(the_dna$Transcript_Pos); ord <- order(tx_pos); the_dna <- the_dna[ord, , drop = FALSE];  #now the data is truly ordered by transcript strand
  tx_base <- toupper(as.character(the_dna$Transcript_Base)); tx_pos<- as.integer(the_dna$Transcript_Pos) 

  ## validate bounds
  start_rown <- as.integer(start_rown); end_rown <- as.integer(end_rown)
  if (start_rown < 1L || end_rown > nrow(the_dna) || start_rown > end_rown) stop("Invalid start_rown/end_rown range")

  my_idx <- start_rown:end_rown
  ref_splice_flank <- tx_base[my_idx]

  if (any(!ref_splice_flank %in% valid_BASES)) return(.empty_df())
  ## ----------------------------
  ## Reference affinity score
  score_mat <- motif$Score_Mat
  if (length(my_idx) != ncol(score_mat)) stop("Window length does not match motif length")
  
  ref_affinity_score <- 0
  for (i in seq_along(my_idx)) {
    b <- ref_splice_flank[i]
    if(b %in% valid_BASES) ref_affinity_score <- ref_affinity_score + score_mat[b, i]
  }
  ## ----------------------------
  ## Canonical core check (fixes core/need typo bugs)
  flank5   <- motif$Flank5_Length
  core_len <- motif$Core_Length
  core_start <- flank5 + 1L
  core_end   <- flank5 + core_len
  if (core_end > length(ref_splice_flank)) stop("Motif core extends beyond provided flank window")
  core_seq <- paste0(ref_splice_flank[core_start:core_end], collapse = "")  #columns of score mat are designed to match the expected bases for donor and accpetor, respectively
  core_key <- motif$Core_Key #columns of score mat are designed to match the expected bases for donor and accpetor, respectively
  isCanonical <- identical(core_seq, core_key)
  ## ----------------------------
  ## Enumerate single-base mutations and score impact
  res_list <- vector("list", length(my_idx) * 3L)  # safe upper bound
  k <- 0L
  for (j in seq_along(my_idx)) {
    t_ref_base <- ref_splice_flank[j]
    if (!t_ref_base %in% valid_BASES) next
    for (t_alt_base in valid_BASES) {
      if (t_alt_base == t_ref_base) next #make sure there is change
      ## Build mutated flank (this was missing in your original code)
      alt_flank <- ref_splice_flank
      alt_flank[j] <- t_alt_base
      alt_affinity_score <- 0
      for (i in seq_along(my_idx)) {
        b <- alt_flank[i]
        if(b %in% valid_BASES) alt_affinity_score <- alt_affinity_score + score_mat[b, i]
      }
      ## Abolish rule: must drop below ref by at least Delta_Score
      abolish <- (alt_affinity_score <= (ref_affinity_score - motif$Delta_Score)) # This is the final biological rule in this program
      ## If canonical splice key and mutation is inside core positions, force classify as splice_site
      in_core <- (j >= core_start && j <= core_end)
      alt_core_seq <- paste0(alt_flank[core_start:core_end], collapse = "")
      if (isCanonical && in_core && alt_core_seq != core_key) { 
        abolish <- TRUE
      }
      if (!abolish) next
      tmp_Class <- paste0(splice_type,"_abolish")
      gpos <- the_dna$Genomic_Pos[my_idx[j]]
      ref_gb <- toupper(the_dna$Genomic_Base[my_idx[j]])
      alt_gb <- if(gene_strand=="-") comp[t_alt_base] else t_alt_base #make it work for both strands
      k <- k + 1L
      mut_str <- paste(gene_chr, gpos, ref_gb, alt_gb, sep = ".")
      res_list[[k]] <- data.frame(Mutation = mut_str,Region = region,Class = tmp_Class,AAChange = "exon_lost" ,stringsAsFactors = FALSE)
    }
  }
  if (k == 0L) return(.empty_df())
  unique(do.call(rbind, res_list[seq_len(k)]))
}

###############################################################################
## detect donor/acceptor abolishing
###############################################################################

find_splice_abolish_hits_for_exon <- function(dna_info, regiontype, gene_chr, gene_strand, consts, buffer_len) {
  if (!(regiontype %in% c("first_exon","last_exon","internal_exon"))) stop("regiontype must be one of: first_exon, last_exon, internal_exon")
  .empty_df <- function() {
    data.frame(Mutation = character(),Region = character(),Class = character(),AAChange = character(),stringsAsFactors = FALSE)
  }
  max_motif_len <- max(consts$Donor_Motif$Motif_Length,consts$Acceptor_Motif$Motif_Length)
  if (buffer_len <= max_motif_len) stop("buffer_len must be greater than max motif length: ", max_motif_len)
  ## order by transcript position
  tx_pos <- as.integer(dna_info$Transcript_Pos)
  ord <- order(tx_pos)
  the_dna <- dna_info[ord, , drop = FALSE]
  ret <- .empty_df()
  ## ----------------------------
  ## Donor (downstream)
  if (regiontype %in% c("internal_exon","first_exon")) {
    motif<-consts$Donor_Motif
    my_start <- nrow(the_dna) - buffer_len - motif$Flank5_Length + 1L
    my_end   <- my_start + motif$Motif_Length - 1L
    if (my_start >= 1L && my_end <= nrow(the_dna)) {
      rlt <- score_affinity_with_mutagenesis(the_dna,my_start,my_end,gene_strand,motif,gene_chr,"Donor","splice")
      if (nrow(rlt)) { ret <- rbind(ret, rlt) }
    }
  }
  ## ----------------------------
  ## Acceptor (upstream)
  if (regiontype %in% c("internal_exon","last_exon")) {
    motif<-consts$Acceptor_Motif
    my_start <- buffer_len - motif$Flank5_Length - (motif$Core_Length-1) #buffer_len contains the core AG
    my_end   <- my_start + motif$Motif_Length - 1L
    if (my_start >= 1L && my_end <= nrow(the_dna)) {
      rlt <- score_affinity_with_mutagenesis(the_dna,my_start,my_end,gene_strand,motif,gene_chr,"Acceptor","splice")
      if (nrow(rlt)) { ret <- rbind(ret, rlt) }
    }
  }
  unique(ret)
}

########################################
## detect donor/acceptor establishing
########################################
find_neo_splice_sites_with_buffer <- function(dna_info, regiontype, buffer_len, gene_chr, gene_strand, consts) {
  ## ---------------------------- basic validation ----------------------------
  gene_strand <- as.character(gene_strand)
  if (!gene_strand %in% c("+","-")) stop("gene_strand must be '+' or '-'")
  if (!regiontype %in% c("exon","intron")) stop("regiontype must be one of: exon, intron")
  stopifnot(is.character(gene_chr), length(gene_chr) == 1L, nzchar(gene_chr))

  .empty_df <- function() {
    data.frame(Mutation = character(),Region = character(),Class = character(),AAChange = character(),stringsAsFactors = FALSE)
  }

  max_motif_len <- max(consts$Donor_Motif$Motif_Length, consts$Acceptor_Motif$Motif_Length)
  if (buffer_len <= max_motif_len) stop("buffer_len must be greater than max motif length")

  ## ---------------------------- motif selection ----------------------------
  if (regiontype == "exon") {
    motif <- consts$Donor_Motif; class_label <- "Neo_Donor"
  } else {
    motif <- consts$Acceptor_Motif; class_label <- "Neo_Acceptor" #in intron, looking for neo-acceptors
  }

  alt_BASES <- c("A","C","G","T")
  comp <- c(A="T", T="A", C="G", G="C", N="N")

  key_seq   <- strsplit(motif$Core_Key, "", fixed = TRUE)[[1]]
  core_len  <- motif$Core_Length
  flank5    <- motif$Flank5_Length
  motif_len <- motif$Motif_Length

  ## ---------------------------- order DNA ----------------------------
  ord <- order(as.integer(dna_info$Transcript_Pos))
  the_dna <- dna_info[ord, , drop = FALSE]

  tx_base <- toupper(the_dna$Transcript_Base)
  gb_base <- toupper(the_dna$Genomic_Base)
  gb_pos  <- as.integer(the_dna$Genomic_Pos)
  L <- length(tx_base)

  if (anyNA(tx_base) || anyNA(gb_base) || anyNA(gb_pos)) stop("dna_info contains NA values")

  ## ---------------------------- interior definition ----------------------------
  my_start <- buffer_len + 1L
  my_end   <- L - buffer_len
  if (my_start >= my_end) stop("Region too short for given buffer_len")

  interior_len <- my_end - my_start + 1L

  ## ---------------------------- score matrix ----------------------------
  score_mat <- as.matrix(motif$Score_Mat)
  storage.mode(score_mat) <- "numeric"

  score_window <- function(winbases) {
    if (length(winbases) != motif_len) return(NA_real_)
    if (any(!winbases %in% rownames(score_mat))) return(NA_real_)
    ret<-0
    for(i in 1:length(winbases)) ret<-ret + score_mat[winbases[i],i]
    ret
  }

  delta_cut <- motif$Delta_Score
  score_cut <- motif$Score_Cutoff

  ## ---------------------------- result storage ----------------------------
  res <- list()
  k <- 0L

  ## ---------------------------- main scan ----------------------------
  for (m in my_start:my_end) {

    ref_base <- tx_base[m]
    if (!ref_base %in% alt_BASES) next

    for (alt_base in setdiff(alt_BASES, ref_base)) {
      ## collect all hits for this SNV
      snv_hits <- list()
      h <- 0L
      ############################
      ## Branch 1: new core
      for (core_i in seq_len(core_len)) {

        core_start <- m - (core_i - 1L)
        core_end   <- core_start + core_len - 1L
        if (core_start < my_start || core_end > my_end) next

        ref_core <- tx_base[core_start:core_end]
        new_core <- ref_core
        new_core[core_i] <- alt_base

        if (paste0(new_core, collapse="") != motif$Core_Key) next
        if (paste0(ref_core, collapse="") == motif$Core_Key) next

        win_s <- core_start - flank5
        win_e <- win_s + motif_len - 1L
        if (win_s < 1L || win_e > L) next

        win_ref <- tx_base[win_s:win_e]
        win_alt <- win_ref
        win_alt[m - win_s + 1L] <- alt_base

        ref_score <- score_window(win_ref)
        alt_score <- score_window(win_alt)

        if (!is.finite(ref_score) || !is.finite(alt_score)) next
        if (alt_score < score_cut) next
        if ((alt_score - ref_score) < delta_cut) next

        dists <- if (class_label == "Neo_Donor") 
          interior_len - (core_start - my_start + 1L) + 1L
        else
          interior_len - (core_end   - my_start + 1L)

        h <- h + 1L
        snv_hits[[h]] <- list(core_start = core_start, Class = paste0(class_label, "_core"), 
                              AAChange = paste0(dists, ifelse(class_label=="Neo_Donor", "_nt_lost","_nt_gain")))
      }

      ############################
      ## Branch 2: flank strengthen
      ############################
      kstart_lo <- max(my_start, m - (motif_len - 1L) + flank5)
      kstart_hi <- min(my_end - core_len + 1L, m + flank5)

      for (kstart in kstart_lo:kstart_hi) {

        if (!identical(tx_base[kstart:(kstart + core_len - 1L)], key_seq)) next

        win_s <- kstart - flank5
        win_e <- win_s + motif_len - 1L
        if (win_s < 1L || win_e > L) next

        idx <- m - win_s + 1L
        if (idx < 1L || idx > motif_len) next

        win_ref <- tx_base[win_s:win_e]
        win_alt <- win_ref
        win_alt[idx] <- alt_base

        ref_score <- score_window(win_ref)
        alt_score <- score_window(win_alt)

        if (!is.finite(ref_score) || !is.finite(alt_score)) next
        if (alt_score < score_cut) next
        if ((alt_score - ref_score) < delta_cut) next

        dist <- if (class_label == "Neo_Donor")
          interior_len - (kstart - my_start + 1L) + 1L
        else
          interior_len - ((kstart + core_len - 1L) - my_start + 1L)

        h <- h + 1L
        snv_hits[[h]] <- list(core_start = kstart, Class      = paste0(class_label, "_flank", abs(kstart - m)),
                              AAChange   = paste0(dist, ifelse(class_label=="Neo_Donor", "_nt_lost","_nt_gain"))
        )
      }

      ## ---------------------------- keep most upstream hit ----------------------------
      if (length(snv_hits)) {
        best <- snv_hits[[ which.min(vapply(snv_hits, `[[`, 0L, "core_start")) ]]
        gb_alt <- if (gene_strand == "-") comp[alt_base] else alt_base
        k <- k + 1L
        res[[k]] <- data.frame(Mutation = paste(gene_chr, gb_pos[m], gb_base[m], gb_alt, sep="."),
          Region = regiontype,Class = best$Class,AAChange = best$AAChange, stringsAsFactors = FALSE
        )
      }
    }
  }

  if (k == 0L) return(.empty_df())
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  unique(out)
}

####################################
## consolidate
####################################
dedup_one_transcript <- function(anno_rlt) {
  ndt <- anno_rlt
  ## ---------------- abolish shadows opposite-side neo ----------------
  ## Donor_abolish shadows Neo_Acceptor; Acceptor_abolish shadows Neo_Donor
  to_keep <- rep(TRUE, nrow(ndt))
  muts <- unique(ndt$Mutation)
  for (mut in muts) {
    idx <- which(ndt$Mutation == mut)
    cls <- as.character(ndt$Class[idx])
    has_donor_abolish <- any(cls == "Donor_abolish", na.rm = TRUE)
    has_acc_abolish   <- any(cls == "Acceptor_abolish", na.rm = TRUE)
    has_neo_acc       <- any(grepl("^Neo_Acceptor", cls), na.rm = TRUE)
    has_neo_don       <- any(grepl("^Neo_Donor", cls), na.rm = TRUE)
    if (isTRUE(has_donor_abolish && has_neo_acc)) {
      to_keep[idx[grepl("^Neo_Acceptor", cls)]] <- FALSE
    }
    if (isTRUE(has_acc_abolish && has_neo_don)) {
      to_keep[idx[grepl("^Neo_Donor", cls)]] <- FALSE
    }
  }

  ndt <- ndt[to_keep, , drop = FALSE]
  rownames(ndt) <- NULL

  ## ---------------- neo shadows abolish (same side) ----------------
  to_keep <- rep(TRUE, nrow(ndt))
  dups <- table(ndt$Mutation)
  dups <- names(dups)[dups > 1]

  if (length(dups) > 0) {
    dups_dt <- ndt[ndt$Mutation %in% dups, ]
    unik_muts <- unique(dups_dt$Mutation)
    for (mut in unik_muts) {
      cls <- as.character(dups_dt$Class[dups_dt$Mutation == mut])
      has_acc_abolish <- any(cls == "Acceptor_abolish", na.rm = TRUE)
      has_neo_acc     <- any(grepl("^Neo_Acceptor", cls), na.rm = TRUE)
      if (isTRUE(has_acc_abolish && has_neo_acc)) {
        to_keep[ndt$Mutation == mut & ndt$Class == "Acceptor_abolish"] <- FALSE
      }
      has_don_abolish <- any(cls == "Donor_abolish", na.rm = TRUE)
      has_neo_don     <- any(grepl("^Neo_Donor", cls), na.rm = TRUE)
      if (isTRUE(has_don_abolish && has_neo_don)) {
        to_keep[ndt$Mutation == mut & ndt$Class == "Donor_abolish"] <- FALSE
      }
    }
  }
  ndt <- ndt[to_keep, , drop = FALSE]
  rownames(ndt) <- NULL

  ## ---------------- second pass: CDS + splice merge ----------------
  to_keep <- rep(TRUE, nrow(ndt))
  res <- list()
  k <- 0L

  coding_classes <- c("missense", "synonymous","stoploss", "startloss", "nonsense")
  splice_pattern <- "(Neo_Donor|Neo_Acceptor|Donor_abolish|Acceptor_abolish)"

  muts <- unique(ndt$Mutation)
  for (mut in muts) {
    idx <- which(ndt$Mutation == mut)
    if (length(idx) != 2L) next #only deal with mutations with 2 annotations
    cls <- as.character(ndt$Class[idx])
    is_coding <- cls %in% coding_classes
    is_splice <- grepl(splice_pattern, cls)
    if (sum(is_coding, na.rm = TRUE) == 1L && sum(is_splice, na.rm = TRUE) == 1L) {
      coding_row <- ndt[idx[is_coding], ]
      splice_row <- ndt[idx[is_splice], ]
      merged <- coding_row
      merged$Region   <- "CDS_splice"
      merged$Class    <- paste0(coding_row$Class, "_", splice_row$Class)
      merged$AAChange <- paste0(coding_row$AAChange, "_", splice_row$AAChange)
      k <- k + 1L
      res[[k]] <- merged
      to_keep[idx] <- FALSE
    }
  }
  if (k > 0L) {
    ndt <- rbind(ndt[to_keep, , drop = FALSE],do.call(rbind, res))
  } else {
    ndt <- ndt[to_keep, , drop = FALSE]
  }
  rownames(ndt) <- NULL

  
  ## ---------------- Neo ATG (ORF) merges with Donor_abolish ----------------
  ## If mutation annotated with an ORF class and Donor_abolish, MERGE them. ORF record takes precedence (keeps ORF row's main identity), but we keep splice info.
  to_keep <- rep(TRUE, nrow(ndt))
  res <- list()
  k <- 0L

  for (mut in unique(ndt$Mutation)) {
    idx <- which(ndt$Mutation == mut)
    if (length(idx) != 2L) next #only deal with mutations with 2 annotations
    cls <- as.character(ndt$Class[idx])
    is_orf <- grepl("^(ruORF|euORF|iuORF)", cls) # ruORF/euORF/iuORF etc.
    is_dab <- grepl("Neo_Donor|Donor_abolish",cls) # donor change
    if (any(is_orf) && any(is_dab)) {
      # pick one ORF row to be the base (priority)
      orf_row <- ndt[idx[which(is_orf)[1]], , drop = FALSE]
      dab_row <- ndt[idx[which(is_dab)[1]], , drop = FALSE]
      merged <- orf_row
      merged$Region   <- "UTR5_splice"
      merged$Class    <- paste0(orf_row$Class, "_", dab_row$Class)
      merged$AAChange <- paste0(orf_row$AAChange, "_", dab_row$AAChange)
      k <- k + 1L
      res[[k]] <- merged
      to_keep[idx[is_dab]] <- FALSE
      to_keep[idx[is_orf]] <- FALSE
    }
  }
  if (k > 0L) {
    ndt <- rbind(ndt[to_keep, , drop = FALSE], do.call(rbind, res))
    rownames(ndt) <- NULL
  } else {
    ndt <- ndt[to_keep, , drop = FALSE]
    rownames(ndt) <- NULL
  }
  
  ndt
}

#####################################
## analyze one transcript
#####################################

annoverse_one_transcript<-function(defaults, nmID, consts,prefix) {
  max_motif_len<-max(consts$Acceptor_Motif$Motif_Length,consts$Donor_Motif$Motif_Length)
  buffer_len<-max_motif_len+3
  my_idx<-which(defaults$nmID==nmID)
  comp <- c(A="T", T="A", C="G", G="C", N="N")
  final<-data.frame(Mutation = character(),Region = character(),Class = character(),AAChange = character(),stringsAsFactors = FALSE)
  if(length(my_idx)==1) {
    gene_chr<-defaults$chr[my_idx]
    gene_name<-defaults$gName[my_idx]
    gene_nmID<-defaults$nmID[my_idx]
    gene_strand<-defaults$strand[my_idx]
    gene_txS<-as.integer(defaults$txS[my_idx])
    gene_txE<-as.integer(defaults$txE[my_idx])
    gene_cdsS<-as.integer(defaults$cdsS[my_idx])
    gene_cdsE<-as.integer(defaults$cdsE[my_idx])
    gene_exonS<-as.integer(unlist(strsplit(as.character(defaults$exonS[my_idx]),split=",")))
    gene_exonE<-as.integer(unlist(strsplit(as.character(defaults$exonE[my_idx]),split=",")))
    ord<-order(gene_exonS); gene_exonS<-gene_exonS[ord]; gene_exonE<-gene_exonE[ord] #ensure genomic strand    
    gene_nEexon<-as.integer(defaults$exonN[my_idx])
    if(gene_nEexon==length(gene_exonS) && gene_nEexon==length(gene_exonE) && gene_nEexon>=1) {
      cds_info<-extract_sequences_for_one_transcript_in_genomic_order(consts$GenomeFa,consts$Tools$samtools,defaults,nmID) 
      ret <- analyze_coding_aachanges(cds_info, gene_chr, gene_strand, consts$CodonTable)
      utr5_alts   <- detect_neo_start_codons_UTR5(cds_info, gene_chr, gene_strand, consts$CodonTable)
      ret<-rbind(ret,utr5_alts)
      if(gene_nEexon>1) {
        for(i in seq_len(gene_nEexon)) {
          my_exonS<-gene_exonS[i]-buffer_len
          my_exonE<-gene_exonE[i]+buffer_len
          if(my_exonS<my_exonE) {
            regiontype<-"internal_exon"
            if(i==1 && gene_strand=="+") regiontype <- "first_exon"
            if(i==1 && gene_strand=="-") regiontype <- "last_exon"
            if(i==gene_nEexon && gene_strand=="+") regiontype <- "last_exon"
            if(i==gene_nEexon && gene_strand=="-") regiontype <- "first_exon"
            the_exon<-extract_DNA_segs_using_samtools(consts$GenomeFa,consts$Tools$samtools,gene_chr,my_exonS,my_exonE,prefix)
            Genomic_Pos<-the_exon$Start:the_exon$End #here we process exon one by one (only one row)
            Genomic_Base<-unlist(strsplit(the_exon$Seqs,split=""))
            Transcript_Pos<-seq_along(Genomic_Pos)
            Transcript_Base<-Genomic_Base
            Region<-rep("exon",length(Genomic_Pos))
            if(gene_strand=="-") {
              Transcript_Base<-comp[Transcript_Base]
              Transcript_Base[is.na(Transcript_Base)] <- "N"
              Transcript_Pos<-rev(Transcript_Pos)
            }
            dna_info<-data.frame(Genomic_Pos=Genomic_Pos,Genomic_Base=Genomic_Base,Transcript_Pos=Transcript_Pos,Transcript_Base=Transcript_Base)
            abolish_alts<-find_splice_abolish_hits_for_exon(dna_info, regiontype, gene_chr, gene_strand, consts, buffer_len)
            ret<-rbind(ret,abolish_alts)
            if(regiontype!="last_exon") {
              regiontype <- "exon"
              neo_splice_alts<-find_neo_splice_sites_with_buffer(dna_info,regiontype,buffer_len,gene_chr,gene_strand,consts)
              ret<-rbind(ret,neo_splice_alts)
            }
          }
          if(i<gene_nEexon) {
            my_intronS<-gene_exonE[i]-buffer_len+1
            my_intronE<-gene_exonS[i+1]+buffer_len-1
            if(my_intronS<my_intronE) {
              regiontype<-"intron"
              the_intron<-extract_DNA_segs_using_samtools(consts$GenomeFa,consts$Tools$samtools,gene_chr,my_intronS,my_intronE,prefix)
              Genomic_Pos<-the_intron$Start:the_intron$End
              Genomic_Base<-unlist(strsplit(the_intron$Seqs,split=""))
              Transcript_Pos<-seq_along(Genomic_Pos)
              Transcript_Base<-Genomic_Base
              Region<-rep("intron",length(Genomic_Pos))
              if(gene_strand=="-") {
                Transcript_Base<-comp[Transcript_Base]
                Transcript_Base[is.na(Transcript_Base)] <- "N"
                Transcript_Pos<-rev(Transcript_Pos)
              }
              dna_info<-data.frame(Genomic_Pos=Genomic_Pos,Genomic_Base=Genomic_Base,Transcript_Pos=Transcript_Pos,Transcript_Base=Transcript_Base)
              neo_splice_alts<-find_neo_splice_sites_with_buffer(dna_info,regiontype,buffer_len,gene_chr,gene_strand,consts)
              ret<-rbind(ret,neo_splice_alts)
            }
          }
        }
      }
    }
    final<-data.frame(Mutation = ret$Mutation, GeneName= gene_name, RefSeq = gene_nmID, Region = ret$Region,Class = ret$Class,AAChange = ret$AAChange,stringsAsFactors = FALSE)
    final<-dedup_one_transcript(final)
  }  
  return(final)
}


#################################
##annoverse for default nmIDs
#################################
annoverse_snv_defaults <-function(build) {
  consts<-initialize_constants(build)
  refseq<-read.table(consts$RefSeq,header=TRUE,sep="\t")
  defaults<-refseq[refseq$isDefault=="Yes",]
  nmIDs<-unique(defaults$nmID)
  nmIDs<-grep("NM_",nmIDs,value=TRUE)  
  for(i in seq_along(nmIDs)) {
    nmID<-nmIDs[i]
    print(paste0(nmID, "/", i))
    rlt<-annoverse_one_transcript(defaults, nmID, consts,prefix=nmID)
    ofn<-paste(nmID,"_anno.out.txt",sep="")
    write.table(rlt,ofn,row.names=F,col.names=TRUE,sep="\t",quote=F)
  }
}

#obtain_donor_acceptor_count_file(build,Donor_Flank5=3,Donor_Flank3=6,Acceptor_Flank5=20,Acceptor_Flank3=2,ofn_donor,ofn_acceptor)

annoverse_snv_defaults("hg38")
#annoverse_snv_defaults("hg19")


