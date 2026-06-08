##################################################################################################
## initialize_constants
##
## PURPOSE
##   Initialize and validate ALL global constants, reference paths,
##   scoring models, and PROJECT-LEVEL ANALYSIS POLICIES used throughout
##   the indel / splice / translation annotation pipeline.
##
##   This function is the SINGLE SOURCE OF TRUTH for:
##     - splice motif definitions (donor / acceptor)
##     - codon translation table
##     - reference file locations
##     - external tool paths
##     - buffer sizes and coordinate expansion rules
##     - project-level admissibility and execution policies
##
##
## ===================== PROJECT-LEVEL POLICY PARAMETERS =====================
##
##   build
##     Genome build identifier.
##
##     Allowed values:
##       - "hg19"
##       - "hg38"
##
##
##   max_indel_length
##     Maximum allowed length (in bp) of Ref or Mut bases for an indel
##     to be analyzed in THIS PROJECT.
##
##     - Applies to BOTH Ref and Mut sequences
##     - Insertions, deletions, and substitutions exceeding this limit
##       are excluded at TERMINAL orchestration layers
##     - This is a HARD project-level constraint
##
##
##   enforce_refbase_check
##     Logical flag controlling whether reference-base validation is
##     REQUIRED when converting genomic indels to transcript-local
##     coordinates.
##
##     - TRUE
##         * Enforce exact reference-base consistency
##         * Safer, slower
##
##     - FALSE
##         * Coordinate-only conversion
##         * Faster, permissive
##
##
## ===================== SPLICE MOTIF CONSTRUCTION POLICY =====================
##
##   Splice motifs are constructed via an INTERNAL helper function:
##
##       build_motif()
##
##   Key design rules:
##
##     - Canonical splice cores are specified EXCLUSIVELY as
##       SINGLE character strings:
##
##           Donor    → "GT"
##           Acceptor → "AG"
##
##     - The core string length MUST equal core_len
##     - Case is normalized internally
##     - The core string is split internally into per-base keys
##
##   NOTE:
##     build_motif() is NOT a public API. It is defined inside
##     initialize_constants() and must NOT be called elsewhere.
##
##
## ===================== BUFFER / WINDOW DEFINITIONS =====================
##
##   Tx_End_Buffer
##     Extra transcript buffer used near transcript ends.
##
##   E_I_end_buffer
##     Genomic buffer added to exon/intron extraction windows.
##     Ensures sufficient flanking sequence for splice motif scanning.
##
##   Splice_Region_Length
##     Length (bp) defining the canonical splice region context.
##
##
## ===================== FAILURE POLICY =====================
##
##   - This function MAY throw via stop() / stopifnot()
##   - Failure here is considered FATAL
##   - No downstream function is expected to recover from
##     an initialize_constants() failure
##
##################################################################################################
initialize_constants <- function(build,
                                 max_indel_length = 50L, enhancer_length = 50000L,
                                 enforce_refbase_check = FALSE) {

  ## ---------------- Basic validation ----------------
  stopifnot(build %in% c("hg19", "hg38"))
  max_indel_length <- as.integer(max_indel_length)
  stopifnot(!is.na(max_indel_length), max_indel_length > 0L)
  stopifnot(is.logical(enforce_refbase_check), length(enforce_refbase_check) == 1L)

  ## ---------------- Tool paths ----------------
  Tools <- list(
    samtools = "/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/apps/samtools"
  )
  stopifnot(file.exists(Tools$samtools))
  stopifnot(file.access(Tools$samtools, 1) == 0)  # executable

  ## ==========================================================================
  ## Donor motif (GT)
  ## ==========================================================================
  Donor_counts <- rbind(
    A = c(69078,133203,20306,344,68,124042,144269,18185,37190,61333,46865,45907,47024,46191,45237),
    C = c(76100,22636,5909,90,2306,6204,16214,12413,31747,39876,52585,55648,50072,51271,50437),
    G = c(38921,24475,169126,208677,140,72594,24918,162563,40437,62928,50038,51226,54130,55092,54561),
    T = c(25093,28878,13851,81,206678,6352,23791,16031,99818,45055,59704,56411,57966,56638,58957)
  )
  colnames(Donor_counts) <- c(
    "N_3","N_2","N_1","G","T",
    "P_1","P_2","P_3","P_4","P_5",
    "P_6","P_7","P_8","P_9","P_10"
  )

  ## ==========================================================================
  ## Acceptor motif (AG)
  ## ==========================================================================
  Acceptor_counts <- rbind(
    A = c(40734,36613,33102,29628,27245,24987,23178,21133,19501,17699,18037,20296,22224,23471,18254,19145,48523,12085,204984,83,52969,50620,52803),
    C = c(53377,54812,55903,56448,57231,58971,57997,57779,58066,54220,59175,61115,66801,68831,71330,61779,56893,133511,129,345,29914,39336,48198),
    G = c(30363,30345,29272,28350,27850,26412,25050,24009,21881,20869,22035,23326,21519,18864,13311,13277,42431,1020,125,204804,98960,40403,48740),
    T = c(80855,83559,87052,90903,93003,94959,99104,102408,105881,112541,106082,100592,94785,94163,102434,111128,57482,58713,91,97,23486,74970,55588)
  )
  colnames(Acceptor_counts) <- c(
    "N_18","N_17","N_16","N_15","N_14","N_13","N_12","N_11","N_10",
    "N_9","N_8","N_7","N_6","N_5","N_4","N_3","N_2","N_1",
    "A","G","P_1","P_2","P_3"
  )


  ## ---------------- Global buffers / limits ----------------
  Tx_End_Buffer         <- 30L
  exon_intron_buffer   <- max_indel_length + max(ncol(Acceptor_counts),ncol(Donor_counts)) + 50L
  Splice_Region_Length <- 10L
  Enhancer_Length <- enhancer_length

  ## ==========================================================================
  ## INTERNAL splice motif constructor
  ## ==========================================================================
  build_motif <- function(counts,
                           f5_len,
                           core_len,
                           f3_len,
                           score_cutoff,
                           delta_abolish,
                           delta_enhance,
                           key_string,
                           core_shift) {

    possible_BASES <- c("A","C","G","T")

    ## ---------------- Core key normalization ----------------
    if (!is.character(key_string) || length(key_string) != 1L)
      stop("key_string must be a character scalar")

    key_string <- toupper(key_string)

    if (nchar(key_string) != core_len)
      stop("key_string length must equal core_len")

    Core_Key_Bases <- strsplit(key_string, "", fixed = TRUE)[[1]]

    ## ---------------- Genome background ----------------
    genome_freq <- c(A = 0.3, C = 0.2, G = 0.2, T = 0.3)
    genome_freq <- genome_freq / sum(genome_freq)

    ## ---------------- Counts validation ----------------
    if (!all(possible_BASES %in% rownames(counts)))
      stop("counts must have rownames A/C/G/T")

    counts <- counts[possible_BASES, , drop = FALSE]

    stopifnot(ncol(counts) == f5_len + core_len + f3_len)
    stopifnot(core_shift >= 0L && core_shift < core_len)

    ## ---------------- PWM construction ----------------
    pseudo <- 0.5
    prob <- (counts + pseudo) /
      (colSums(counts) + length(possible_BASES) * pseudo)

    score_mat <- prob
    for (b in possible_BASES)
      score_mat[b, ] <- log(score_mat[b, ] / genome_freq[b])

    if (any(!is.finite(score_mat)))
      stop("non-finite values in motif score matrix")

    ## ambiguity base contributes neutral score
    score_mat <- rbind(score_mat, N = rep(0, ncol(score_mat)))
    score_mat <- score_mat[c("A","C","G","T","N"), , drop = FALSE]

    list(
      Score_Mat      = score_mat,
      F5_Len         = f5_len,
      Core_Length    = core_len,
      F3_Len         = f3_len,
      Motif_Length   = f5_len + core_len + f3_len,
      Score_Cutoff   = score_cutoff,
      Delta_Abolish  = delta_abolish,
      Delta_Enhance  = delta_enhance,
      Core_Shift     = core_shift,
      Core_Key       = key_string,
      Core_Key_Bases = Core_Key_Bases	  
    )
  }


  ## ==========================================================================
  ## Codon table
  ## ==========================================================================
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

  ## ==========================================================================
  ## Reference files
  ## ==========================================================================
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
    cmd <- sprintf("%s faidx %s",
                   shQuote(Tools$samtools),
                   shQuote(GenomeFa[[build]]))
    rc <- system(cmd)
    if (rc != 0) stop("samtools faidx failed with exit code ", rc)
  }

  ## ==========================================================================
  ## Final constants object
  ## ==========================================================================
  list(
    ## motifs
    Donor_Motif = build_motif(
      Donor_counts,
      f5_len        = 3,
      core_len      = 2,
      f3_len        = 10,
      score_cutoff  = 3,
      delta_abolish = 1,
      delta_enhance = 2,
      key_string    = "GT",
      core_shift    = 0
    ),

    Acceptor_Motif = build_motif(
      Acceptor_counts,
      f5_len        = 18,
      core_len      = 2,
      f3_len        = 3,
      score_cutoff  = 3,
      delta_abolish = 1,
      delta_enhance = 2,
      key_string    = "AG",
      core_shift    = 1
    ),

    ## translation
    CodonTable = CodonTable,

    ## references / tools
    RefSeq   = RefSeq[[build]],
    GenomeFa = GenomeFa[[build]],
    Tools    = Tools,

    ## buffers
    Tx_End_Buffer        = Tx_End_Buffer,
    Splice_Region_Length = Splice_Region_Length,
    E_I_end_buffer       = exon_intron_buffer,

    ## project-level policy
    Max_Indel_Length       = max_indel_length,
    Enforce_Refbase_Check  = enforce_refbase_check,
	Enhancer_Length        = Enhancer_Length
  )
}



##################### reverse complement a DNA segment, reture a string 
revcomp_string <- function(myseq) {
  stopifnot(is.character(myseq), length(myseq) == 1L)
  myseq <- toupper(myseq)
  comp <- c(A="T", T="A", C="G", G="C", N="N")
  s <- strsplit(myseq, "", fixed = TRUE)[[1]]
  s <- rev(comp[ifelse(s %in% names(comp), s, "N")])
  return(paste0(s, collapse = ""))
}

##################### reverse complement a DNA segment, reture an array of letters
revcomp_array <- function(myseq) {
  comp <- c(A="T", T="A", C="G", G="C", N="N")
  if (length(myseq) > 1L) {
    s <- toupper(as.character(myseq))
  } else {
    s <- strsplit(toupper(myseq), "", fixed = TRUE)[[1]]
  }
  rev(comp[ifelse(s %in% names(comp), s, "N")])
}

##############################################################
## extract_DNA_segments_using_samtools
##
## PURPOSE
##   Extract genomic DNA sequences for one or more genomic regions
##   using `samtools faidx`.
##
## DESIGN CONTRACT (FAILURE-SAFE)
##   - NEVER throws
##   - ALWAYS returns a list
##   - Success/failure indicated by scalar `Status`
##   - On failure, sequences are empty strings
##
## INPUTS
##   my_chrs   : character vector of chromosome names (e.g. "chr1")
##   my_starts : integer-like numeric vector of 1-based inclusive start positions
##   my_ends   : integer-like numeric vector of 1-based inclusive end positions
##
##   my_chrs, my_starts, and my_ends must have identical length.
##   Scalar inputs are allowed and treated as length-1 requests.
##
##   prefix_for_tempfile : character scalar used as prefix for temp files
##   consts    : list containing:
##                 consts$GenomeFa        (path to indexed genome fasta)
##                 consts$Tools$samtools  (path to samtools executable)
##
## ===================== COORDINATE VALIDATION POLICY =====================
##
##   - my_starts and my_ends MUST be integer-valued genomic coordinates
##   - Numeric inputs are accepted IF AND ONLY IF they are integer-like
##     (e.g. 100 is allowed, 100.5 is rejected)
##   - Internally, all coordinates are normalized to integer storage mode
##
##   RATIONALE:
##     samtools faidx requires integer genomic coordinates.
##     Silent truncation (e.g. as.integer(100.9) → 100) is NOT allowed,
##     as it can produce biologically incorrect sequence extraction.
##
## ===================== OUTPUT =====================
##
##   list(
##     Status  = TRUE / FALSE,     # scalar logical
##     Message = character(1),     # empty on success
##     Data    = data.frame(
##                 Chr   = character,
##                 Start = integer,
##                 End   = integer,
##                 Seqs  = character     # empty strings on failure
##               )
##   )
##
## FAILURE POLICY
##   Any error at any stage results in:
##     - Status = FALSE
##     - Message = non-empty diagnostic string
##     - Data$Seqs = ""
##
##############################################################
extract_DNA_segments_using_samtools <- function(my_chrs,my_starts,my_ends,prefix_for_tempfile,consts) {
  n <- length(my_starts)
  ## ---------------- Default failure return (shape invariant) ----------------
  fail_ret <- list(
    Status  = FALSE,
    Message = "in extract_DNA_segments_using_samtools",
    Data = data.frame(
      Chr   = my_chrs,
      Start = my_starts,
      End   = my_ends,
      Seqs  = rep("", n),
      stringsAsFactors = FALSE
    )
  )
  ## ---------------- Basic validation ----------------
  if (length(my_chrs) != length(my_starts)) return(fail_ret)
  if (length(my_ends)  != length(my_starts)) return(fail_ret)
  if (n == 0L) return(fail_ret)
  if (!is.character(my_chrs)) return(fail_ret)
  if (any(is.na(my_chrs) | is.na(my_starts) | is.na(my_ends))) return(fail_ret)

  ## ---------------- Integer coordinate validation ----------------
  .is_integerish_vec <- function(x) { is.numeric(x) && !anyNA(x) && all(x == as.integer(x)) }
  if (!.is_integerish_vec(my_starts) || !.is_integerish_vec(my_ends)) return(fail_ret)

  ## normalize to integer storage mode
  my_starts <- as.integer(my_starts)
  my_ends   <- as.integer(my_ends)

  if (any(my_starts > my_ends)) return(fail_ret)

  ## ---------------- Tool availability ----------------
  if (is.null(consts$GenomeFa) || is.null(consts$Tools$samtools)) return(fail_ret)

  genomefa <- consts$GenomeFa
  samtools <- consts$Tools$samtools

  ## ---------------- Prepare temporary files ----------------
  uniq_prefix <- paste0(
    prefix_for_tempfile, "_",
    Sys.getpid(), "_",
    as.integer(runif(1, 1, 1e9))
  )
  region_file <- paste0(uniq_prefix, ".regions")
  fa_file     <- paste0(uniq_prefix, ".fa")
  on.exit(unlink(c(region_file, fa_file)), add = TRUE)

  ## ---------------- Write region specification ----------------
  regions <- paste0(my_chrs, ":", my_starts, "-", my_ends)
  writeLines(regions, region_file)

  ## ---------------- Call samtools faidx ----------------
  cmd <- paste(
    shQuote(samtools),
    "faidx",
    shQuote(genomefa),
    "-r",
    shQuote(region_file),
    ">",
    shQuote(fa_file)
  )

  rc <- suppressWarnings(system(cmd))
  if (rc != 0L || !file.exists(fa_file) || file.info(fa_file)$size == 0L)
    return(fail_ret)

  ## ---------------- Parse FASTA output ----------------
  lines <- readLines(fa_file, warn = FALSE)
  header_idx <- which(startsWith(lines, ">"))

  ## Each requested region must produce exactly one FASTA record
  if (length(header_idx) != n) return(fail_ret)

  seqs <- vapply(
    seq_along(header_idx),
    function(i) {
      s1 <- header_idx[i] + 1L
      s2 <- if (i < length(header_idx))
        header_idx[i + 1L] - 1L
      else
        length(lines)
      if (s1 <= s2) paste(lines[s1:s2], collapse = "") else ""
    },
    character(1)
  )

  if (any(nchar(seqs) == 0L)) return(fail_ret)

  ## ---------------- Success return ----------------
  list(
    Status  = TRUE,
    Message = "",
    Data = data.frame(
      Chr   = my_chrs,
      Start = my_starts,
      End   = my_ends,
      Seqs  = toupper(seqs),
      stringsAsFactors = FALSE
    )
  )
}


##################################################################################################
## extract_mRNA_in_transcript_orientation
##
## PURPOSE
##   Construct the transcript-oriented sequence for a single RefSeq transcript
##   (nmID), with per-base annotation suitable for:
##
##     - exon region extraction
##     - splice motif scanning near transcript ends
##     - codon-frame–aware translation analysis
##
##   Sequence is assembled by concatenating extracted genomic segments in
##   ascending genomic order, then (if strand == "-") reversing and
##   reverse-complementing to obtain transcript 5' → 3' orientation.
##
##
## ===================== IMPORTANT SEMANTICS =====================
##
##   EXTRACTION MODEL (as implemented here)
##     - Extract EXON segments only (NOT introns between exons)
##     - Additionally apply symmetric transcript-end buffering:
##         first segment start  = exonS[1] - consts$Tx_End_Buffer
##         last  segment end    = exonE[last] + consts$Tx_End_Buffer
##
##     Buffer bases outside annotated exons remain:
##       Region   = "Buffer"
##       G_ExonNo = -9
##       T_CodonFrame = -9
##
##
## ===================== INPUT =====================
##
##   nmID : character scalar
##
##   refseq : data.frame with at minimum:
##       nmID, chr, strand,
##       cdsS, cdsE,
##       exonS, exonE
##
##     Coordinates follow UCSC / RefSeq convention:
##       - cdsS, exonS are 0-based, half-open (start is excluded)
##       - cdsE, exonE are 0-based, end-exclusive
##
##   consts : list
##     Must contain:
##       - consts$GenomeFa
##       - consts$Tools$samtools
##       - consts$Tx_End_Buffer (non-negative integer)
##
##   prefix_for_tempfile : character scalar
##
##
## ===================== OUTPUT =====================
##
##   list(Status, Message, Data)
##
##   Data columns (row-aligned, transcript-oriented):
##
##     G_Pos        : integer (1-based, inclusive)
##     G_Base       : character (A/C/G/T/N)
##     Region       : "UTR5","CDS","UTR3","Buffer"
##     T_Pos_Local  : integer (1-based, contiguous)
##     T_Base       : transcript-oriented base
##     T_CodonFrame : 0/1/2 for CDS bases; -9 otherwise
##     G_ExonNo     : exon number (1-based); -9 otherwise
##
##
## ===================== FAILURE POLICY =====================
##
##   - NEVER throws
##   - Any failure returns Status = FALSE with empty Data
##
##################################################################################################
extract_mRNA_in_transcript_orientation <- function(nmID,refseq,consts,prefix_for_tempfile) {
  default_fail <- list(Status = FALSE, Message = "", Data = data.frame())
  fail_ret <- function() default_fail

  ## ---------------- Basic validation ----------------
  if (!is.character(nmID) || length(nmID) != 1L) return(fail_ret())
  if (!is.data.frame(refseq)) return(fail_ret())
  if (is.null(consts$GenomeFa) || is.null(consts$Tools$samtools)) return(fail_ret())
  if (is.null(consts$Tx_End_Buffer)) return(fail_ret())

  tx_end_buffer <- consts$Tx_End_Buffer
  if (!is.numeric(tx_end_buffer) || length(tx_end_buffer) != 1L || is.na(tx_end_buffer)) return(fail_ret())
  if (tx_end_buffer < 0) return(fail_ret())
  tx_end_buffer <- as.integer(tx_end_buffer)

  idx <- which(refseq$nmID == nmID)
  if (length(idx) != 1L) return(fail_ret())
  row <- refseq[idx, , drop = FALSE]

  ## ---------------- Transcript metadata ----------------
  chr    <- as.character(row$chr)
  strand <- as.character(row$strand)

  cdsS <- as.integer(row$cdsS) + 1L
  cdsE <- as.integer(row$cdsE)

  exonS <- as.integer(strsplit(row$exonS, ",", fixed = TRUE)[[1]]) + 1L
  exonE <- as.integer(strsplit(row$exonE, ",", fixed = TRUE)[[1]])

  if (length(exonS) == 0L || length(exonE) == 0L) return(fail_ret())
  if (length(exonS) != length(exonE)) return(fail_ret())
  if (any(is.na(exonS) | is.na(exonE))) return(fail_ret())

  ## Ensure exon order in genomic space
  ord_ex <- order(exonS, exonE, na.last = NA)
  exonS <- exonS[ord_ex]
  exonE <- exonE[ord_ex]
  if (any(exonS > exonE)) return(fail_ret())

  ## ---------------- Define extraction segments: exons + end buffers ----------------
  my_starts <- exonS
  my_ends   <- exonE

  my_starts[1] <- my_starts[1] - tx_end_buffer
  my_ends[length(my_ends)] <- my_ends[length(my_ends)] + tx_end_buffer

  ## clamp to valid genomic coordinates (1-based inclusive)
  if (my_starts[1] < 1L) my_starts[1] <- 1L
  if (any(my_starts > my_ends)) return(fail_ret())

  ## ---------------- Extract genomic DNA (source of truth) ----------------
  ext <- extract_DNA_segments_using_samtools(
    my_chrs   = rep(chr, length(my_starts)),
    my_starts = my_starts,
    my_ends   = my_ends,
    prefix_for_tempfile = prefix_for_tempfile,
    consts = consts
  )
  if (!isTRUE(ext$Status)) return(fail_ret())
  if (!is.data.frame(ext$Data) || nrow(ext$Data) == 0L) return(fail_ret())

  ## ---------------- Assemble per-base G_Pos / G_Base row-by-row ----------------
  ext_df <- ext$Data
  ext_df$Start <- as.integer(ext_df$Start)
  ext_df$End   <- as.integer(ext_df$End)

  ord <- order(ext_df$Start, ext_df$End, na.last = NA)
  ext_df <- ext_df[ord, , drop = FALSE]

  G_Pos  <- integer(0)
  G_Base <- character(0)

  for (i in seq_len(nrow(ext_df))) {
    s <- ext_df$Start[i]
    e <- ext_df$End[i]
    if (is.na(s) || is.na(e) || s > e) return(fail_ret())

    seq_i <- toupper(ext_df$Seqs[i])
    bases_i <- strsplit(seq_i, "", fixed = TRUE)[[1]]
    bases_i[!(bases_i %in% c("A","C","G","T"))] <- "N"

    pos_i <- seq.int(s, e)
    if (length(bases_i) != length(pos_i)) return(fail_ret())

    G_Pos  <- c(G_Pos,  pos_i)
    G_Base <- c(G_Base, bases_i)
  }

  n <- length(G_Pos)
  if (n == 0L || length(G_Base) != n) return(fail_ret())

  ## ---------------- Initialize fields ----------------
  Region       <- rep("Buffer", n)
  G_ExonNo     <- rep(-9L, n)
  T_CodonFrame <- rep(-9L, n)

  ## ---------------- Exon labeling (original exon boundaries only) ----------------
  for (j in seq_along(exonS)) {
    idx_ex <- G_Pos >= exonS[j] & G_Pos <= exonE[j]
    G_ExonNo[idx_ex] <- j
  }

  ## ---------------- Region classification (exonic only) ----------------
  is_exonic <- G_ExonNo > 0L
  Region[is_exonic & G_Pos <  cdsS]                  <- "UTR5"
  Region[is_exonic & G_Pos >= cdsS & G_Pos <= cdsE]  <- "CDS"
  Region[is_exonic & G_Pos >  cdsE]                  <- "UTR3"

  ## ---------------- Codon frame assignment ----------------
  cds_idx <- which(Region == "CDS")
  if (length(cds_idx) > 0L) {
    T_CodonFrame[cds_idx] <- (seq_along(cds_idx) - 1L) %% 3L
  }

  ## ---------------- Strand correction ----------------
  if (strand == "-") {
    comp <- c(A="T", T="A", C="G", G="C", N="N")

    G_Pos        <- rev(G_Pos)
    G_Base       <- rev(G_Base)
    G_ExonNo     <- rev(G_ExonNo)
    Region       <- rev(Region)
    T_CodonFrame <- rev(T_CodonFrame)
    Region[is_exonic & G_Pos <  cdsS]                  <- "UTR3"
    Region[is_exonic & G_Pos >= cdsS & G_Pos <= cdsE]  <- "CDS"
    Region[is_exonic & G_Pos >  cdsE]                  <- "UTR5"
    T_Base <- unname(comp[G_Base])
  } else {
    T_Base <- G_Base
  }

  ## ---------------- Transcript-local coordinates ----------------
  T_Pos_Local <- seq_len(n)

  ## ---------------- Assemble output ----------------
  Data <- data.frame(
    G_Pos        = G_Pos,
    G_Base       = G_Base,
    Region       = Region,
    T_Pos_Local  = T_Pos_Local,
    T_Base       = T_Base,
    T_CodonFrame = T_CodonFrame,
    G_ExonNo     = G_ExonNo,
    stringsAsFactors = FALSE
  )

  list(Status = TRUE, Message = "", Data = Data)
}




####################################################################################################
## extract_one_exon_or_intron_region_with_buffer
##
## PURPOSE
##   Extract a SINGLE exon or intron together with surrounding buffer
##   sequence, and return a transcript-oriented (5' → 3') base-level
##   representation suitable for:
##
##     - splice motif scanning
##     - frame-aware inspection
##     - low-level sequence comparison
##
##   This function is intentionally LOW-LEVEL and STRICT:
##     - It performs NO splice interpretation
##     - It performs NO biological classification
##     - It ONLY labels regions according to explicit genomic roles
##
##
## ===================== CRITICAL INVARIANT =====================
##
##   The input mRNA data frame is assumed to be VALIDATED UPSTREAM:
##
##     - Exon numbering is correct and contiguous
##     - Intron boundaries are correct
##     - For any intron, the immediate downstream exon (Exon_3)
##       is fully present and correctly annotated in mRNA
##
##   Under this invariant:
##     - Exon_3 is always well-defined
##     - The extraction window always fully contains:
##         Intron + Exon_3 + buffer
##
##
## ===================== REGION LABELING POLICY =====================
##
##   region_type == "Intron":
##
##     [ Buffer ] [ Intron ] [ Exon_3 ] [ Buffer ]
##
##   region_type == "Exon":
##
##     [ Buffer ] [ Exon ] [ Buffer ]
##
##
## ===================== COORDINATE & ORIENTATION =====================
##
##   - Genomic coordinates are 1-based, inclusive
##   - Output is ALWAYS transcript-oriented (5' → 3')
##
##   Strand handling:
##
##     gene_strand == "+"
##       - Genomic order preserved
##
##     gene_strand == "-"
##       - Genomic order REVERSED
##       - Bases are REVERSE-COMPLEMENTED
##       - Transcript-local coordinates are recomputed as 1..N
##
##   IMPORTANT:
##     All per-base fields remain ROW-ALIGNED after strand correction.
##
##
## ===================== OUTPUT =====================
##
##   list(
##     Status  = logical scalar,
##     Message = character scalar,
##     Data    = data.frame
##   )
##
##   Data columns:
##     G_Pos        : genomic position
##     G_Base       : genomic base (A/C/G/T/N)
##     Region       : {"Buffer","Exon","Intron","Exon_3"}
##     T_Pos_Local  : transcript-local coordinate (1-based)
##     T_Base       : transcript-oriented base
##     T_CodonFrame : codon frame copied from mRNA, or -9
##     G_ExonNo     : exon number (-9L for Buffer)
##
##
## ===================== FAILURE POLICY =====================
##
##   - The function NEVER throws
##   - Return shape is invariant
##   - Any error results in Status = FALSE
##
####################################################################################################
extract_one_exon_or_intron_region_with_buffer <- function(region_chr,region_start,region_end,region_number,region_type,tmpfile_prefix,gene_strand,mRNA,consts) {
  default_msg <- "in extract_one_exon_or_intron_region_with_buffer"
  fail_ret <- function(msg = default_msg) list(Status = FALSE, Message = msg, Data = data.frame())

  tryCatch({
    ## ---------------- Validation ----------------
    if (!region_type %in% c("Exon", "Intron")) return(fail_ret())
    if (!gene_strand %in% c("+", "-")) return(fail_ret())
    if (!is.integer(region_start) || !is.integer(region_end)) return(fail_ret())
    if (region_start > region_end) return(fail_ret())
    if (!is.data.frame(mRNA) || nrow(mRNA) == 0L) return(fail_ret())
    if (is.null(consts$E_I_end_buffer)) return(fail_ret())

    buffer_len <- as.integer(consts$E_I_end_buffer)

    ## ---------------- Downstream exon (intron only) ----------------
    if (region_type == "Intron") {
      exon3_no  <- region_number + 1L
      exon3_idx <- which(mRNA$G_ExonNo == exon3_no & mRNA$Region != "Buffer")
      if (!length(exon3_idx)) return(fail_ret())
      exon3_start <- min(mRNA$G_Pos[exon3_idx])
      exon3_end   <- max(mRNA$G_Pos[exon3_idx])
    }

    ## ---------------- Extraction window ----------------
    left_buffer   <- min(buffer_len, region_start - 1L)
    extract_start <- region_start - left_buffer
    extract_end   <- if (region_type == "Exon") region_end + buffer_len else exon3_end + buffer_len

    ## ---------------- Genomic DNA extraction ----------------
    ext <- extract_DNA_segments_using_samtools(my_chrs  = region_chr,my_starts = extract_start,my_ends   = extract_end,prefix_for_tempfile = tmpfile_prefix,consts = consts)
    if (!isTRUE(ext$Status)) return(fail_ret(ext$Message))

    seq <- toupper(ext$Data$Seqs[1])
    G_Base <- strsplit(seq, "", fixed = TRUE)[[1]]
    G_Base[!(G_Base %in% c("A","C","G","T"))] <- "N"

    G_Pos <- extract_start:extract_end
    n     <- length(G_Pos)

    ## ---------------- Initialize fields ----------------
    Region       <- rep("Buffer", n)
    G_ExonNo     <- rep(-9L, n)
    T_CodonFrame <- rep(-9L, n)

    ## ---------------- Core region labeling ----------------
    core_start_idx <- left_buffer + 1L

    if (region_type == "Exon") {
      core_end_idx <- core_start_idx + (region_end - region_start)
      Region[core_start_idx:core_end_idx]   <- "Exon"
      G_ExonNo[core_start_idx:core_end_idx] <- region_number
    } else {
      intron_len <- region_end - region_start + 1L
      intron_end_idx <- core_start_idx + intron_len - 1L

      Region[core_start_idx:intron_end_idx]   <- "Intron"
      G_ExonNo[core_start_idx:intron_end_idx] <- region_number

      exon3_len <- exon3_end - exon3_start + 1L
      exon3_start_idx <- intron_end_idx + 1L
      exon3_end_idx   <- exon3_start_idx + exon3_len - 1L

      Region[exon3_start_idx:exon3_end_idx]   <- "Exon_3"
      G_ExonNo[exon3_start_idx:exon3_end_idx] <- region_number + 1L
    }

    ## ---------------- Codon frame mapping (genomic) ----------------
    idx <- match(G_Pos, mRNA$G_Pos)
    hit <- !is.na(idx)
    T_CodonFrame[hit] <- mRNA$T_CodonFrame[idx[hit]]

    ## ---------------- Strand correction (CRITICAL FIX) ----------------
    if (gene_strand == "-") {
      comp <- c(A="T", T="A", C="G", G="C", N="N")

      ## reverse genomic order to transcript order
      G_Pos        <- rev(G_Pos)
      G_Base       <- rev(G_Base)
      Region       <- rev(Region)
      G_ExonNo     <- rev(G_ExonNo)
      T_CodonFrame <- rev(T_CodonFrame)

      ## reverse-complement bases
      T_Base <- comp[G_Base]
    } else {
      T_Base <- G_Base
    }

    ## transcript-local coordinate is ALWAYS recomputed
    T_Pos_Local <- seq_len(n)

    ## ---------------- Assemble output ----------------
    Data <- data.frame(
      G_Pos        = G_Pos,
      G_Base       = G_Base,
      Region       = Region,
      T_Pos_Local  = T_Pos_Local,
      T_Base       = T_Base,
      T_CodonFrame = T_CodonFrame,
      G_ExonNo     = G_ExonNo,
      stringsAsFactors = FALSE
    )

    list(Status = TRUE, Message = "", Data = Data)

  }, error = function(e) {
    fail_ret()
  })
}

##################################################################################################
## scan_splice_motif_in_mut_and_wt_DNA_segment_with_sufficient_buffer_in_both_ends
##
## PURPOSE
##   Scan splice motif binding in reference and mutant DNA contigs that are
##   ALREADY oriented in transcript space (5' → 3'), under a STRICT execution
##   policy.
##
##   Motif scanning is PERMISSIVE and performed on BOTH contigs using the
##   motif-specific score cutoff (motif$Score_Cutoff).
##
##   This function reports ONLY raw motif hit positions, scores, and the
##   exact DNA strings producing those hits.
##
##   CORE BASES are rendered in UPPER CASE and are surrounded by a SINGLE
##   BLANK SPACE on both sides for human readability.
##   Flanking bases remain contiguous and upper case.
##
##   Example:
##     Donor   : NNNGTNNNN   →   NNN GT NNNN
##     Acceptor: NNAGNNN     →   NN AG NNN
##
##   NO biological interpretation is performed here.
##
##
## ===================== CRITICAL POSITION SEMANTICS =====================
##
##   Motif hits are reported in TRANSCRIPT-LOCAL COORDINATES and correspond
##   to the ANCHOR BASE of the motif core:
##
##       hit_position = (start_of_motif_window)
##                      + motif$F5_Len
##                      + motif$Core_Shift
##
##   Canonical anchor definitions:
##
##     Donor motif ("GT"):
##       - Reported position : the 'G' of "GT"
##
##     Acceptor motif ("AG"):
##       - Reported position : the 'G' of "AG"
##
##
## ===================== STATUS SEMANTICS (STRICT) =====================
##
##   Status = TRUE
##     → motif scan WAS EXECUTED successfully.
##       Zero hits are allowed.
##
##   Status = FALSE
##     → motif scan was NOT executed due to validation or precondition
##       failure.
##
##
## ===================== INPUT =====================
##
##   ref_vec : character vector
##     Reference DNA bases (A/C/G/T/N), transcript-oriented (5' → 3')
##
##   mut_vec : character vector
##     Mutant DNA bases (A/C/G/T/N), transcript-oriented (5' → 3')
##
##   motif : list
##     Splice motif object created by build_motif(), containing:
##       - Motif_Length
##       - F5_Len
##       - Core_Length
##       - Core_Shift
##       - Core_Key_Bases
##       - Score_Mat (rows A/C/G/T/N)
##       - Score_Cutoff
##
##
## ===================== OUTPUT =====================
##
##   list(
##     Status  = logical scalar,
##     Message = character scalar,
##     Data    = list
##   )
##
##   Data list fields (valid only if Status == TRUE):
##
##     Ref_Hit_Pos     : integer vector
##     Ref_Hit_Score   : numeric vector (rounded to 1 decimal)
##     Ref_Hit_String  : character vector (motif window with spaced core)
##
##     Mut_Hit_Pos     : integer vector
##     Mut_Hit_Score   : numeric vector (rounded to 1 decimal)
##     Mut_Hit_String  : character vector (motif window with spaced core)
##
##
## ===================== FAILURE POLICY =====================
##
##   - NEVER throws
##   - Return shape is invariant
##   - Message is "" on success, fixed default string on failure
##
##################################################################################################
scan_splice_motif_in_mut_and_wt_DNA_segment_with_sufficient_buffer_in_both_ends <- function(
  ref_vec,
  mut_vec,
  motif
) {
  default_msg <- "in scan_splice_motif_in_mut_and_wt_DNA_segment_with_sufficient_buffer_in_both_ends"

  fail_ret <- function()
    list(Status = FALSE, Message = default_msg, Data = list())

  ## ===================== Motif validation =====================
  required_fields <- c(
    "Motif_Length","F5_Len","Core_Length","Core_Shift",
    "Core_Key_Bases","Score_Mat","Score_Cutoff"
  )
  if (!all(required_fields %in% names(motif))) return(fail_ret())
  if (!is.matrix(motif$Score_Mat) || !is.numeric(motif$Score_Mat)) return(fail_ret())
  if (is.null(rownames(motif$Score_Mat))) return(fail_ret())
  if (!all(c("A","C","G","T","N") %in% rownames(motif$Score_Mat))) return(fail_ret())
  if (any(!is.finite(motif$Score_Mat))) return(fail_ret())

  motif_len <- as.integer(motif$Motif_Length)
  F5        <- as.integer(motif$F5_Len)
  CS        <- as.integer(motif$Core_Shift)
  CL        <- as.integer(motif$Core_Length)
  if (any(is.na(c(motif_len, F5, CS, CL))) || motif_len < 1L) return(fail_ret())

  CKB <- motif$Core_Key_Bases
  if (length(CKB) == 1L && nchar(CKB) > 1L)
    CKB <- strsplit(CKB, "", fixed = TRUE)[[1]]
  if (length(CKB) != CL) return(fail_ret())

  valid_bases  <- rownames(motif$Score_Mat)
  score_cutoff <- motif$Score_Cutoff

  ## ===================== Input validation =====================
  if (!is.character(ref_vec) || !is.character(mut_vec)) return(fail_ret())
  if (!(all(nchar(ref_vec) == 1L) && all(nchar(mut_vec) == 1L))) return(fail_ret())

  len_ref <- length(ref_vec)
  len_mut <- length(mut_vec)
  min_len <- min(len_ref, len_mut)
  if (min_len < motif_len) return(fail_ret())

  ## ===================== Identify identical flanks =====================
  left_identical <- 0L
  for (i in seq_len(min_len)) {
    if (ref_vec[i] != mut_vec[i]) break
    left_identical <- i
  }
  if (left_identical == min_len) return(fail_ret())

  right_identical <- 0L
  max_right <- min_len - left_identical
  if (max_right > 0L) {
    for (i in seq_len(max_right)) {
      if (ref_vec[len_ref - i + 1L] != mut_vec[len_mut - i + 1L]) break
      right_identical <- i
    }
  }
  if (left_identical < 1L || right_identical < 1L) return(fail_ret())

  ## ===================== Search windows =====================
  search_start <- max(left_identical - motif_len + 1L, 1L)
  search_end_ref <- min(len_ref - motif_len + 1L, len_ref - right_identical + 1L)
  search_end_mut <- min(len_mut - motif_len + 1L, len_mut - right_identical + 1L)
  if (search_end_ref < search_start && search_end_mut < search_start)
    return(fail_ret())

  ## ===================== Deduplication helper =====================
  dedup_hits <- function(pos, score, str) {
    if (!length(pos))
      return(list(pos = integer(0), score = numeric(0), string = character(0)))
    o <- order(pos, -score)
    keep <- !duplicated(pos[o])
    list(
      pos    = pos[o][keep],
      score  = score[o][keep],
      string = str[o][keep]
    )
  }

  ## ===================== Motif scanning =====================
  scan_one <- function(vec, search_end) {
    hit_pos    <- integer(0)
    hit_score  <- numeric(0)
    hit_string <- character(0)

    for (j in seq.int(search_start, search_end)) {
      key_start <- j + F5
      key_end   <- key_start + CL - 1L
      if (key_end > length(vec)) next
      if (!all(vec[key_start:key_end] == CKB)) next

      score <- 0
      for (ii in seq_len(motif_len)) {
        b <- vec[j + ii - 1L]
        if (!(b %in% valid_bases)) b <- "N"
        score <- score + motif$Score_Mat[b, ii]
      }

      if (score >= score_cutoff) {
        win <- vec[j:(j + motif_len - 1L)]

        core_start <- F5 + 1L
        core_end   <- F5 + CL

        hit_string <- c(
          hit_string,
          paste0(
            paste0(win[1:(core_start - 1L)], collapse = ""),
            " ",
            paste0(win[core_start:core_end], collapse = ""),
            " ",
            paste0(win[(core_end + 1L):length(win)], collapse = "")
          )
        )

        hit_pos   <- c(hit_pos, key_start + CS)
        hit_score <- c(hit_score, round(score, 1))
      }
    }
    dedup_hits(hit_pos, hit_score, hit_string)
  }

  ref_hits <- if (search_end_ref >= search_start)
    scan_one(ref_vec, search_end_ref)
  else
    list(pos = integer(0), score = numeric(0), string = character(0))

  mut_hits <- if (search_end_mut >= search_start)
    scan_one(mut_vec, search_end_mut)
  else
    list(pos = integer(0), score = numeric(0), string = character(0))

  ## ===================== Success =====================
  list(
    Status  = TRUE,
    Message = "",
    Data    = list(
      Ref_Hit_Pos    = ref_hits$pos,
      Ref_Hit_Score  = ref_hits$score,
      Ref_Hit_String = ref_hits$string,
      Mut_Hit_Pos    = mut_hits$pos,
      Mut_Hit_Score  = mut_hits$score,
      Mut_Hit_String = mut_hits$string
    )
  )
}


##################################################################################################
## make_mut_vec_full_in_local_coordinates
##
## PURPOSE
##   Construct reference and mutant DNA base vectors in transcript-local
##   coordinates around a substitution, indel, or complex insertion event.
##
##   This is a LOW-LEVEL SEQUENCE PREPARATION utility. It performs:
##     - NO splice interpretation
##     - NO translation interpretation
##     - NO biological classification
##
##   It simply slices and assembles transcript-oriented base vectors
##   required by downstream logic.
##
##
## ===================== COORDINATE & EVENT SEMANTICS =====================
##
##   Coordinates are transcript-local and 1-based.
##
##   Let:
##     T_indel_start, T_indel_end define the mutation interval.
##
##   Cases:
##
##     1. T_indel_start == T_indel_end
##        - SNV
##        - ref and mut bases typically length 1
##
##     2. T_indel_start < T_indel_end
##        - MNV, deletion, or complex replacement
##        - ref and mut bases MAY be empty or non-empty
##
##     3. T_indel_start > T_indel_end
##        - insertion or complex insertion
##        - ref bases MUST be empty (character(0))
##        - mut bases MAY be empty or non-empty
##
##
## ===================== REFERENCE CHECK POLICY =====================
##
##   If to_perform_reference_check == TRUE:
##
##     - Only applies when T_indel_start <= T_indel_end
##     - Requires:
##         * T_ref_bases length matches (end - start + 1)
##         * T_ref_bases content EXACTLY matches transcript reference
##
##   Reference checking is SKIPPED for insertion cases
##   (T_indel_start > T_indel_end).
##
##
## ===================== INPUT =====================
##
##   this_exon_data_frame : data.frame
##     Must contain:
##       - T_Pos_Local   (1-based, contiguous)
##       - T_Base
##       - T_CodonFrame
##
##   T_indel_start, T_indel_end : integer
##     Transcript-local coordinates (1-based)
##
##   T_ref_bases : character vector
##     Reference bases at mutation core (may be empty)
##
##   T_mut_bases : character vector
##     Mutant bases at mutation core (may be empty)
##
##   to_perform_reference_check : logical
##     Whether to enforce reference consistency
##
##
## ===================== OUTPUT =====================
##
##   list(
##     Status  = logical scalar,
##     Message = character scalar,
##     Data    = list
##   )
##
##   Data fields:
##
##     ref_up          : transcript bases upstream of mutation
##     ref_core        : reference bases at mutation core
##     alt_core        : mutant bases at mutation core
##     ref_down        : transcript bases downstream of mutation
##
##     ref_up_coding   : CDS-only upstream bases
##     ref_down_coding : CDS-only downstream bases
##
##
## ===================== FAILURE POLICY =====================
##
##   - The function NEVER throws
##   - Return shape is invariant
##   - Message is:
##       * "" on success
##       * a fixed default string on failure
##
##################################################################################################
make_mut_vec_full_in_local_coordinates <- function(this_exon_data_frame,T_indel_start,T_indel_end,T_ref_bases,T_mut_bases,to_perform_reference_check = FALSE) {
  default_msg <- "in make_mut_vec_full_in_local_coordinates"
  fail_ret <- function() list(Status = FALSE, Message = default_msg, Data = list())

  ## ---------------- Defensive validation ----------------
  if (!is.data.frame(this_exon_data_frame) || nrow(this_exon_data_frame) < 1L) return(fail_ret())
  required_cols <- c("T_Pos_Local", "T_Base", "T_CodonFrame")
  if (!all(required_cols %in% names(this_exon_data_frame))) return(fail_ret())
  if (!all(diff(this_exon_data_frame$T_Pos_Local) == 1L)) return(fail_ret())
  if (!is.numeric(T_indel_start) || !is.numeric(T_indel_end)) return(fail_ret())
  T_indel_start <- as.integer(T_indel_start)
  T_indel_end   <- as.integer(T_indel_end)
  n <- nrow(this_exon_data_frame)
  ## allow start > end (insertion)
  if (T_indel_start < 1L || T_indel_end > n) return(fail_ret())

  ## base vectors must be character (empty allowed)
  if (!is.character(T_ref_bases) || !is.character(T_mut_bases)) return(fail_ret())

  ## ---------------- Event-type-specific constraints ----------------

  if (T_indel_start > T_indel_end) {
    ## insertion / complex insertion, we do not take note of the reference bases here. if yes, it is error
    if (length(T_ref_bases) != 0L)
      return(fail_ret())
    ## reference check does NOT apply here

  } else {
    ## substitution / deletion / complex replacement
    if (to_perform_reference_check) {
      expected_len <- T_indel_end - T_indel_start + 1L
      if (length(T_ref_bases) != expected_len)
        return(fail_ret())

      ref_from_tx <- this_exon_data_frame$T_Base[T_indel_start:T_indel_end]
      if (!identical(as.character(T_ref_bases), as.character(ref_from_tx)))
        return(fail_ret())
    }
  }

  ## ---------------- Construct transcript-level vectors ----------------
  ref_up <- if (T_indel_start > 1L) this_exon_data_frame$T_Base[seq_len(T_indel_start - 1L)] else character(0)
  ref_down <- if (T_indel_end < n) this_exon_data_frame$T_Base[seq.int(T_indel_end + 1L, n)] else character(0)

  ## ---------------- CDS-only context ----------------
  ref_up_coding <- this_exon_data_frame$T_Base[this_exon_data_frame$T_Pos_Local < T_indel_start & this_exon_data_frame$T_CodonFrame >= 0]
  ref_down_coding <- this_exon_data_frame$T_Base[this_exon_data_frame$T_Pos_Local > T_indel_end & this_exon_data_frame$T_CodonFrame >= 0]
  ## ---------------- Assemble output ----------------
  Data <- list(
    ref_up          = ref_up,
    ref_core        = T_ref_bases,
    alt_core        = T_mut_bases,
    ref_down        = ref_down,
    ref_up_coding   = ref_up_coding,
    ref_down_coding = ref_down_coding
  )
  list(Status = TRUE, Message = "", Data = Data)
}




##################################################################################################
## translate_dna_frame0_vec
##
## PURPOSE
##   Translate a DNA base vector into an amino acid sequence assuming
##   the first base is in reading frame 0.
##
##   Translation proceeds codon-by-codon from the 5' end and stops at
##   the last complete codon. Any trailing bases that do not form a
##   full codon are ignored.
##
##   This is a LOW-LEVEL TRANSLATION utility. It performs:
##     - NO start-codon detection
##     - NO stop-codon truncation
##     - NO splice interpretation
##
##
## ===================== ASSUMPTIONS =====================
##
##   - dna_vec represents a contiguous coding sequence already
##     oriented in the correct transcriptional direction.
##   - The first base of dna_vec corresponds to frame 0.
##   - consts$CodonTable defines the genetic code to use.
##
##
## ===================== CODON HANDLING POLICY =====================
##
##   - Codons not present in the codon table (e.g. containing
##     ambiguous bases such as 'N') are translated as "X".
##   - Translation does NOT terminate at stop codons; stop symbols
##     are returned exactly as defined in the codon table.
##
##
## ===================== STATUS SEMANTICS (STRICT) =====================
##
##   Status = TRUE
##     → translation executed successfully AND
##       at least one complete codon was translated.
##
##   Status = FALSE
##     → translation NOT executed due to:
##         - missing or malformed codon table
##         - invalid dna_vec input
##         - fewer than 3 bases provided
##         - zero complete codons available
##
##
## ===================== INPUT =====================
##
##   dna_vec : character vector
##     DNA bases (A/C/G/T/N), oriented 5' → 3'
##
##   consts : list
##     Must contain:
##       - consts$CodonTable : named character vector
##         Names are codons (e.g. "ATG"), values are amino acids
##
##
## ===================== OUTPUT =====================
##
##   list(
##     Status  = logical scalar,
##     Message = character scalar,
##     Data    = character vector
##   )
##
##   Data:
##     Character vector of translated amino acids.
##     Length equals floor(length(dna_vec) / 3).
##
##
## ===================== FAILURE POLICY =====================
##
##   - The function NEVER throws
##   - Return shape is invariant
##   - Message is:
##       * "" on success
##       * a fixed default string on failure
##
##################################################################################################
translate_dna_frame0_vec <- function(dna_vec, consts) {
  default_msg <- "in translate_dna_frame0_vec"
  ## Unified return constructor (used for ALL returns)
  fail_ret <- function(Status = FALSE, Data = character(0)) { list(Status  = Status,Message = if (Status) "" else default_msg,Data    = Data) }
  ## ---------------- Validation ----------------
  if (is.null(consts$CodonTable)) return(fail_ret())
  if (!is.character(consts$CodonTable)) return(fail_ret())
  if (is.null(names(consts$CodonTable))) return(fail_ret())
  if (!is.character(dna_vec)) return(fail_ret())
  dna_vec <- toupper(dna_vec)
  ## STRICT: fewer than one codon is a failure
  if (length(dna_vec) < 3L) return(fail_ret())

  n_codons <- length(dna_vec) %/% 3L
  if (n_codons == 0L) return(fail_ret())
  ## ---------------- Codon construction ----------------
  dna_vec <- dna_vec[seq_len(n_codons * 3L)]
  codons <- paste0(
    dna_vec[seq(1L, by = 3L, length.out = n_codons)],
    dna_vec[seq(2L, by = 3L, length.out = n_codons)],
    dna_vec[seq(3L, by = 3L, length.out = n_codons)]
  )
  ## ---------------- Translation ----------------
  aa <- consts$CodonTable[codons]
  aa[is.na(aa)] <- "X"
  fail_ret(Status = TRUE, Data = unname(as.character(aa)))
}



## obtain_indel_notation_with_padding(c(1,2,3,3,5,6),c(1,2,3,4,5,6))
## obtain_indel_notation_with_padding(c(1,2,3,5,6),c(1,2,3,4,5,6))
## obtain_indel_notation_with_padding(c(1,2,3,4,5,6),c(1,2,3,4,4,5,6))
## obtain_indel_notation_with_padding(c(1,2,3,4,4,4,5,6),c(1,2,3,4,4,5,6))
##################################################################################################
## obtain_indel_notation_with_padding
##
## PURPOSE
##   Normalize an indel or substitution (SNV/MNV) event given reference
##   and mutant base vectors.
##
##
## ===================== POLICY =====================
##
##   Event type:
##     - INDEL        : length(ref_vec) != length(mut_vec)
##     - SUBSTITUTION : length(ref_vec) == length(mut_vec)
##
##   Padding rules:
##
##     - INDEL:
##         * Require ≥1 identical base on BOTH 5' and 3' ends
##         * EXACTLY 1 padding base on each end is INCLUDED:
##             - left pad  = last common base on 5' side
##             - right pad = first common base on 3' side
##
##     - SUBSTITUTION:
##         * NO padding requirement
##         * Minimal differing core is extracted
##
##   Failure rules:
##
##     - ref_vec == mut_vec → Status = FALSE
##     - insufficient padding for indel → Status = FALSE
##     - extracted substitution core identical → Status = FALSE
##
##
## ===================== OUTPUT =====================
##
##   list(
##     Status  = logical scalar,
##     Message = character scalar,
##     Data    = list
##   )
##
##   Data fields (valid only if Status == TRUE):
##
##     Class       : "insertion" | "deletion" | "substitution"
##     Change      : "REF.MUT"
##     Ref_Bases   : character vector
##     Mut_Bases   : character vector
##     Core_Start  : integer (1-based index in input vectors)
##
##################################################################################################
obtain_indel_notation_with_padding <- function(ref_vec, mut_vec) {
  default_msg <- "in obtain_indel_notation_with_padding"
  ret <- function(Status = FALSE, Data = list()) { list(Status  = Status,Message = if (Status) "" else default_msg,Data    = Data) }

  ## ---------------- Validation ----------------
  if (!is.atomic(ref_vec) || !is.atomic(mut_vec)) return(ret())
  if (length(ref_vec) == 0L || length(mut_vec) == 0L) return(ret())
  ref_vec <- as.character(ref_vec)
  mut_vec <- as.character(mut_vec)
  ## identical sequences → no event
  if (identical(paste0(ref_vec, collapse = ""), paste0(mut_vec, collapse = ""))) return(ret())
  len_ref <- length(ref_vec)
  len_mut <- length(mut_vec)

  ## =============================================================================================
  ## SUBSTITUTION (SNV / MNV) — NO PADDING
  ## =============================================================================================
  if (len_ref == len_mut) {
    n <- len_ref
    ## leftmost difference
    left_same <- 0L
    for (i in seq_len(n)) {
      if (ref_vec[i] != mut_vec[i]) break
      left_same <- i
    }
    ## rightmost difference without crossing left boundary
    right_same <- 0L
    max_right <- n - left_same
    if (max_right > 0L) {
      for (i in seq_len(max_right)) {
        if (ref_vec[n - i + 1L] != mut_vec[n - i + 1L]) break
        right_same <- i
      }
    }
    core_start <- left_same + 1L
    core_end   <- n - right_same
    if (core_start > core_end) return(ret())
    Ref_Bases <- ref_vec[core_start:core_end]
    Mut_Bases <- mut_vec[core_start:core_end]
    my_ref <- paste0(Ref_Bases, collapse = "")
    my_mut <- paste0(Mut_Bases, collapse = "")
    ## required by policy
    if (my_ref == my_mut) return(ret())
    return(ret(
      Status = TRUE,
      Data = list(
        Class      = "substitution",
        Change     = paste(my_ref, my_mut, sep = "."),
        Ref_Bases  = Ref_Bases,
        Mut_Bases  = Mut_Bases,
        Core_Start = core_start
      )
    ))
  }

  ## =============================================================================================
  ## INDEL — REQUIRE 1-PADDING ON BOTH ENDS (AND INCLUDE THEM)
  ## =============================================================================================
  min_len <- min(len_ref, len_mut)

  ## ---- 5' identical padding ----
  left_pad <- 0L
  for (i in seq_len(min_len)) {
    if (ref_vec[i] != mut_vec[i]) break
    left_pad <- i
  }
  if (left_pad < 1L) return(ret())

  ## ---- 3' identical padding ----
  right_pad <- 0L
  max_right <- min(len_ref - left_pad, len_mut - left_pad)
  if (max_right < 1L) return(ret())

  for (i in seq_len(max_right)) {
    if (ref_vec[len_ref - i + 1L] != mut_vec[len_mut - i + 1L]) break
    right_pad <- i
  }
  if (right_pad < 1L) return(ret())

  ## ---- include exactly one pad base on each end ----
  start_idx_ref <- left_pad
  end_idx_ref   <- len_ref - right_pad + 1L
  start_idx_mut <- left_pad
  end_idx_mut   <- len_mut - right_pad + 1L
  if (start_idx_ref < 1L || end_idx_ref > len_ref || start_idx_ref > end_idx_ref) return(ret())
  if (start_idx_mut < 1L || end_idx_mut > len_mut || start_idx_mut > end_idx_mut) return(ret())
  Ref_Bases <- ref_vec[start_idx_ref:end_idx_ref]
  Mut_Bases <- mut_vec[start_idx_mut:end_idx_mut]
  my_ref <- paste0(Ref_Bases, collapse = "")
  my_mut <- paste0(Mut_Bases, collapse = "")
  if (my_ref == my_mut) return(ret())
  Class <- if (len_mut > len_ref) "insertion" else "deletion"
  ret(
    Status = TRUE,
    Data = list(
      Class      = Class,
      Change     = paste(my_ref, my_mut, sep = "."),
      Ref_Bases  = Ref_Bases,
      Mut_Bases  = Mut_Bases,
      Core_Start = start_idx_ref
    )
  )
}


##################################################################################################
## compute_orf_lengths_stop_anchored
##
## PURPOSE
##   Compute ORF length metrics for reference and mutant mRNA sequences
##   in a STOP-ANCHORED, FRAME-CORRECT manner, using EXPLICIT CDS
##   boundaries provided in REF mRNA local coordinates.
##
##   This function is PURELY STRUCTURAL:
##     - NO variant classification
##     - NO splice logic
##     - NO biological interpretation
##
##   It reports three ORF-length metrics:
##
##     1) WT_orf
##        Amino-acid length of the REFERENCE ORF from canonical CDS start
##        (First_Base_of_CDS) to the canonical stop codon
##        (Last_Base_of_CDS).
##
##     2) Mut_orf5
##        Amino-acid length of the MUTANT ORF starting from the SAME
##        canonical CDS start (First_Base_of_CDS), translated in frame 0,
##        until the FIRST stop codon encountered.
##
##        If no stop is encountered, length extends to the end of the
##        translated sequence.
##
##     3) Mut_orf3
##        Amino-acid length of the LONGEST possible ORF in the MUTANT
##        sequence that TERMINATES at the CANONICAL stop codon.
##
##        Key properties:
##          - The canonical stop codon is defined using
##            Last_Base_of_CDS from REF.
##          - Frame is derived from the canonical stop position in MUT.
##          - Translation is STOP-ANCHORED, not START-ANCHORED.
##          - The ORF start is defined as the LAST 'M' (start codon)
##            strictly between:
##              (previous stop, canonical stop).
##
##        Special cases:
##          - Canonical stop lost            → Mut_orf3 = 0
##          - No start codon found upstream  → Mut_orf3 = 0
##
##
## ===================== INPUT =====================
##
##   mut_mRNA : character vector
##     Mutant mRNA bases (5' → 3'), transcript-local coordinates.
##
##   ref_mRNA : character vector
##     Reference mRNA bases (5' → 3'), transcript-local coordinates.
##
##   First_Base_of_CDS : integer scalar
##     Transcript-local coordinate (1-based) of the first CDS base.
##
##   Last_Base_of_CDS : integer scalar
##     Transcript-local coordinate (1-based, inclusive) of the last CDS base
##     (i.e., last base of the canonical stop codon).
##
##   consts : list
##     Must contain CodonTable usable by translate_dna_frame0_vec().
##
##
## ===================== OUTPUT =====================
##
##   list(
##     Status  = logical scalar,
##     Message = character scalar,
##     Data    = list
##   )
##
##   Data fields (valid only if Status == TRUE):
##
##     WT_orf   : integer scalar
##     Mut_orf5 : integer scalar
##     Mut_orf3 : integer scalar
##
##
## ===================== FAILURE POLICY =====================
##
##   - NEVER throws
##   - Return shape invariant
##   - Message empty on success, fixed default string on failure
##
##################################################################################################
compute_orf_lengths_stop_anchored <- function(mut_mRNA,ref_mRNA,First_Base_of_CDS,Last_Base_of_CDS,consts) {

  default_msg <- "in compute_orf_lengths_stop_anchored"
  ret <- function(Status = FALSE, Data = list())
    list(Status = Status, Message = if (Status) "" else default_msg, Data = Data)

  ## ---------------- Validation ----------------
  if (!is.character(mut_mRNA) || !is.character(ref_mRNA)) return(ret())
  if (any(nchar(mut_mRNA) != 1L) || any(nchar(ref_mRNA) != 1L)) return(ret())

  len_ref <- length(ref_mRNA)
  len_mut <- length(mut_mRNA)
  min_len <- min(len_ref, len_mut)

  if (!is.finite(First_Base_of_CDS) || !is.finite(Last_Base_of_CDS)) return(ret())
  if (First_Base_of_CDS < 1L ||
      Last_Base_of_CDS < First_Base_of_CDS ||
      Last_Base_of_CDS > len_ref) return(ret())
  if (((Last_Base_of_CDS - First_Base_of_CDS + 1L) %% 3L) != 0L) return(ret())

  ## ---------------- WT_orf ----------------
  WT_orf <- (Last_Base_of_CDS - First_Base_of_CDS + 1L) %/% 3L - 1L
  if (WT_orf < 0L) return(ret())

  ## ---------------- Mut_orf5 (start-anchored) ----------------
  mut_ORF5_nt <- mut_mRNA[First_Base_of_CDS:len_mut]
  mut_tr5 <- translate_dna_frame0_vec(mut_ORF5_nt, consts)
  if (!isTRUE(mut_tr5$Status)) return(ret())

  mut_pep5 <- mut_tr5$Data
  stop_pos5 <- which(mut_pep5 == "*")
  Mut_orf5 <- if (length(stop_pos5)) min(stop_pos5) - 1L else length(mut_pep5)

  ## ---------------- Right-identical suffix ----------------
  Left_Identical <- 0L
  for (i in seq_len(min_len)) {
    if (mut_mRNA[i] != ref_mRNA[i]) break
    Left_Identical <- i
  }

  Right_Identical <- 0L
  max_right <- min_len - Left_Identical
  if (max_right > 0L) {
    for (i in seq_len(max_right)) {
      if (ref_mRNA[len_ref - i + 1L] != mut_mRNA[len_mut - i + 1L]) break
      Right_Identical <- i
    }
  }

  ## ---------------- Canonical stop preservation ----------------
  ref_UTR3_len <- len_ref - Last_Base_of_CDS
  canonical_preserved <- isTRUE(Right_Identical >= (ref_UTR3_len + 3L))

  Mut_orf3 <- 0L

  if (canonical_preserved) {

    ## Canonical stop NT start in REF
    ref_stop_nt_start <- Last_Base_of_CDS - 2L

    ## Map into MUT using length delta
    mut_stop_nt_start <- ref_stop_nt_start + (len_mut - len_ref)

    if (mut_stop_nt_start >= 1L && mut_stop_nt_start <= len_mut - 2L) {

      ## Derive frame from stop codon
      stop_frame <- (mut_stop_nt_start - 1L) %% 3L
      mut_seq_frame <- mut_mRNA[(stop_frame + 1L):len_mut]

      mut_trF <- translate_dna_frame0_vec(mut_seq_frame, consts)
      if (isTRUE(mut_trF$Status)) {

        pepF <- mut_trF$Data
        stop_pos <- which(pepF == "*")

        if (length(stop_pos)) {
          last_stop <- max(stop_pos)
          prev_stop <- if (length(stop_pos) >= 2L)
            stop_pos[length(stop_pos) - 1L] else -9L

          if (last_stop > prev_stop + 1L) {
            M_pos <- which(
              pepF == "M" &
                seq_along(pepF) > prev_stop &
                seq_along(pepF) < last_stop
            )
            if (length(M_pos)) {
              best_M <- min(M_pos)   ## earliest gives longest ORF
              Mut_orf3 <- last_stop - best_M
            }
          }
        }
      }
    }
  }

  ## ---------------- Success ----------------
  ret(
    Status = TRUE,
    Data = list(
      WT_orf   = as.integer(WT_orf),
      Mut_orf5 = as.integer(Mut_orf5),
      Mut_orf3 = as.integer(Mut_orf3)
    )
  )
}


##################################################################################################
## analyze_translation_effect_indel_after_start_codon_with_UTR
##
## PURPOSE
##   Analyze the amino-acid–level consequence of an indel after splicing,
##   assuming the alteration occurs downstream of the canonical start codon
##   and may extend into the 3' UTR.
##
##   This function compares reference and mutant mRNA sequences, derives
##   open reading frames (ORFs), translates them assuming frame 0 at the
##   annotated CDS start, and classifies the resulting peptide-level effect.
##
##   This is a MID-LEVEL TRANSLATION INTERPRETATION utility:
##     - splicing has already been resolved upstream
##     - no splice-site analysis is performed here
##     - UTR5 events and start-codon–proximal events are explicitly excluded
##       and deferred to other analysis stages
##
##
## ===================== ORF LENGTH EXTENSION =====================
##
##   This function now ALSO reports transcript-level ORF capacities:
##
##     WT_orf
##       Length (AA) of the reference ORF from canonical start to canonical stop.
##
##     Mut_orf3
##       Length (AA) from canonical start to canonical stop in the mutant.
##       0 if canonical stop is lost or unreachable.
##
##     Mut_orf5
##       Length (AA) of the longest mutant ORF that terminates at ANY stop
##       upstream of (or at) the canonical stop.
##
##   ORF lengths are computed ONCE at function entry and appended to all
##   successful returns. Classification logic is unchanged.
##
##
## ===================== OUTPUT =====================
##
##   Data fields (valid only if Status == TRUE):
##
##     Class
##     AAChange
##     Epitope
##     AAPos
##     WT_orf
##     Mut_orf3
##     Mut_orf5
##
##################################################################################################
analyze_translation_effect_indel_after_start_codon_with_UTR <- function(mut_mRNA,ref_mRNA,First_Base_of_CDS,Last_Base_of_CDS,consts) {
  default_msg <- "in analyze_translation_effect_indel_after_start_codon_with_UTR"
  ret <- function(Status = FALSE, Data = list()) list(Status = Status, Message = if (Status) "" else default_msg, Data = Data)

  ## ---------------- Validation ----------------
  len_ref <- length(ref_mRNA)
  len_mut <- length(mut_mRNA)
  min_len <- min(len_ref, len_mut)

  if (len_ref < 100L) return(ret())
  if (is.na(First_Base_of_CDS) || First_Base_of_CDS < 1L || First_Base_of_CDS > len_ref) return(ret())
  if (is.na(Last_Base_of_CDS) || Last_Base_of_CDS < First_Base_of_CDS || Last_Base_of_CDS > len_ref || ((Last_Base_of_CDS - First_Base_of_CDS + 1L) %% 3L) != 0L) return(ret())

  ## ---------------- ORF length computation ----------------
  orf <- compute_orf_lengths_stop_anchored(mut_mRNA,ref_mRNA,First_Base_of_CDS,Last_Base_of_CDS,consts)
  if (!isTRUE(orf$Status)) return(ret())

  WT_orf   <- orf$Data$WT_orf
  Mut_orf3 <- orf$Data$Mut_orf3
  Mut_orf5 <- orf$Data$Mut_orf5

  ## ---------------- ORF essentially lost ----------------
  if (len_mut <= First_Base_of_CDS + 5L) {
    return(ret(
      Status = TRUE,
      Data = list(
        Class    = "orf_loss",
        AAChange = "orf_loss",
        Epitope  = "no_epitope",
        AAPos    = -9L,
        WT_orf   = WT_orf,
        Mut_orf3 = Mut_orf3,
        Mut_orf5 = Mut_orf5
      )
    ))
  }

  ## ---------------- First nucleotide difference ----------------
  Left_Identical <- 0L
  for (i in seq_len(min_len)) {
    if (mut_mRNA[i] != ref_mRNA[i]) break
    Left_Identical <- i
  }

  ## CDS-start–proximal or upstream → defer
  if (Left_Identical <= First_Base_of_CDS + 5L) return(ret())

  ## ---------------- Pure UTR3 ----------------
  if (Left_Identical >= Last_Base_of_CDS) {
    return(ret(
      Status = TRUE,
      Data = list(
        Class    = "utr3",
        AAChange = "neutral",
        Epitope  = "no_epitope",
        AAPos    = -9L,
        WT_orf   = WT_orf,
        Mut_orf3 = Mut_orf3,
        Mut_orf5 = Mut_orf5
      )
    ))
  }

  ## ---------------- ORF translation ----------------
  ref_ORF <- ref_mRNA[First_Base_of_CDS:len_ref]
  mut_ORF <- mut_mRNA[First_Base_of_CDS:len_mut]

  ref_tr <- translate_dna_frame0_vec(ref_ORF, consts)
  mut_tr <- translate_dna_frame0_vec(mut_ORF, consts)
  if (!isTRUE(ref_tr$Status) || !isTRUE(mut_tr$Status)) return(ret())

  ref_pep <- ref_tr$Data
  mut_pep <- mut_tr$Data

  if (length(ref_pep) < 3L || length(mut_pep) < 3L) {
    return(ret(
      Status = TRUE,
      Data = list(
        Class    = "orf_loss",
        AAChange = "orf_loss",
        Epitope  = "no_epitope",
        AAPos    = -9L,
        WT_orf   = WT_orf,
        Mut_orf3 = Mut_orf3,
        Mut_orf5 = Mut_orf5
      )
    ))
  }

  ## ---------------- Canonical stop validation ----------------
  ref_stop_pos <- which(ref_pep == "*")
  if (!length(ref_stop_pos)) return(ret())

  annotated_len <- (Last_Base_of_CDS - First_Base_of_CDS + 1L) %/% 3L
  if (min(ref_stop_pos) < annotated_len) return(ret())

  ## ---------------- AA-level divergence anchor ----------------
  aa_left <- 0L
  max_cmp <- min(length(ref_pep), length(mut_pep))
  for (i in seq_len(max_cmp)) {
    if (ref_pep[i] != mut_pep[i]) break
    aa_left <- i
  }
  if (aa_left < 1L) return(ret())

  AAPos <- aa_left
  in_frame <- ((length(mut_ORF) - length(ref_ORF)) %% 3L) == 0L
  stops_mut <- which(mut_pep == "*")

  ## ---------------- Frameshift ----------------
  if (!in_frame) {
    stop_pos <- if (length(stops_mut) && any(stops_mut >= AAPos))
      min(stops_mut[stops_mut >= AAPos]) else length(mut_pep)

    Epitope <- paste0(mut_pep[AAPos:stop_pos], collapse = "")
    fs_len  <- sum(mut_pep[AAPos:stop_pos] != "*")

    return(ret(
      Status = TRUE,
      Data = list(
        Class    = "frameshift",
        AAChange = paste0(AAPos, ".", ref_pep[AAPos], ".fs_", fs_len, "_AA"),
        Epitope  = Epitope,
        AAPos    = AAPos,
        WT_orf   = WT_orf,
        Mut_orf3 = Mut_orf3,
        Mut_orf5 = Mut_orf5
      )
    ))
  }

  ## ---------------- In-frame core comparison ----------------
  aa_right <- 0L
  for (i in seq_len(max_cmp - aa_left)) {
    if (ref_pep[length(ref_pep) - i + 1L] != mut_pep[length(mut_pep) - i + 1L]) break
    aa_right <- i
  }

  core_ref <- ref_pep[AAPos:(length(ref_pep) - aa_right)]
  core_mut <- mut_pep[AAPos:(length(mut_pep) - aa_right)]

  if (AAPos + length(core_ref) <= length(ref_pep)) {
    core_ref <- c(core_ref, ref_pep[AAPos + length(core_ref)])
    core_mut <- c(core_mut, mut_pep[AAPos + length(core_mut)])
  }

  ref_has_stop <- "*" %in% core_ref
  mut_has_stop <- "*" %in% core_mut

  ## ---------------- Stop gain ----------------
  if (!ref_has_stop && mut_has_stop) {
    stop_pos <- min(which(mut_pep == "*"))
    ref_aa   <- if (stop_pos <= length(ref_pep)) ref_pep[stop_pos] else "X"

    return(ret(
      Status = TRUE,
      Data = list(
        Class    = "stop_gain",
        AAChange = paste0(stop_pos, ".", ref_aa, ">*"),
        Epitope  = "no_epitope",
        AAPos    = stop_pos,
        WT_orf   = WT_orf,
        Mut_orf3 = Mut_orf3,
        Mut_orf5 = Mut_orf5
      )
    ))
  }

  ## ---------------- Stop loss ----------------
  if (ref_has_stop && !mut_has_stop) {
    stop_pos <- if (length(stops_mut) && any(stops_mut >= AAPos))
      min(stops_mut[stops_mut >= AAPos]) else length(mut_pep)

    Epitope <- paste0(mut_pep[AAPos:stop_pos], collapse = "")
    fs_len  <- sum(mut_pep[AAPos:stop_pos] != "*")

    return(ret(
      Status = TRUE,
      Data = list(
        Class    = "stop_loss",
        AAChange = paste0(AAPos, ".", ref_pep[AAPos], ".fs_", fs_len, "_AA"),
        Epitope  = Epitope,
        AAPos    = AAPos,
        WT_orf   = WT_orf,
        Mut_orf3 = Mut_orf3,
        Mut_orf5 = Mut_orf5
      )
    ))
  }

  ## ---------------- Neutral / missense / insertion / deletion ----------------
  ref_str <- paste0(core_ref, collapse = "")
  mut_str <- paste0(core_mut, collapse = "")

  if (ref_str == mut_str) {
    Class <- "neutral"
    AAChange <- "neutral"
  } else if (nchar(mut_str) > nchar(ref_str)) {
    Class <- "insertion"
    AAChange <- paste0(AAPos, ".", ref_str, ">", mut_str)
  } else if (nchar(mut_str) < nchar(ref_str)) {
    Class <- "deletion"
    AAChange <- paste0(AAPos, ".", ref_str, ">", mut_str)
  } else {
    Class <- "missense"
    AAChange <- paste0(AAPos, ".", ref_str, ">", mut_str)
  }

  ret(
    Status = TRUE,
    Data = list(
      Class    = Class,
      AAChange = AAChange,
      Epitope  = "no_epitope",
      AAPos    = AAPos,
      WT_orf   = WT_orf,
      Mut_orf3 = Mut_orf3,
      Mut_orf5 = Mut_orf5
    )
  )
}


##################################################################################################
## analyze_translation_effect_in_neo_exon
##
## PURPOSE
##   Analyze the amino-acid–level translation effect of an indel that
##   results in the creation of a neo-exon, after splicing has already
##   been applied upstream.
##
##   The function constructs a mutant mRNA by concatenating:
##
##     [ upstream transcript (non-Buffer) ] +
##     [ neo-exon ] +
##     [ downstream transcript (non-Buffer) ]
##
##   It then compares this mutant mRNA against a reference mRNA derived
##   from the provided mRNA annotation table (excluding "Buffer") and
##   determines peptide-level consequences starting from the canonical
##   CDS start.
##
##   This function is a MID-LEVEL TRANSLATION INTERPRETATION utility.
##   It assumes:
##     - splicing decisions have already been made upstream
##     - the neo-exon sequence is finalized and transcript-oriented (5'→3')
##     - translation starts at the annotated CDS start
##
##   This function explicitly does NOT analyze:
##     - indels in UTR5
##     - indels overlapping or proximal to the start codon
##     - splice-site creation or loss
##
##
## ===================== REFERENCE mRNA DEFINITION =====================
##
##   local_mRNA = mRNA[mRNA$Region != "Buffer", ]
##   ref_mRNA   = local_mRNA$T_Base
##
##   CDS boundaries in this local coordinate system:
##
##     First_Base_of_CDS = (# UTR5 bases) + 1
##     Last_Base_of_CDS  = (# UTR5 + # CDS bases)
##
##
## ===================== MUTANT mRNA CONSTRUCTION =====================
##
##   Buffer bases in upstream/downstream are DROPPED.
##   Buffer context is NOT translated.
##
##
## ===================== ORF METRICS =====================
##
##   WT_orf
##     - Length (AA) of the reference ORF
##
##   Mut_orf5
##     - Length (AA) of the mutant ORF starting from the canonical START
##       codon to the first STOP codon (novel or canonical)
##     - 0 if canonical START is lost
##
##   Mut_orf3
##     - Length (AA) of the LONGEST mutant ORF that terminates at the
##       canonical STOP codon
##     - 0 if no valid upstream START exists
##     - 0 if canonical-start ORF already accounts for this length
##
##
## ===================== OUTPUT =====================
##
##   list(
##     Status  = logical scalar,
##     Message = character scalar,
##     Data    = list
##   )
##
##   Data fields (valid only if Status == TRUE):
##
##     AAPos     : integer scalar
##     Class     : character scalar
##     AAChange  : character scalar
##     Epitope   : character scalar
##     WT_orf    : integer scalar
##     Mut_orf5  : integer scalar
##     Mut_orf3  : integer scalar
##
##
## ===================== FAILURE POLICY =====================
##
##   - NEVER throws
##   - Return shape is invariant
##   - Message is propagated from inner calls when applicable
##
##################################################################################################
analyze_translation_effect_in_neo_exon <- function(the_upstream,neo_exon,the_downstream,mRNA,consts) {
  default_msg <- "in analyze_translation_effect_in_neo_exon"

  ret <- function(Status = FALSE, Message = default_msg, Data = list()) {
    list(
      Status  = Status,
      Message = if (Status) "" else Message,
      Data    = Data
    )
  }

  ## ---------------- Validation ----------------
  if (!is.data.frame(mRNA) || nrow(mRNA) == 0L) return(ret())
  if (is.null(mRNA$T_Base) || is.null(mRNA$Region)) return(ret())
  if (!is.character(neo_exon)) return(ret())
  if (length(neo_exon) > 0L && any(nchar(neo_exon) != 1L)) return(ret())

  ## ------------------------------------------------------------
  ## Construct reference mRNA (drop Buffer)
  ## ------------------------------------------------------------
  local_mRNA <- mRNA[mRNA$Region != "Buffer", , drop = FALSE]
  if (!nrow(local_mRNA)) return(ret())

  ref_mRNA <- as.character(local_mRNA$T_Base)

  First_Base_of_CDS <- sum(local_mRNA$Region == "UTR5") + 1L
  Last_Base_of_CDS  <- sum(local_mRNA$Region %in% c("UTR5", "CDS"))

  if (First_Base_of_CDS < 1L ||
      Last_Base_of_CDS < First_Base_of_CDS ||
      Last_Base_of_CDS > length(ref_mRNA) ||
      ((Last_Base_of_CDS - First_Base_of_CDS + 1L) %% 3L) != 0L)
    return(ret())

  ## ------------------------------------------------------------
  ## Construct mutant mRNA (drop Buffer)
  ## ------------------------------------------------------------
  mut_mRNA <- character(0)

  if (is.data.frame(the_upstream) && nrow(the_upstream)) {
    if (is.null(the_upstream$T_Base) || is.null(the_upstream$Region)) return(ret())
    mut_mRNA <- c(mut_mRNA, the_upstream$T_Base[the_upstream$Region != "Buffer"])
  }

  mut_mRNA <- c(mut_mRNA, neo_exon)

  if (is.data.frame(the_downstream) && nrow(the_downstream)) {
    if (is.null(the_downstream$T_Base) || is.null(the_downstream$Region)) return(ret())
    mut_mRNA <- c(mut_mRNA, the_downstream$T_Base[the_downstream$Region != "Buffer"])
  }

  ## ------------------------------------------------------------
  ## Neutral / no-change shortcut
  ## ------------------------------------------------------------
  if (length(mut_mRNA) == length(ref_mRNA) &&
      identical(paste0(mut_mRNA, collapse = ""),
                paste0(ref_mRNA, collapse = ""))) {

    ## compute WT ORF length
    cds_len <- Last_Base_of_CDS - First_Base_of_CDS + 1L
    WT_orf  <- cds_len %/% 3L
    Data = list(
      AAPos    = -9L,
      Class    = "neutral",
      AAChange = "neutral",
      Epitope  = "no_epitope",
      WT_orf   = WT_orf,
      Mut_orf5 = WT_orf,
      Mut_orf3 = 0L
    )
    return(ret(Status = TRUE,Data = Data))
  }

  ## ------------------------------------------------------------
  ## Identify leftmost difference
  ## ------------------------------------------------------------
  min_len <- min(length(ref_mRNA), length(mut_mRNA))
  Left_Identical <- 0L
  for (i in seq_len(min_len)) {
    if (ref_mRNA[i] != mut_mRNA[i]) break
    Left_Identical <- i
  }

  ## ------------------------------------------------------------
  ## CASE: indel sufficiently downstream of CDS start
  ## ------------------------------------------------------------
  if (Left_Identical > First_Base_of_CDS + 5L) {

    rlt <- analyze_translation_effect_indel_after_start_codon_with_UTR(mut_mRNA,ref_mRNA,First_Base_of_CDS,Last_Base_of_CDS,consts)

    if (!isTRUE(rlt$Status))
      return(ret(Status = FALSE, Message = rlt$Message))

    return(ret(
      Status = TRUE,
      Data = list(
        AAPos     = rlt$Data$AAPos,
        Class     = rlt$Data$Class,
        AAChange  = rlt$Data$AAChange,
        Epitope   = rlt$Data$Epitope,
        WT_orf    = rlt$Data$WT_orf,
        Mut_orf5  = rlt$Data$Mut_orf5,
        Mut_orf3  = rlt$Data$Mut_orf3
      )
    ))
  }

  ## ------------------------------------------------------------
  ## CASE: UTR5 or start-codon–proximal; big TODO; some fake defaults
  ## ------------------------------------------------------------
  ret(
    Status  = TRUE,
    Message = "UTR5/start-codon cases not implemented",
    Data    = list(
      AAPos     = 100L,
      Class     = "unresolved",
      AAChange  = "unresolved",
      Epitope   = "no_epitope",
      WT_orf    = 999L,
      Mut_orf5  = 0L,
      Mut_orf3  = 0L
    )
  )
}


##################################################################################################
## convert_indel_notation_from_genomic_to_transcript
##
## PURPOSE
##   Convert a genomic substitution (SNV / MNV) or deletion (defined in
##   genomic coordinates and genomic strand orientation) into a
##   transcript-local (pre-mRNA) indel representation.
##
##   This function exists primarily to:
##     - resolve minus-strand issues
##     - map genomic positions to transcript-local coordinates
##     - optionally validate that the provided REF bases match the host
##       transcript sequence
##
##   INSERTIONS ARE EXPLICITLY DISALLOWED at this layer to avoid ambiguity
##   and inconsistent downstream representations.
##
##
## ===================== CRITICAL POLICY (STRICT) =====================
##
##   ONLY events satisfying:
##
##       indel_G_start <= indel_G_end
##
##   are allowed.
##
##   This corresponds to:
##     - SNV                 (start == end)
##     - MNV                 (end - start + 1 > 1)
##     - deletion
##     - complex replacement
##
##   Events with:
##
##       indel_G_start > indel_G_end
##
##   (i.e. insertions) are REJECTED with Status = FALSE.
##
##
## ===================== ASSUMPTIONS =====================
##
##   - host_gene_region_data_frame is ALREADY ordered in transcript
##     orientation (5' → 3').
##
##   - host_gene_region_data_frame contains at least:
##       G_Pos        : genomic coordinate (1-based)
##       T_Pos_Local  : transcript-local coordinate (1-based)
##       T_Base       : transcript-oriented base (A/C/G/T/N)
##
##   - indel_G_refbase and indel_G_mutbase are provided in genomic
##     forward-strand notation.
##
##
## ===================== STATUS SEMANTICS =====================
##
##   Status = TRUE
##     → event is fully contained within the host transcript region,
##       transcript-local coordinates were successfully derived, and
##       (if requested) reference bases are consistent.
##
##   Status = FALSE
##     → ANY of the following:
##         - insertion event (Start_G > End_G)
##         - chromosome mismatch
##         - containment failure
##         - reference-base mismatch
##         - malformed input
##         - REF length inconsistent with genomic span
##
##
## ===================== OUTPUT =====================
##
##   list(
##     Status  = logical scalar,
##     Message = character scalar,
##     Data    = list
##   )
##
##   Data fields (valid only if Status == TRUE):
##
##     T_indel_chr     : character scalar
##     T_indel_start   : integer scalar (transcript-local, 1-based)
##     T_indel_end     : integer scalar (transcript-local, inclusive)
##     T_indel_refbase : character vector
##     T_indel_mutbase : character vector
##
##################################################################################################
convert_indel_notation_from_genomic_to_transcript <- function(
  indel_chr,
  indel_G_start,
  indel_G_end,
  indel_G_refbase,
  indel_G_mutbase,
  host_gene_region_data_frame,
  host_gene_chr,
  host_gene_strand,
  to_perform_reference_check = FALSE
) {
  default_msg <- "in convert_indel_notation_from_genomic_to_transcript"
  ret <- function(Status = FALSE, Data = list()) {
    list(
      Status  = Status,
      Message = if (Status) "" else default_msg,
      Data    = Data
    )
  }

  ## ---------------- Validation ----------------
  if (!is.data.frame(host_gene_region_data_frame)) return(ret())
  required_cols <- c("G_Pos", "T_Pos_Local", "T_Base")
  if (!all(required_cols %in% names(host_gene_region_data_frame))) return(ret())

  if (!is.character(indel_chr) || !is.character(host_gene_chr)) return(ret())
  if (!is.character(indel_G_refbase) || !is.character(indel_G_mutbase)) return(ret())
  if (!host_gene_strand %in% c("+", "-")) return(ret())
  if (indel_chr != host_gene_chr) return(ret())

  indel_G_start <- as.integer(indel_G_start)
  indel_G_end   <- as.integer(indel_G_end)
  if (is.na(indel_G_start) || is.na(indel_G_end)) return(ret())

  ## ---------------- STRICT POLICY ----------------
  ## allow SNV / MNV / deletion (Start_G <= End_G), disallow insertions
  if (indel_G_start > indel_G_end) return(ret())

  ## ---------------- Span / REF-length consistency ----------------
  span <- indel_G_end - indel_G_start + 1L
  if (span < 1L) return(ret())
  if (nchar(indel_G_refbase) != span) return(ret())

  ## ---------------- Strand correction ----------------
  T_refbase <- indel_G_refbase
  T_mutbase <- indel_G_mutbase
  if (host_gene_strand == "-") {
    T_refbase <- revcomp_string(indel_G_refbase)
    T_mutbase <- revcomp_string(indel_G_mutbase)
  }

  ## ---------------- Genomic containment ----------------
  indel_G_Pos <- indel_G_start:indel_G_end
  if (!all(indel_G_Pos %in% host_gene_region_data_frame$G_Pos)) return(ret())

  idx <- host_gene_region_data_frame$G_Pos %in% indel_G_Pos
  genomic_ref_in_transcript <- paste0(
    host_gene_region_data_frame$T_Base[idx],
    collapse = ""
  )

  ## ---------------- Reference validation ----------------
  if (to_perform_reference_check && genomic_ref_in_transcript != T_refbase)
    return(ret())

  tp <- host_gene_region_data_frame$T_Pos_Local[idx]

  ## ---------------- Success ----------------
  ret(
    Status = TRUE,
    Data = list(
      T_indel_chr     = indel_chr,
      T_indel_start   = min(tp),
      T_indel_end     = max(tp),
      T_indel_refbase = if (nzchar(T_refbase))
        strsplit(T_refbase, "", fixed = TRUE)[[1]]
      else
        character(0),
      T_indel_mutbase = if (nzchar(T_mutbase))
        strsplit(T_mutbase, "", fixed = TRUE)[[1]]
      else
        character(0)
    )
  )
}

##################################################################################################
## check_mutant_mRNA_for_indels_with_splice_consideration
##
## PURPOSE
##   Perform unified splice-impact analysis for substitutions (SNVs/MNVs)
##   and deletions overlapping exonic or intronic regions, under a STRICT
##   project policy that DISALLOWS insertions.
##
##   This function evaluates BOTH:
##     - loss of canonical splice site affinity (abolish)
##     - creation of novel splice site affinity (establish)
##
##   Region-specific behavior:
##     Exon   → donor site analysis (neo-donor / donor abolish)
##     Intron → acceptor site analysis (neo-acceptor / acceptor abolish)
##
##   Translation consequences of neo-exons are evaluated using
##   analyze_translation_effect_in_neo_exon().
##
##
## ===================== CRITICAL POLICY (STRICT) =====================
##
##   ONLY events satisfying:
##
##       indels$Start_G <= indels$End_G
##
##   are eligible for analysis in THIS FUNCTION.
##
##
## ===================== STRAND / EXON ORDER SEMANTICS =====================
##
##   mRNA$G_ExonNo is ALWAYS in GENOMIC order.
##   Upstream / downstream context is STRAND-AWARE.
##
##
## ===================== OUTPUT =====================
##
##   Data columns:
##     Marker
##     SearchFor
##     Region
##     Class
##     AAChange
##     Affinity_WT
##     Affinity_Mut
##     Splice_Core_WT
##     Splice_Core_Mut
##     WT_orf
##     Mut_orf5
##     Mut_orf3
##
##################################################################################################
check_mutant_mRNA_for_indels_with_splice_consideration <- function(indels, region_df, region_type, gene_chr, gene_strand, mRNA, consts) {
  default_msg <- "in check_mutant_mRNA_for_indels_with_splice_consideration"
  AFF_NA <- -9.9

  ret <- function(Status = FALSE, Message = default_msg, Data = data.frame())
    list(Status = Status, Message = if (Status) "" else Message, Data = Data)

  ## ---------------- Validation ----------------
  if (!region_type %in% c("Exon", "Intron")) return(ret())
  if (!is.data.frame(region_df) || nrow(region_df) == 0L)
    return(ret(Status = TRUE, Data = data.frame()))
  if (!all(c("T_Pos_Local","T_Base","Region","G_ExonNo") %in% names(region_df)))
    return(ret())
  if (!gene_strand %in% c("+","-")) return(ret())
  if (!is.data.frame(mRNA) || nrow(mRNA) == 0L) return(ret())

  required_indel_cols <- c("Chr","Start_G","End_G","Ref","Mut","Marker")
  if (!all(required_indel_cols %in% names(indels))) return(ret())

  ## ---------------- Motif selection ----------------
  if (region_type == "Exon") {
    motif     <- consts$Donor_Motif
    SearchFor <- "NeoDonor"
  } else {
    motif     <- consts$Acceptor_Motif
    SearchFor <- "NeoAcceptor"
  }

  ## ---------------- Output table ----------------
  out <- data.frame(
    Marker          = character(),
    SearchFor       = character(),
    Region          = character(),
    Class           = character(),
    AAChange        = character(),
    Affinity_WT     = numeric(),
    Affinity_Mut    = numeric(),
    Splice_Core_WT  = character(),
    Splice_Core_Mut = character(),
    WT_orf          = integer(),
    Mut_orf5        = integer(),
    Mut_orf3        = integer(),
    stringsAsFactors = FALSE
  )

  .append_row <- function(df, ...) rbind(df, data.frame(..., stringsAsFactors = FALSE))

  .max_score_at_pos <- function(pos, score, p) {
    if (!length(pos)) return(AFF_NA)
    x <- score[pos == p]
    if (!length(x)) return(AFF_NA)
    suppressWarnings(max(x, na.rm = TRUE))
  }

  .string_at_pos <- function(pos, str, p) {
    if (!length(pos)) return(" ")
    x <- str[pos == p]
    if (!length(x)) " " else x[1L]
  }

  ## ---------------- Canonical context ----------------
  core_tpos <- region_df$T_Pos_Local[region_df$Region == region_type]
  if (!length(core_tpos)) return(ret(Status = TRUE, Data = out))

  exon_nos <- mRNA$G_ExonNo[mRNA$G_ExonNo > 0]
  total_exon_number <- max(exon_nos)

  if (region_type == "Exon") {
    T_acceptor_AG <- min(core_tpos) - 1L
    T_donor_GT    <- max(core_tpos) + 1L
    Canonic_Pos   <- T_donor_GT
    exon_no_core  <- min(region_df$G_ExonNo[region_df$G_ExonNo > 0])
  } else {
    T_acceptor_AG <- max(core_tpos)
    Canonic_Pos   <- T_acceptor_AG
    exon3_tpos    <- region_df$T_Pos_Local[region_df$Region == "Exon_3"]
    end_of_exon3  <- max(exon3_tpos)
    exon_no_core  <- min(region_df$G_ExonNo[region_df$Region == "Exon_3" & region_df$G_ExonNo > 0])
  }

  ## ---------------- Strand-aware context ----------------
  if (gene_strand == "+") {
    up_range   <- if (exon_no_core > 1L) seq_len(exon_no_core - 1L) else integer(0)
    down_range <- if (exon_no_core < total_exon_number) seq.int(exon_no_core + 1L, total_exon_number) else integer(0)
  } else {
    up_range   <- if (exon_no_core < total_exon_number) seq.int(exon_no_core + 1L, total_exon_number) else integer(0)
    down_range <- if (exon_no_core > 1L) seq_len(exon_no_core - 1L) else integer(0)
  }

  the_upstream   <- mRNA[mRNA$G_ExonNo %in% up_range   & mRNA$Region != "Buffer", , drop = FALSE]
  the_downstream <- mRNA[mRNA$G_ExonNo %in% down_range & mRNA$Region != "Buffer", , drop = FALSE]

  ## ======================================================
  ## Iterate indels
  ## ======================================================
  for (i in seq_len(nrow(indels))) {

    Marker <- indels$Marker[i]
    if (indels$Chr[i] != gene_chr) next
    if (indels$Start_G[i] > indels$End_G[i]) next
    if (indels$Ref[i] == indels$Mut[i]) next

    indel_T <- convert_indel_notation_from_genomic_to_transcript(
      indels$Chr[i], indels$Start_G[i], indels$End_G[i],
      indels$Ref[i], indels$Mut[i],
      region_df, gene_chr, gene_strand,
      consts$Enforce_Refbase_Check
    )
    if (!isTRUE(indel_T$Status)) next

    mut_info <- make_mut_vec_full_in_local_coordinates(
      region_df,
      indel_T$Data$T_indel_start,
      indel_T$Data$T_indel_end,
      indel_T$Data$T_indel_refbase,
      indel_T$Data$T_indel_mutbase,
      consts$Enforce_Refbase_Check
    )
    if (!isTRUE(mut_info$Status)) return(ret(FALSE, mut_info$Message))

    ref_vec <- c(mut_info$Data$ref_up, mut_info$Data$ref_core, mut_info$Data$ref_down)
    mut_vec <- c(mut_info$Data$ref_up, mut_info$Data$alt_core, mut_info$Data$ref_down)

    hits <- scan_splice_motif_in_mut_and_wt_DNA_segment_with_sufficient_buffer_in_both_ends(ref_vec, mut_vec, motif)
	
    if (!isTRUE(hits$Status)) return(ret(FALSE, hits$Message))
	if(isTRUE(hits$Status) && !length(hits$Data$Ref_Hit_Pos) && !length(hits$Data$Mut_Hit_Pos)) next

    WT_Canonic_Score <- .max_score_at_pos(hits$Data$Ref_Hit_Pos, hits$Data$Ref_Hit_Score, Canonic_Pos)
    MUT_Global_Max <- if (length(hits$Data$Mut_Hit_Score)) max(hits$Data$Mut_Hit_Score) else AFF_NA

    WT_String <- .string_at_pos(hits$Data$Ref_Hit_Pos, hits$Data$Ref_Hit_String, Canonic_Pos)

    ## ---------------- Abolish ----------------
    if (WT_Canonic_Score > AFF_NA &&
        MUT_Global_Max < WT_Canonic_Score - motif$Delta_Abolish) {

      out <- .append_row(
        out,
        Marker          = Marker,
        SearchFor       = SearchFor,
        Region          = if (region_type == "Exon") "splice_donor" else "splice_acceptor",
        Class           = if (region_type == "Exon") "donor_abolish" else "acceptor_abolish",
        AAChange        = "exon_loss",
        Affinity_WT     = WT_Canonic_Score,
        Affinity_Mut    = MUT_Global_Max,
        Splice_Core_WT  = WT_String,
        Splice_Core_Mut = " ",
        WT_orf          = -9L,
        Mut_orf5        = -9L,
        Mut_orf3        = -9L
      )
      next
    }

    ## ---------------- Establish ----------------
    thr <- max(motif$Score_Cutoff, WT_Canonic_Score + motif$Delta_Enhance)
    good_mut_pos <- hits$Data$Mut_Hit_Pos[hits$Data$Mut_Hit_Score >= thr]

    for (p in good_mut_pos) {

      MUT_String <- .string_at_pos(
        hits$Data$Mut_Hit_Pos, hits$Data$Mut_Hit_String, p
      )

      if (region_type == "Exon") {
        if (p <= T_acceptor_AG + 1L) next
        neo_exon  <- mut_vec[(T_acceptor_AG + 1L):(p - 1L)]
        RegionTag <- "exon"
      } else {
        if (p >= end_of_exon3) next
        neo_exon  <- mut_vec[(p + 1L):end_of_exon3]
        RegionTag <- if (p <= T_acceptor_AG) "intron_retention" else "exon_truncation"
      }

      tr <- analyze_translation_effect_in_neo_exon(the_upstream, neo_exon, the_downstream, mRNA, consts)

      if (!isTRUE(tr$Status)) return(ret(FALSE, tr$Message))

      out <- .append_row(
        out,
        Marker          = Marker,
        SearchFor       = SearchFor,
        Region          = RegionTag,
        Class           = if (region_type == "Exon") "neo_donor" else "neo_acceptor",
        AAChange        = tr$Data$AAChange,
        Affinity_WT     = WT_Canonic_Score,
        Affinity_Mut    = .max_score_at_pos(hits$Data$Mut_Hit_Pos, hits$Data$Mut_Hit_Score, p),
        Splice_Core_WT  = WT_String,
        Splice_Core_Mut = MUT_String,
        WT_orf          = tr$Data$WT_orf,
        Mut_orf5        = tr$Data$Mut_orf5,
        Mut_orf3        = tr$Data$Mut_orf3
      )
    }
  }

  ret(Status = TRUE, Data = out)
}

##################################################################################################
## check_mutant_mRNA_for_indels_fully_inside_exon
##
## PURPOSE
##   Analyze small sequence variants that are fully contained within a
##   single exon, under the assumption of exon retention and NO splice
##   alteration.
##
##   Supported mutation classes include:
##     - SNVs
##     - MNVs
##     - small deletions / small replacements (project-defined)
##
##   For each qualifying mutation:
##     1) construct the mutant exon sequence
##     2) splice it into transcript context
##     3) evaluate amino-acid–level translation impact
##
##
## ===================== PROJECT-LEVEL MUTATION POLICY =====================
##
##   Although lower-level utilities in this codebase support general
##   insertions (T_indel_start > T_indel_end), THIS FUNCTION DOES NOT.
##
##   Mutation admissibility is governed EXCLUSIVELY by consts:
##
##     - consts$Max_Indel_Length
##     - consts$Enforce_Refbase_Check
##
##
## ===================== TRANSCRIPT CONTEXT POLICY (STRAND-AWARE) =====================
##
##   this_exon_number is a GENOMIC exon index.
##
##   mRNA$G_ExonNo is ALWAYS in GENOMIC order (>0 for real exons).
##
##   Upstream / downstream rules (EDGE-SAFE):
##
##     If strand == "+":
##       upstream   = exons 1 .. (this_exon_number - 1)
##       downstream = exons (this_exon_number + 1) .. total_exon_number
##
##     If strand == "-":
##       upstream   = exons (this_exon_number + 1) .. total_exon_number
##       downstream = exons 1 .. (this_exon_number - 1)
##
##
## ===================== OUTPUT =====================
##
##   Data columns:
##     Marker
##     SearchFor
##     Region
##     Class
##     AAChange
##     Affinity_WT
##     Affinity_Mut
##     WT_orf
##     Mut_orf5
##     Mut_orf3
##
##################################################################################################
check_mutant_mRNA_for_indels_fully_inside_exon <- function(indels,region_df,this_exon_number,gene_chr,gene_strand,mRNA,consts) {
  default_msg <- "in check_mutant_mRNA_for_indels_fully_inside_exon"
  AFF_NA <- -9.9

  ## ---------------- unified return helpers ----------------
  empty_df <- function() {
    data.frame(
      Marker       = character(),
      SearchFor    = character(),
      Region       = character(),
      Class        = character(),
      AAChange     = character(),
      Affinity_WT  = numeric(),
      Affinity_Mut = numeric(),
      WT_orf       = integer(),
      Mut_orf5     = integer(),
      Mut_orf3     = integer(),
      stringsAsFactors = FALSE
    )
  }

  fail_ret <- function(msg = default_msg) list(Status = FALSE, Message = msg, Data = empty_df())
  ok_ret <- function(df) list(Status = TRUE, Message = "", Data = df)

  ## ---------------- Validation ----------------
  if (!is.data.frame(mRNA) || nrow(mRNA) == 0L) return(fail_ret())
  if (!gene_strand %in% c("+","-")) return(fail_ret())
  need_cols <- c("Chr","Start_G","End_G","Ref","Mut","Marker")
  if (!is.data.frame(indels) || !all(need_cols %in% names(indels))) return(fail_ret())
  exon_nos <- mRNA$G_ExonNo[mRNA$G_ExonNo > 0]
  if (!length(exon_nos)) return(fail_ret())
  total_exon_number <- max(exon_nos)
  if (this_exon_number < 1L || this_exon_number > total_exon_number) return(fail_ret())

  ## ---------------- strand-aware exon ranges ----------------
  if (gene_strand == "+") {
    up_range   <- if (this_exon_number > 1L)
      seq_len(this_exon_number - 1L) else integer(0)
    down_range <- if (this_exon_number < total_exon_number)
      seq.int(this_exon_number + 1L, total_exon_number) else integer(0)
  } else {
    up_range   <- if (this_exon_number < total_exon_number) seq.int(this_exon_number + 1L, total_exon_number) else integer(0)
    down_range <- if (this_exon_number > 1L) seq_len(this_exon_number - 1L) else integer(0)
  }

  the_upstream <- if (length(up_range)) mRNA[mRNA$G_ExonNo %in% up_range, , drop = FALSE] else data.frame()
  the_downstream <- if (length(down_range)) mRNA[mRNA$G_ExonNo %in% down_range, , drop = FALSE] else data.frame()

  ## ---------------- extract exon core ----------------
  exon_only <- region_df[region_df$Region == "Exon", , drop = FALSE]
  if (nrow(exon_only) == 0L) return(ok_ret(empty_df()))

  exon_only <- exon_only[order(exon_only$T_Pos_Local), , drop = FALSE]
  exon_only$T_Pos_Local <- seq_len(nrow(exon_only))

  ## ---------------- output container ----------------
  out_rows <- list()
  out_n <- 0L

  ## ===================== iterate mutations =====================
  for (i in seq_len(nrow(indels))) {

    Marker <- indels$Marker[i]
    if (indels$Chr[i] != gene_chr) next
    if (is.na(indels$Ref[i]) || is.na(indels$Mut[i])) next
    if (indels$Ref[i] == indels$Mut[i]) next

    Start_G <- as.integer(indels$Start_G[i])
    End_G   <- as.integer(indels$End_G[i])
    if (is.na(Start_G) || is.na(End_G)) next
    if (Start_G > End_G) next

    if (nchar(indels$Ref[i]) >= consts$Max_Indel_Length) next
    if (nchar(indels$Mut[i]) >= consts$Max_Indel_Length) next

    span <- End_G - Start_G + 1L
    if (nchar(indels$Ref[i]) != span) next
    if (!all(Start_G:End_G %in% exon_only$G_Pos)) next

    indel_T <- convert_indel_notation_from_genomic_to_transcript(
      indels$Chr[i], Start_G, End_G,
      indels$Ref[i], indels$Mut[i],
      exon_only, gene_chr, gene_strand,
      consts$Enforce_Refbase_Check
    )
    if (!isTRUE(indel_T$Status))
      return(fail_ret(indel_T$Message))

    mut_info <- make_mut_vec_full_in_local_coordinates(
      exon_only,
      indel_T$Data$T_indel_start,
      indel_T$Data$T_indel_end,
      indel_T$Data$T_indel_refbase,
      indel_T$Data$T_indel_mutbase,
      consts$Enforce_Refbase_Check
    )
    if (!isTRUE(mut_info$Status))
      return(fail_ret(mut_info$Message))

    neo_exon <- c(
      mut_info$Data$ref_up,
      mut_info$Data$alt_core,
      mut_info$Data$ref_down
    )

    tr <- analyze_translation_effect_in_neo_exon(the_upstream, neo_exon, the_downstream, mRNA, consts)
    if (!isTRUE(tr$Status))
      return(fail_ret(tr$Message))

    out_n <- out_n + 1L
    out_rows[[out_n]] <- data.frame(
      Marker       = Marker,
      SearchFor    = "Coding",
      Region       = "exon",
      Class        = tr$Data$Class,
      AAChange     = tr$Data$AAChange,
      Affinity_WT  = AFF_NA,
      Affinity_Mut = AFF_NA,
      WT_orf       = tr$Data$WT_orf,
      Mut_orf5     = tr$Data$Mut_orf5,
      Mut_orf3     = tr$Data$Mut_orf3,
      stringsAsFactors = FALSE
    )
  }

  if (out_n == 0L)
    ok_ret(empty_df())
  else
    ok_ret(do.call(rbind, out_rows))
}



##################################################################################################
## validate_nmID_for_pipeline
##
## PURPOSE
##   Validate that a given RefSeq transcript (nmID), defined in refseq_df,
##   is structurally and numerically compatible with THIS indel / splice /
##   translation pipeline.
##
##   This function enforces ALL transcript-level invariants that are
##   ASSUMED (but not re-checked) throughout downstream analysis.
##
##
## ===================== INPUT =====================
##
##   nmID : character scalar
##     RefSeq transcript ID (e.g. "NM_000546")
##
##   refseq_df : data.frame
##     RefSeq annotation table containing (at minimum) the following columns:
##
##       nmID
##       txS, txE
##       cdsS, cdsE
##       exonS, exonE
##       chr
##       strand
##
##   consts : list
##     Pipeline constants object.
##     Currently not consumed directly, but accepted to:
##       - enforce uniform call signatures
##       - allow future policy-based validation extensions
##
##
## ===================== PIPELINE ASSUMPTIONS ENFORCED =====================
##
##   1. Transcript coordinates
##      - txS < txE
##      - cdsS < cdsE
##      - CDS fully contained within transcript
##
##   2. Exon structure
##      - exonS and exonE have equal length
##      - exon starts < exon ends
##      - exons are non-overlapping
##      - exons are strictly ordered in genomic space
##
##   3. Intron semantics
##      - Every intron has a well-defined downstream exon (Exon_3)
##      - Transcript contains at least one exon
##
##   4. CDS / exon consistency
##      - CDS start and end fall within exon boundaries
##      - CDS length (sum over coding exons) is divisible by 3
##
##   5. Strand validity
##      - strand ∈ { "+", "-" }
##
##
## ===================== OUTPUT =====================
##
##   list(
##     Status  = logical scalar,
##     Message = character scalar
##   )
##
##   Status = TRUE
##     → nmID is SAFE to use in this pipeline
##
##   Status = FALSE
##     → nmID violates one or more invariants
##       Message reports ALL detected failures
##
##
## ===================== FAILURE POLICY =====================
##
##   - NEVER throws
##   - Return shape is invariant
##   - Message is deterministic and human-readable
##
##################################################################################################
validate_nmID_for_pipeline <- function(nmID, refseq, consts) {

  fail <- function(msgs) {
    list(
      Status  = FALSE,
      Message = paste(msgs, collapse = "; ")
    )
  }

  ok <- function() list(Status = TRUE, Message = "")

  ## ---------------- Basic validation ----------------
  if (!is.character(nmID) || length(nmID) != 1L)
    return(fail("nmID must be a character scalar"))

  if (!is.data.frame(refseq))
    return(fail("refseq must be a data.frame"))

  if (!"nmID" %in% names(refseq))
    return(fail("refseq missing column: nmID"))

  idx <- which(refseq$nmID == nmID)
  if (length(idx) == 0L)
    return(fail(paste0("nmID not found in refseq: ", nmID)))

  if (length(idx) > 1L)
    return(fail(paste0("nmID is not unique in refseq: ", nmID)))

  row <- refseq[idx, , drop = FALSE]

  required_cols <- c("txS","txE","cdsS","cdsE","exonS","exonE","chr","strand")
  if (!all(required_cols %in% names(row)))
    return(fail(paste(
      "missing required columns:",
      paste(setdiff(required_cols, names(row)), collapse = ",")
    )))

  msgs <- character(0)

  ## ---------------- Transcript coordinates ----------------
  txS  <- as.integer(row$txS)
  txE  <- as.integer(row$txE)
  cdsS <- as.integer(row$cdsS)
  cdsE <- as.integer(row$cdsE)

  if (is.na(txS) || is.na(txE) || txS >= txE)
    msgs <- c(msgs, "invalid transcript coordinates: txS >= txE")

  if (is.na(cdsS) || is.na(cdsE) || cdsS >= cdsE)
    msgs <- c(msgs, "invalid CDS coordinates: cdsS >= cdsE")

  if (cdsS < txS || cdsE > txE)
    msgs <- c(msgs, "CDS not fully contained within transcript")

  ## ---------------- Strand ----------------
  strand <- as.character(row$strand)
  if (!strand %in% c("+", "-"))
    msgs <- c(msgs, "strand must be '+' or '-'")

  ## ---------------- Exon parsing ----------------
  exonS <- suppressWarnings(as.integer(strsplit(row$exonS, ",", fixed = TRUE)[[1]]))
  exonE <- suppressWarnings(as.integer(strsplit(row$exonE, ",", fixed = TRUE)[[1]]))

  if (length(exonS) == 0L || length(exonE) == 0L)
    msgs <- c(msgs, "empty exon definition")

  if (length(exonS) != length(exonE))
    msgs <- c(msgs, "exonS and exonE length mismatch")

  if (length(exonS) < 1L)
    msgs <- c(msgs, "transcript must contain at least one exon")

  ## ---------------- Exon geometry ----------------
  for (i in seq_along(exonS)) {
    if (is.na(exonS[i]) || is.na(exonE[i]) || exonS[i] >= exonE[i])
      msgs <- c(msgs, paste0("invalid exon coordinates at exon ", i))
  }

  ## Strict ordering / non-overlap
  if (length(exonS) > 1L) {
    ord <- order(exonS)
    exonS <- exonS[ord]
    exonE <- exonE[ord]
    if (any(exonS[-1] < exonE[-length(exonE)]))
      msgs <- c(msgs, "exons overlap or are not strictly ordered")
  }

  ## ---------------- CDS–exon compatibility ----------------
  cds_start_in_exon <- any(cdsS >= exonS & cdsS < exonE)
  cds_end_in_exon   <- any(cdsE > exonS & cdsE <= exonE)

  if (!cds_start_in_exon)
    msgs <- c(msgs, "CDS start does not fall within any exon")

  if (!cds_end_in_exon)
    msgs <- c(msgs, "CDS end does not fall within any exon")

  ## ---------------- CDS length (EXON-AWARE) ----------------
  ## NOTE:
  ## CDS length MUST be computed as the sum of coding bases across exons.
  ## cdsE - cdsS is NOT valid for multi-exon transcripts.
  cds_len <- 0L
  for (i in seq_along(exonS)) {
    cs <- max(exonS[i], cdsS)
    ce <- min(exonE[i], cdsE)
    if (cs < ce)
      cds_len <- cds_len + (ce - cs)
  }

  if (cds_len <= 0L || (cds_len %% 3L) != 0L)
    msgs <- c(msgs, "CDS length (sum over coding exons) is not a positive multiple of 3")

  ## ---------------- Final decision ----------------
  if (length(msgs))
    return(fail(msgs))

  ok()
}

##################################################################################################
## annoverse_indel_one_transcript
##
## PURPOSE
##   Annotate indels for ONE RefSeq transcript (nmID).
##
##   This is a TERMINAL, gene-level orchestration function. It drives
##   the complete indel annotation workflow for a single transcript by
##   coordinating:
##
##     - optional transcript validation
##     - transcript construction
##     - exon / intron region extraction with buffering
##     - splice-impact analysis
##     - pure coding-impact analysis (exons only)
##
##
## ===================== EXECUTION / ERROR POLICY =====================
##
##   - This is a TERMINAL function.
##   - Inner functions are assumed to return:
##       list(Status, Message, Data)
##
##   - If ANY inner call returns Status == FALSE:
##       1) its Message is printed immediately
##       2) execution ABORTS for this transcript
##       3) Status == FALSE is returned with that Message
##
##
## ===================== OUTPUT =====================
##
##   Data columns (ORDER GUARANTEED):
##
##     Marker
##     GeneName
##     RefSeq
##     Region
##     SearchFor
##     Class
##     Affinity_WT
##     Affinity_Mut
##     Splice_Core_WT
##     Splice_Core_Mut
##     AAChange
##     WT_orf
##     Mut_orf5
##     Mut_orf3
##
##################################################################################################
annoverse_indel_one_transcript <- function(the_indels,nmID,refseq,consts,validate_transcript = TRUE) {
  default_msg <- "in annoverse_indel_one_transcript"

  ret <- function(Status = FALSE, Message = default_msg, Data = data.frame())
    list(Status = Status, Message = if (Status) "" else Message, Data = Data)

  ## ------------------------------------------------------------
  ## Locate transcript
  ## ------------------------------------------------------------
  idx <- which(refseq$nmID == nmID)
  if (length(idx) != 1L) {
    msg <- paste0("check existence of ", nmID)
    print(msg)
    return(ret(FALSE, msg, data.frame()))
  }

  GeneName <- as.character(refseq$gName[idx])
  RefSeqID <- nmID

  ## ---------------- FINAL output schema ----------------
  final <- data.frame(
    Marker          = character(),
    GeneName        = character(),
    RefSeq          = character(),
    Region          = character(),
    SearchFor       = character(),
    Class           = character(),
    Affinity_WT     = numeric(),
    Affinity_Mut    = numeric(),
    Splice_Core_WT  = character(),
    Splice_Core_Mut = character(),
    AAChange        = character(),
    WT_orf          = integer(),
    Mut_orf5        = integer(),
    Mut_orf3        = integer(),
    stringsAsFactors = FALSE
  )

  abort_with_message <- function(msg) {
    print(msg)
    ret(FALSE, msg, final)
  }

  ## ------------------------------------------------------------
  ## OPTIONAL transcript validation
  ## ------------------------------------------------------------
  if (isTRUE(validate_transcript)) {
    v <- validate_nmID_for_pipeline(nmID, refseq, consts)
    if (!isTRUE(v$Status)) {
      print(v$Message)
      return(ret(FALSE, v$Message, final))
    }
  }

  ## ------------------------------------------------------------
  ## Transcript metadata
  ## ------------------------------------------------------------
  gene_chr    <- refseq$chr[idx]
  gene_strand <- refseq$strand[idx]
  gene_txS    <- as.integer(refseq$txS[idx]) + 1L
  gene_txE    <- as.integer(refseq$txE[idx])

  exon_starts <- as.integer(strsplit(refseq$exonS[idx], ",", fixed = TRUE)[[1]]) + 1L
  exon_ends   <- as.integer(strsplit(refseq$exonE[idx], ",", fixed = TRUE)[[1]])

  ord <- order(exon_starts)
  exon_starts <- exon_starts[ord]
  exon_ends   <- exon_ends[ord]
  exon_count  <- length(exon_starts)

  if (exon_count < 1L || length(exon_ends) != exon_count)
    return(abort_with_message(paste0("check validity of exon definition for ", nmID)))

  ## ------------------------------------------------------------
  ## Filter indels to transcript span
  ## ------------------------------------------------------------
  indels <- the_indels[
    the_indels$Chr == gene_chr &
    the_indels$Start_G >= gene_txS &
    the_indels$End_G   <= gene_txE,
    , drop = FALSE
  ]
  if (!nrow(indels)) return(ret(TRUE, "", final))

  ## ------------------------------------------------------------
  ## Extract mRNA
  ## ------------------------------------------------------------
  tmpfile_prefix <- paste0("tmp_", nmID)
  mRNA <- extract_mRNA_in_transcript_orientation(nmID, refseq, consts, tmpfile_prefix)
  if (!isTRUE(mRNA$Status)) return(abort_with_message(mRNA$Message))
  mRNA <- mRNA$Data

  ## ------------------------------------------------------------
  ## Iterate exons / introns
  ## ------------------------------------------------------------
  for (i in seq_len(exon_count)) {

    ## ===================== EXON =====================
    exon_start <- exon_starts[i]
    exon_end   <- exon_ends[i]

    exon_indels <- indels[indels$Start_G >= exon_start - (consts$Max_Indel_Length - 1L) & indels$End_G   <= exon_end   + (consts$Max_Indel_Length - 1L), , drop = FALSE]

    if (nrow(exon_indels)) {
	  region_type <- "Exon"
	  this_exon_number <- i
	  indels <- exon_indels

      exon_df <- extract_one_exon_or_intron_region_with_buffer(gene_chr, exon_start, exon_end, this_exon_number, region_type,tmpfile_prefix, gene_strand, mRNA, consts)
      if (!isTRUE(exon_df$Status)) return(abort_with_message(exon_df$Message))
      region_df <- exon_df$Data
      splice_rlt <- check_mutant_mRNA_for_indels_with_splice_consideration(indels, region_df, region_type, gene_chr, gene_strand, mRNA, consts)
	  
      if (!isTRUE(splice_rlt$Status)) return(abort_with_message(splice_rlt$Message))

      coding_rlt <- check_mutant_mRNA_for_indels_fully_inside_exon(indels, region_df, this_exon_number,gene_chr, gene_strand, mRNA, consts)
	  
	  
      if (!isTRUE(coding_rlt$Status)) return(abort_with_message(coding_rlt$Message))

      if (nrow(splice_rlt$Data)) {
        splice_rlt$Data$GeneName <- GeneName
        splice_rlt$Data$RefSeq   <- RefSeqID
        splice_rlt$Data$WT_orf   <- 0L
        splice_rlt$Data$Mut_orf5 <- 0L
        splice_rlt$Data$Mut_orf3 <- 0L
        splice_rlt$Data <- splice_rlt$Data[, names(final), drop = FALSE]
        final <- rbind(final, splice_rlt$Data)
      }

      if (nrow(coding_rlt$Data)) {
        coding_rlt$Data$GeneName        <- GeneName
        coding_rlt$Data$RefSeq          <- RefSeqID
        coding_rlt$Data$Splice_Core_WT  <- " "
        coding_rlt$Data$Splice_Core_Mut <- " "
        coding_rlt$Data <- coding_rlt$Data[, names(final), drop = FALSE]
        final <- rbind(final, coding_rlt$Data)
      }
    }
    ## ===================== INTRON =====================
    if (i < exon_count) {
      intron_start <- exon_ends[i] + 1L
      intron_end   <- exon_starts[i + 1L] - 1L
      intron_indels <- indels[indels$Start_G >= intron_start - (consts$Max_Indel_Length - 1L) &indels$End_G <= intron_end + (consts$Max_Indel_Length - 1L),, drop = FALSE]
      if (nrow(intron_indels)) {
	    region_type <- "Intron"
	    this_exon_number <- i
	    indels <- intron_indels
        intron_df <- extract_one_exon_or_intron_region_with_buffer(gene_chr, intron_start, intron_end, this_exon_number, region_type, tmpfile_prefix, gene_strand, mRNA, consts)
        if (!isTRUE(intron_df$Status)) return(abort_with_message(intron_df$Message))
		region_df<-intron_df$Data
        splice_rlt <- check_mutant_mRNA_for_indels_with_splice_consideration(indels,region_df , region_type,gene_chr, gene_strand, mRNA, consts)
		
        if (!isTRUE(splice_rlt$Status)) return(abort_with_message(splice_rlt$Message))

        if (nrow(splice_rlt$Data)) {
          splice_rlt$Data$GeneName <- GeneName
          splice_rlt$Data$RefSeq   <- RefSeqID
          splice_rlt$Data$WT_orf   <- 0L
          splice_rlt$Data$Mut_orf5 <- 0L
          splice_rlt$Data$Mut_orf3 <- 0L
          splice_rlt$Data <- splice_rlt$Data[, names(final), drop = FALSE]
          final <- rbind(final, splice_rlt$Data)
        }
      }
    }
  }

  final <- final[order(final$Marker), , drop = FALSE]
  ret(Status = TRUE, Data = final)
}


##################################################################################################
## coarsegrained_gene_level_annotation_snvindel
##
## PURPOSE
##   Perform a coarse-grained, gene-level RefSeq annotation for SNVs and indels.
##   Each indel is mapped to RefSeq transcript genomic spans using BOTH its
##   start and end positions, with chromosome matching enforced.
##
##   This function assigns:
##     - RefSeq : one or more RefSeq transcript IDs (nmID)
##     - Region : coarse region label ("Gene", "Enhancer", or "Intergenic")
##
##   It is intentionally conservative and lightweight, serving as an EARLY
##   annotation stage before transcript-level or splice/translation analysis.
##
##
## INPUT
##   indels_all : data.frame
##     Required columns:
##       - Chr     : chromosome name (character)
##       - Start_G : genomic start position (integer / numeric)
##       - End_G   : genomic end position (integer / numeric)
##       - Ref     : reference allele (character)
##       - Mut     : mutant allele (character)
##
##   consts : list
##     Constants object returned by initialize_constants().
##     This function does not directly consume consts fields, but the argument
##     is required for pipeline consistency and future policy enforcement.
##
##   refseq : data.frame
##     RefSeq transcript annotation table with columns:
##       - chr       : chromosome name (character; lowercase)
##       - txS       : transcript genomic start
##       - txE       : transcript genomic end
##       - cdsS      : CDS start
##       - cdsE      : CDS end
##       - exonS     : exon start positions
##       - exonE     : exon end positions
##       - gName     : gene symbol
##       - nmID      : RefSeq transcript ID (e.g., NM_*)
##       - isDefault : character flag, values "Yes" or "No"
##
##   useDefault : logical, default = TRUE
##     Transcript selection policy when multiple RefSeq transcripts overlap
##     the same indel.
##
##     - TRUE:
##         * For each gene (gName), select ONLY transcripts with
##           isDefault == "Yes"
##         * If multiple genes overlap, record ALL their default nmID values
##           joined by ";"
##
##     - FALSE:
##         * Record ALL overlapping nmID values across all genes, joined by ";"
##
##   enhancer_length : integer, default = 0L
##     Length (bp) to extend transcript boundaries on both sides to define a
##     putative enhancer region.
##
##       - Gene region:
##           txS <= Start_G AND txE >= End_G
##
##       - Enhancer region:
##           (txS - enhancer_length) <= Start_G AND
##           (txE + enhancer_length) >= End_G
##         but NOT fully contained within [txS, txE]
##
##
## OUTPUT
##   A named list with three elements:
##
##     $Status : logical scalar
##       - TRUE  → annotation completed successfully
##       - FALSE → validation or execution problem occurred
##
##     $Message : character scalar
##       - Empty string if Status == TRUE
##       - Non-empty diagnostic message if Status == FALSE
##
##     $Data : data.frame or NULL
##       - If Status == TRUE:
##           The input indels_all with TWO additional columns:
##             - RefSeq : RefSeq transcript ID(s) or "Intergenic"
##             - Region : "Gene", "Enhancer", or "Intergenic"
##
##       - If Status == FALSE:
##           NULL
##
##
## MATCHING POLICY
##   - Chromosome must match exactly:
##         refseq$chr == indels_all$Chr
##
##   - Both Start_G AND End_G must satisfy the containment rule.
##
##
## FAILURE POLICY
##   - This function NEVER throws or calls stop().
##   - ALL errors are captured internally.
##   - Any error results in:
##         Status  = FALSE
##         Message = descriptive error string
##         Data    = NULL
##
##################################################################################################
coarsegrained_gene_level_annotation_snvindel <- function(indels_all,consts,refseq,useDefault = TRUE,enhancer_length = 0L) {

  ## ---------------- Initialize return object ----------------
  ret <- list(
    Status  = TRUE,
    Message = "",
    Data    = NULL
  )

  ## ---------------- Error-safe execution ----------------
  tryCatch({

    ## ---- Validate indels_all ----
    if (!is.data.frame(indels_all)) {
      ret$Status  <- FALSE
      ret$Message <- "indels_all must be a data.frame"
      return(ret)
    }

    required_indel_cols <- c("Chr", "Start_G", "End_G", "Ref", "Mut")
    missing_indel_cols <- setdiff(required_indel_cols, colnames(indels_all))
    if (length(missing_indel_cols) > 0L) {
      ret$Status  <- FALSE
      ret$Message <- paste(
        "indels_all missing required columns:",
        paste(missing_indel_cols, collapse = ", ")
      )
      return(ret)
    }

    if (!is.numeric(indels_all$Start_G) || !is.numeric(indels_all$End_G)) {
      ret$Status  <- FALSE
      ret$Message <- "Start_G and End_G must be numeric"
      return(ret)
    }

    ## ---- Validate consts (lightweight) ----
    if (!is.list(consts)) {
      ret$Status  <- FALSE
      ret$Message <- "consts must be a list returned by initialize_constants()"
      return(ret)
    }

    ## ---- Validate refseq ----
    if (!is.data.frame(refseq)) {
      ret$Status  <- FALSE
      ret$Message <- "refseq must be a data.frame"
      return(ret)
    }

    required_refseq_cols <- c(
      "chr", "txS", "txE", "cdsS", "cdsE",
      "exonS", "exonE", "gName", "nmID", "isDefault"
    )
    missing_refseq_cols <- setdiff(required_refseq_cols, colnames(refseq))
    if (length(missing_refseq_cols) > 0L) {
      ret$Status  <- FALSE
      ret$Message <- paste(
        "refseq missing required columns:",
        paste(missing_refseq_cols, collapse = ", ")
      )
      return(ret)
    }

    if (!all(refseq$isDefault %in% c("Yes", "No"))) {
      ret$Status  <- FALSE
      ret$Message <- "refseq$isDefault must contain only 'Yes' or 'No'"
      return(ret)
    }

    if (!is.logical(useDefault) || length(useDefault) != 1L) {
      ret$Status  <- FALSE
      ret$Message <- "useDefault must be a logical scalar"
      return(ret)
    }

    enhancer_length <- as.integer(enhancer_length)
    if (is.na(enhancer_length) || enhancer_length < 0L) {
      ret$Status  <- FALSE
      ret$Message <- "enhancer_length must be a non-negative integer"
      return(ret)
    }

    ## ---------------- Core annotation logic ----------------
    out <- indels_all
    out$RefSeq <- NA_character_
    out$Region <- NA_character_

    for (i in seq_len(nrow(out))) {

      chr_i <- out$Chr[i]
      s_pos <- out$Start_G[i]
      e_pos <- out$End_G[i]

      ## -------- Gene region hits --------
      gene_hits <- refseq[
        refseq$chr == chr_i &
          refseq$txS <= s_pos &
          refseq$txE >= e_pos,
        ,
        drop = FALSE
      ]

      if (nrow(gene_hits) > 0L) {

        if (useDefault) {
          gene_by_gene <- split(gene_hits, gene_hits$gName)
          nmids <- character(0)

          for (g in names(gene_by_gene)) {
            g_df <- gene_by_gene[[g]]
            def_df <- g_df[g_df$isDefault == "Yes", , drop = FALSE]
            if (nrow(def_df) > 0L)
              nmids <- c(nmids, def_df$nmID)
          }

          if (length(nmids) > 0L) {
            out$RefSeq[i] <- paste(unique(nmids), collapse = ";")
            out$Region[i] <- "Gene"
            next
          }

        } else {
          out$RefSeq[i] <- paste(unique(gene_hits$nmID), collapse = ";")
          out$Region[i] <- "Gene"
          next
        }
      }

      ## -------- Enhancer region hits --------
      if (enhancer_length > 0L) {

        enh_hits <- refseq[
          refseq$chr == chr_i &
            (refseq$txS - enhancer_length) <= s_pos &
            (refseq$txE + enhancer_length) >= e_pos,
          ,
          drop = FALSE
        ]

        ## remove true gene hits
        enh_hits <- enh_hits[
          !(enh_hits$txS <= s_pos & enh_hits$txE >= e_pos),
          ,
          drop = FALSE
        ]

        if (nrow(enh_hits) > 0L) {

          if (useDefault) {
            enh_by_gene <- split(enh_hits, enh_hits$gName)
            nmids <- character(0)

            for (g in names(enh_by_gene)) {
              g_df <- enh_by_gene[[g]]
              def_df <- g_df[g_df$isDefault == "Yes", , drop = FALSE]
              if (nrow(def_df) > 0L)
                nmids <- c(nmids, def_df$nmID)
            }

            if (length(nmids) > 0L) {
              out$RefSeq[i] <- paste(unique(nmids), collapse = ";")
              out$Region[i] <- "Enhancer"
              next
            }

          } else {
            out$RefSeq[i] <- paste(unique(enh_hits$nmID), collapse = ";")
            out$Region[i] <- "Enhancer"
            next
          }
        }
      }

      ## -------- Intergenic --------
      out$RefSeq[i] <- "Intergenic"
      out$Region[i] <- "Intergenic"
    }

    ret$Data <- out

  }, error = function(e) {
    ret$Status  <- FALSE
    ret$Message <- conditionMessage(e)
    ret$Data    <- NULL
  })

  ret
}

##################################################################################################
## annoverse_snvindel_annotation
##
## PURPOSE
##   Perform end-to-end annotation of SNVs and indels using the ANNO-VERSE
##   pipeline, while GUARANTEEING that every input indel appears in the
##   final output (including those not analyzed at transcript level).
##
##   This function is the TOP-LEVEL orchestration layer and coordinates:
##
##     1) Reading input indels from file
##     2) Initializing build-specific constants
##     3) Coarse-grained gene / enhancer / intergenic annotation
##     4) Transcript-level annotation for GENE-region indels only
##     5) Placeholder annotation for ENHANCER-region indels
##     6) Explicit recovery of NOT-ANALYZED indels via table lookup
##     7) Merging all results and writing them to an output file
##
##   An indel may legitimately appear in MULTIPLE output rows if it maps to
##   multiple RefSeq transcripts; this is expected behavior.
##
##
## INPUT
##   infile : character scalar
##     Path to a tab-delimited file containing indels to annotate.
##     Required columns:
##       - Marker   : unique identifier for each indel (REQUIRED for lookup)
##       - Chr
##       - Start_G
##       - End_G
##       - Ref
##       - Mut
##
##   outfile : character scalar
##     Path to the output file (tab-delimited). Existing files are overwritten.
##
##   build : character scalar
##     Genome build identifier.
##     Allowed values:
##       - "hg19"
##       - "hg38"
##
##
## EXECUTION / ERROR POLICY
##   - This function NEVER throws or calls stop().
##   - All failures are reported via:
##       list(Status, Message, Data)
##
##
## OUTPUT
##   A named list with three elements:
##
##     $Status : logical scalar
##       - TRUE  → annotation completed successfully
##       - FALSE → fatal error occurred
##
##     $Message : character scalar
##       - Empty string if Status == TRUE
##       - Non-empty diagnostic message if Status == FALSE
##
##     $Data : data.frame or NULL
##       - Final merged annotation table written to outfile
##
##
## DEPENDENCIES
##   The following functions must already exist and follow the
##   list(Status, Message, Data) contract:
##
##     - initialize_constants()
##     - coarsegrained_gene_level_annotation_snvindel()
##     - annoverse_indel_one_transcript()
##
##################################################################################################
annoverse_snvindel_annotation <- function(infile,
                                          outfile,
                                          build) {

  ## ---------------- Initialize return object ----------------
  ret <- list(Status = TRUE, Message = "", Data = NULL)

  ## ---------------- Error-safe execution ----------------
  tryCatch({

    ## ---- Validate inputs ----
    if (!is.character(infile) || length(infile) != 1L || !file.exists(infile)) {
      ret$Status  <- FALSE
      ret$Message <- "infile must be an existing file"
      return(ret)
    }

    if (!is.character(outfile) || length(outfile) != 1L) {
      ret$Status  <- FALSE
      ret$Message <- "outfile must be a character scalar"
      return(ret)
    }

    if (!build %in% c("hg19", "hg38")) {
      ret$Status  <- FALSE
      ret$Message <- "build must be 'hg19' or 'hg38'"
      return(ret)
    }

    ## ---- Load indels ----
    indels_all <- tryCatch(
      read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(indels_all)) {
      ret$Status  <- FALSE
      ret$Message <- "failed to read indel input file"
      return(ret)
    }

    required_cols <- c("Marker", "Chr", "Start_G", "End_G", "Ref", "Mut")
    missing_cols <- setdiff(required_cols, colnames(indels_all))
    if (length(missing_cols) > 0L) {
      ret$Status  <- FALSE
      ret$Message <- paste("input indel file missing columns:",
                           paste(missing_cols, collapse = ", "))
      return(ret)
    }

    ## ---- Initialize constants ----
    consts <- initialize_constants(build = build)

    ## ---- Load RefSeq annotation ----
    refseq <- tryCatch(
      read.table(consts$RefSeq, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(refseq)) {
      ret$Status  <- FALSE
      ret$Message <- "failed to load RefSeq annotation table"
      return(ret)
    }

    ## ---- Coarse-grained gene-level annotation ----
    cg <- coarsegrained_gene_level_annotation_snvindel(
      indels_all = indels_all,
      consts     = consts,
      refseq     = refseq,
      useDefault = TRUE
    )
    if (!isTRUE(cg$Status)) {
      ret$Status  <- FALSE
      ret$Message <- cg$Message
      return(ret)
    }
    cg_df <- cg$Data

    ## ---------------- Transcript-level annotation (GENE only) ----------------
    gene_rows <- cg_df[cg_df$Region == "Gene" & !is.na(cg_df$RefSeq), , drop = FALSE]

    refseq_ids <- unique(unlist(strsplit(
      gene_rows$RefSeq,
      ";",
      fixed = TRUE
    )))

    tx_results <- list()

    for (nmID in refseq_ids) {

      rlt <- annoverse_indel_one_transcript(
        the_indels = indels_all,
        nmID       = nmID,
        refseq     = refseq,
        consts     = consts,
        validate_transcript = TRUE
      )

      if (!isTRUE(rlt$Status)) {
        ret$Status  <- FALSE
        ret$Message <- rlt$Message
        return(ret)
      }

      if (nrow(rlt$Data)) {
        tx_results[[nmID]] <- rlt$Data
      }
    }

    tx_df <- if (length(tx_results)) do.call(rbind, tx_results) else data.frame()

    ## ---------------- ENHANCER placeholders ----------------
    enhancer_rows <- cg_df[cg_df$Region == "Enhancer", , drop = FALSE]

    enhancer_df <- data.frame()
    if (nrow(enhancer_rows)) {
      enhancer_df <- data.frame(
        Marker          = enhancer_rows$Marker,
        GeneName        = enhancer_rows$RefSeq,
        RefSeq          = enhancer_rows$RefSeq,
        Region          = "Enhancer",
        SearchFor       = "Enhancer",
        Class           = "Enhancer",
        Affinity_WT     = -9,
        Affinity_Mut    = -9,
        Splice_Core_WT  = enhancer_rows$Ref,
        Splice_Core_Mut = enhancer_rows$Mut,
        AAChange        = "notanalyzed",
        WT_orf          = -9L,
        Mut_orf5        = -9L,
        Mut_orf3        = -9L,
        stringsAsFactors = FALSE
      )
    }

    ## ---------------- Combine analyzed results ----------------
    combined_df <- rbind(tx_df, enhancer_df)

    ## ---------------- Recover NOT-ANALYZED indels ----------------
    analyzed_markers <- if (nrow(combined_df))
      unique(combined_df$Marker) else character(0)

    missing_rows <- indels_all[!(indels_all$Marker %in% analyzed_markers), , drop = FALSE]

    notanalyzed_df <- data.frame()
    if (nrow(missing_rows)) {
      notanalyzed_df <- data.frame(
        Marker          = missing_rows$Marker,
        GeneName        = "notanalyzed",
        RefSeq          = "notanalyzed",
        Region          = "notanalyzed",
        SearchFor       = "notanalyzed",
        Class           = "notanalyzed",
        Affinity_WT     = -9,
        Affinity_Mut    = -9,
        Splice_Core_WT  = missing_rows$Ref,
        Splice_Core_Mut = missing_rows$Mut,
        AAChange        = "notanalyzed",
        WT_orf          = -9L,
        Mut_orf5        = -9L,
        Mut_orf3        = -9L,
        stringsAsFactors = FALSE
      )
    }

    ## ---------------- Final merge and write output ----------------
    final_df <- rbind(combined_df, notanalyzed_df)

    write.table(
      final_df,
      file = outfile,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )

    ret$Data <- final_df

  }, error = function(e) {
    ret$Status  <- FALSE
    ret$Message <- conditionMessage(e)
    ret$Data    <- NULL
  })

  ret
}



annoverse_snvindel_annotation("for_testing.txt","result.txt","hg38")


