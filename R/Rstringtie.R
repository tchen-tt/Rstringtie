#' @title R wrapper to Run stringtie tool
#' @description R weapper to run Stringtie enables reconstruction of a
#'  transcriptome from RNA-seq reads
#'
#' @param bam A bam file, must be a SAM, BAM or CRAM file with RNA-Seq read
#'  alignments sorted by their genomic location.
#' @param reference Use a reference annotation file to guide assembly process
#' @param outfile Output GTF file
#' @param params Other parameter
#'
#' @export
stringtieAssembly <- function(bam, reference, outfile, params = "") {
  reference <- paste("-G", reference, sep = " ")
  outfile <- paste0("-o", outfile, sep = " ")
  cmd <- sprintf("%s %s %s %s",
                 bam,
                 reference,
                 outfile,
                 params)
  return(invisible(lapply(cmd, .Stringtiebin)))
}

#' @title R wrapper to Run stringtie tool
#' @description Transcript merge mode
#'
#' @param reference Use a reference annotation file to guide assembly process
#' @param gtffile Gtf files with stringtie assemble transcripts [character vector]
#' @param outfile Outfile for the merged transcripts GTF
#' @param params Other parameter
#' @export
stringtieMerge <- function(reference, gtffile, outfile, params = "") {
  programtags <- "--merge"
  reference <- paste("-G", reference, sep = " ")
  gtffile <- paste(gtffile, collapse = " ")
  outfile <- paste("-o", outfile, sep = " ")
  cmd <- sprintf("%s %s %s %s %s",
                 programtags,
                 reference,
                 gtffile,
                 outfile,
                 params)
  return(invisible(lapply(cmd, .Stringtiebin)))
}

#' @title R wrapper to Run gffcompare tool
#' @description The function can be used to compare, merge, annotate and estimate
#' accuracy of one or more GFF files (the “query” files), when compared with a
#' reference annotation (also provided as GFF).
#' @param reference Use a reference annotation file to guide compare assembly process.
#' @param gtffile Gtf files with gffcompare annotation transcripts.
#' @param outfile Outfile for the annotation transcripts GTF.
#' @param params Other parameter
#' @export

gffcompareAnno <- function(reference, gtffile, outfile, params = "") {
  reference <- paste("-r", reference, sep = " ")
  outfile <- gsub(pattern = "[.]gtf$", replacement = "", x = outfile)
  gtffile <- paste(gtffile, collapse = " ")
  outfile <- paste("-o", outfile, sep = " ")
  cmd <- sprintf("%s %s %s %s",
                 reference,
                 outfile,
                 gtffile,
                 params)
  return(invisible(lapply(cmd, .gffcompareBin)))
}

#' @title Preparing the genome annotation object
#' @description GTF file, many exons appear multiple times, once for each transcripts
#' that contains them, need to 'collapse' this information to define exon counting bins.
#' @param gtffile GTF file.
#' @param singleGens Methods of the disjoint parts, default is TRUE
#' @param transposne A GRanges object with transposones data
#' @param minoverlap Minimum overlap wit findoverlap, default is 10.
#' @importFrom rtracklayer import.gff
#' @importFrom GenomicFeatures exonicParts
#' @importFrom GenomicFeatures makeTxDbFromGRanges
#' @importFrom GenomicRanges strand
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors mcols
#' @export
prepareAnno <- function(gtffile, singleGens = TRUE, transposone = NULL, minoverlap = 10) {
  gtfGr <- rtracklayer::import.gff(con = gtffile)
  message("Remove transcripts which have missing strand information.")
  gtfGr <- gtfGr[!GenomicRanges::strand(gtfGr) == "*"]
  txdb <- GenomicFeatures::makeTxDbFromGRanges(gr = gtfGr)
  exonicParts <- GenomicFeatures::exonicParts(txdb = txdb, 
                                              linked.to.single.gene.only = singleGens)
  exonrank <- split(x = exonicParts$exonic_part, 
                    f = exonicParts$gene_id, drop = TRUE)
  gestrand <- split(x = GenomicRanges::strand(exonicParts), 
                    f = exonicParts$gene_id, drop = TRUE)
  exonicpart <- base::lapply(names(exonrank), FUN = function(x) {
    if(unique(as.character(gestrand[[x]])) == "-") {
      exonrank[[x]] <- rev(exonrank[[x]])
    } else {
      exonrank[[x]] <- exonrank[[x]]
    }
  })
  names(exonicpart) <- names(exonrank)
  exonicpart <- exonicpart[unique(exonicParts$gene_id)]
  exonicParts$exonic_part <- unlist(exonicpart)
  
  if (!is.null(transposone)) {
    if (!c("names") %in% colnames(S4Vectors::mcols(transposone)) || !is(transposone, "GRanges")) {
      stop("Transposone must be a Granges object which cotain names column.")
    }
    overlaps <- IRanges::findOverlaps(query = exonicParts, subject = transposone, 
                                      minoverlap = minoverlap)
    repeats <- split(x = transposone$names[S4Vectors::subjectHits(overlaps)], 
                     f = S4Vectors::queryHits(overlaps))
    repeats <- lapply(X = repeats, FUN = function(x) paste(x, collapse = ","))
    exonicParts$transposone <- "none"
    exonicParts$transposone[as.numeric(names(repeats))] <- unlist(repeats)
  }
  
  return(exonicParts)
}

#' @title Counting reads
#' @description Counting the number of reads that overlap with each of the exon
#' counting bings defined in the flatted GFF file.
#' @param annotation Flatted GFF file.
#' @param bamfile BAM files path
#' @param strandSpecific an integer vector indicating if strand-specific read 
#' counting should be performed. See also \link[Rsubread]{featureCounts}
#' @param allowMultiOverlap logical indicating if a read is allowed to be 
#' assigned to more than one feature (or meta-feature) if it is found to 
#' overlap with more than one feature (or meta-feature). TRUE by default.
#' @param isPairedEnd A logical scalar or a logical vector, indicating whether 
#' libraries contain paired-end reads or not.
#' @param ... Additional arguments. See in \link[Rsubread]{featureCounts}
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom Rsubread featureCounts
#' @export
countAnno <- function(annotation, bamfile, 
                      isPairedEnd = TRUE, 
                      strandSpecific = 0, 
                      allowMultiOverlap = TRUE,
                      ...) {
  #bam <- BamFileList(bamfile)
  #count <- GenomicAlignments::summarizeOverlaps(
  #          features = annotation,
  #          reads = bam,
  #          singleEnd = singleEnd,
  #          ...)
  #return(count)
  annframe <- as.data.frame(annotation)
  names(annframe)[names(annframe) == "seqnames"] <- "Chr"
  annframe$GeneId <- paste(annotation$gene_id, annotation$exonic_part, sep = ":")
  hints <- Rsubread::featureCounts(files = bamfile,
                                   annot.ext = annframe,
                                   strandSpecific = strandSpecific,
                                   allowMultiOverlap = allowMultiOverlap,
                                   isPairedEnd = isPairedEnd,
                                   useMetaFeatures = FALSE,
                                   isGTFAnnotationFile = FALSE,
                                   ...)
  se <- SummarizedExperiment::SummarizedExperiment(list(counts = as.matrix(hints$counts)), 
                                                   rowRanges = annotation)
  return(se)
}


#' Differential expression with flatten exons
#' @description Diffenential expression with exons
#' @param object RangedSummarizedExperiment object
#' @param design Expreiment design
#' @param ...  Ohter arguments
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData
#' @export
diffExons <- function(object, design = vector(), ...) {
  stopifnot(nrow(colData(object)) == length(design))
  colData(object)$condition <- design
  if (is(object, "RangedSummarizedExperiment")) {
    dds <- DEXSeq::DEXSeqDataSetFromSE(SE = object)
    dds <- DEXSeq::estimateSizeFactors(dds)
    dds <- DEXSeq::estimateDispersions(dds, ...)
    dds <- DEXSeq::testForDEU(object = dds, ...)
    dds <- DEXSeq::estimateExonFoldChanges(object = dds, ...)
    res <- DEXSeq::DEXSeqResults(dds)
    
    #reorder
    object <- sort(object)
    res <- res[order(res$genomicData), ]
    
    rowData(object)[names(res)[grep("log2fold",names(res))]] <- res[[grep("log2fold",names(res))]]
    rowData(object)$exonBaseMean <- res$exonBaseMean
    rowData(object)$pvalue <- res$pvalue
    return(object)
  }
}




