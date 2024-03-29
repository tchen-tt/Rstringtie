% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rstringtie.R
\name{stringtieAssembly}
\alias{stringtieAssembly}
\title{R wrapper to Run stringtie tool}
\usage{
stringtieAssembly(bam, reference, outfile, params = "")
}
\arguments{
\item{bam}{A bam file, must be a SAM, BAM or CRAM file with RNA-Seq read
alignments sorted by their genomic location.}

\item{reference}{Use a reference annotation file to guide assembly process}

\item{outfile}{Output GTF file}

\item{params}{Other parameter}
}
\description{
R weapper to run Stringtie enables reconstruction of a
 transcriptome from RNA-seq reads
}
