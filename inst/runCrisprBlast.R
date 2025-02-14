# #' @title Run CRISPR Blast
# #' @description Run CRISPR Blast.
# #' 
# #' @param spacers Character vector of DNA sequences. Must be equal in length.
# #' @param nuclease Cas9 or Cas12a (Cpf1)
# #' @param canonical Should only canonical PAM sequences be considered? TRUE by default.
# #' @param species Either human, mouse, rat or fruitfly.
# #' @param what Bowtie index to be used. Either "dna" (whole genome), "rna" (RNA transcripts only) or "exon" (annotated exons only).
# #' @return \strong{runBot} returns spacer alignment data, including genomic coordinates and sequence, and position of mismatches relative to \code{pam.site}.
# #' @name runCrisprBlast
# NULL
# #' @export
# #' @importFrom dplyr rename
# runCrisprBlast <- function(spacers, 
#                         nuclease=c("Cas9", "Cas12a"),
#                         canonical=TRUE,
#                         species = c("human","mouse", "rat", "fruitfly"), 
#                         what=c("dna", "rna", "exon")
# ){  
#     nuclease <- match.arg(nuclease)
#     species  <- match.arg(species)
#     what <- match.arg(what)
#     possiblePams <- .getPossiblePams(nuclease=nuclease, canonical=canonical)
#     if (nuclease=="Cas9"){
#         spacers.with.pam <- lapply(possiblePams, function(x) paste0(spacers, x))
#     }
#     if (nuclease=="Cas12a"){
#         spacers.with.pam <- lapply(possiblePams, function(x) paste0(x, spacers))
#     }
#     spacers.with.pam <- do.call(c, spacers.with.pam)
#     pam.len    <- .getPamLength(nuclease=nuclease)
#     spacer.len <- unique(nchar(spacers))
#     if (length(spacer.len)>1){
#       stop("Spacers must all have the same length.")
#     }
#     aln <- runBlast(sequences=spacers.with.pam, species=species)
#     # if (nrow(aln)>0){
#     #     pam.indices <- .extractPamIndices(aln$target, nuclease=nuclease)
#     #     bad <- which(aln$mm1 %in% pam.indices | aln$mm2 %in% pam.indices | aln$mm3 %in% pam.indices)
#     #     if (length(bad)>0){
#     #         aln <- aln[-bad,,drop=FALSE]
#     #     }
#     # } else {
#     #     return(NULL)
#     # }
#     # if (nrow(aln)>0){
#     #     aln$pam    <- .extractPamFromProtospacer(aln$target, nuclease=nuclease)
#     #     aln$target <- .extractSpacerFromProtospacer(aln$target, nuclease=nuclease)
#     #     aln$query  <- .extractSpacerFromProtospacer(aln$query, nuclease=nuclease)
#     #     aln$canonical <- aln$pam %in% .getPossiblePams(nuclease=nuclease, canonical=TRUE)
#     #     aln <- aln[order(aln$query, aln$n.mismatches),,drop=FALSE]
#     # } else {
#     #     return(NULL)
#     # }
#     # aln <- dplyr::rename(aln, spacer="query")
#     # aln <- dplyr::rename(aln, protospacer="target")
#     # aln$pam.site <- .getBowtiePamSite(pos=aln$pos, 
#     #     strand=aln$strand, 
#     #     spacer.len=spacer.len,
#     #     nuclease=nuclease
#     # )
#     # aln$pos <- NULL
#     # #For Cas12a, need to add PAM len to mismatc position
#     # if (nuclease=="Cas12a"){
#     #     aln[, c("mm1", "mm2", "mm3")] <- aln[, c("mm1", "mm2", "mm3")]-pam.len
#     # }
#     # cols <- c("spacer", "protospacer", "pam",
#     #     "chr", "pam.site", "strand", 
#     #     "n.mismatches", "mm1", "mm2", "mm3", 
#     #     "mmnuc1", "mmnuc2", "mmnuc3", 
#     #     "canonical"
#     # )
#     # aln <- aln[,cols,drop=FALSE]
#     return(aln)
# }





# testCrisprBlast <- function(){
#   spacers <- c("GGAAGTCTGGAGTCTCCAGG", "GTGGACTTCCAGCTACGGCG", "GTGTGCCGAGGTGTGCTGCG", 
#         "TCAATGGTCACAGTAGCGCT", "GCAGATAAAACCAAAAACCG", "GGAATCAAATCACAAATCCG", 
#         "GGATGCCCCCAACAAATAAT", "GGTGCAGCATCCCAACCAGG", "AAGCCTGTGAATCTTCACCG"
#   )
#   runCrisprBlast(spacers, v=2)
# }


# testCrisprBlast2 <- function(){
#   spacers <- c("GGAAGTCTGGAGTCTCCAGG", "GTGGACTTCCAGCTACGGCG", "GTGTGCCGAGGTGTGCTGCG", 
#         "TCAATGGTCACAGTAGCGCT", "GCAGATAAAACCAAAAACCG", "GGAATCAAATCACAAATCCG", 
#         "GGATGCCCCCAACAAATAAT", "GGTGCAGCATCCCAACCAGG", "AAGCCTGTGAATCTTCACCG"
#   )
#   runCrisprBlast(spacers, v=2, nuclease="Cas12a")
# }





