#' @title Perform short sequence alignment with BLAST
#' @description Perform short sequence alignment with BLAST
#' 
#' @param sequences Character vector of DNA sequences.
#' @param blast_index String specifying path to a BLAST index.
#' @param oppositeStrandOnly If TRUE, only hits for which the target is 
#'     complementary to the query sequences will be returned. 
#' @param codingOnly If TRUE, only hits for which the target is 
#'     located in a coding gene are returned. 
#' @param word_size Integer specifying the BLAST word size parameter.
#'     11 by default. 
#' @param style Either "gencode" (default) or "ensembl". 
#' @param blast_bin String specifying the path of the BLAST binary program. 
#' 
#' @return A data.frame of the alignments. 
#' 
#' @author Jean-Philippe Fortin 
#' 
#' @export 
runBlast <- function(sequences,
                     blast_index,
                     oppositeStrandOnly=TRUE,
                     codingOnly=TRUE,
                     word_size=11,
                     style=c("gencode", "ensembl"),
                     blastn_bin=NULL
){
    style <- match.arg(style)
    results <- lapply(sequences,
                      .runCoreBlast,
                      blastn_bin=blastn_bin,
                      blast_index=blast_index,
                      word_size=word_size)
    results <- lapply(results,
                      .parseBlastResults, 
                      style=style,
                      oppositeStrandOnly=oppositeStrandOnly,
                      codingOnly=codingOnly)
    names(results) <- sequences
    return(results)
}



.runCoreBlast <- function(sequences,
                          blast_index,
                          word_size=11,
                          blastn_bin=NULL
){

    .writeFasta <- function(sequences, file){
          col1 <- paste0(">sirna_", 1:length(sequences))
          col2 <- sequences
          dat <- cbind(as.vector(matrix(c(col1,col2), ncol=2, byrow=TRUE)))
          write.table(dat, quote=FALSE, col.names=FALSE, row.names=FALSE, file=file)
    }

    
    
    inputFile  <- paste0(tempfile(), ".fa")
    outputFile <- paste0(tempfile(), ".txt")
    .writeFasta(sequences, file=inputFile)
    #root <- "blastn -task blastn-short -query "
    #root <- "blastn -query "
    root <- paste0(blastn_bin, " -query ")
    cmd <- paste(root, inputFile,
                 " -db ", blast_index, 
                 " -word_size ",word_size,
                 " -evalue 100 > ",
                 outputFile)
    system(cmd)
    results <- readLines(outputFile)
    return(results)
}



#' @importFrom stringr str_extract
.parseBlastResults <- function(results,
                               style=c("gencode", "ensembl"),
                               oppositeStrandOnly=TRUE,
                               codingOnly=TRUE
){
    style <- match.arg(style)
    # wh <- which(grepl(">", results))
    aln <- which(grepl('^>', results))
    # append ending line of final match result
    final_result <- results[tail(aln, 1):length(results)]
    aln_end <- which(grepl('^Sbjct', final_result))
    aln_end <- tail(aln_end, 1)
    aln <- c(aln, tail(aln, 1)+aln_end)
    
    isEmpty <- sum(grepl("No hits found", results))>0


    

    .processBlock_gencode <- function(block){
        desc <-  gsub("^>", "", block[1])
        desc <- strsplit(desc, split="\\|")
        tx.id <- unlist(lapply(desc, `[[`, 1))
        tx.id <- str_extract(tx.id, "ENS[A-Z]*T[0-9]+")
        tx.name <- unlist(lapply(desc, `[[`, 5))
        gene.id <- unlist(lapply(desc, `[[`, 2))
        gene.id <- str_extract(gene.id, "ENS[A-Z]*G[0-9]+")
        gene.symbol <- unlist(lapply(desc, `[[`, 6))
        biotype <- unlist(lapply(desc, `[[`, 8))
        wh <- which(grepl("Strand=", block))
        strand <- ifelse(block[wh]==" Strand=Plus/Minus", "+-", "++")
        wh <- which(grepl("Score =", block))
        score <- str_extract(block[wh], "Score = [0-9]+\\.{0,1}[0-9]*")
        score <- gsub("Score = ", "", score)
        score <- as.numeric(score)
        evalue <- str_extract(block[wh], "Expect = [0-9]+[e\\.]?-?[0-9]*")
        evalue <- gsub("Expect = ", "", evalue)
        evalue <- as.numeric(evalue)
        wh <- which(grepl("Identities =", block))
        identities <- str_extract(block[wh], "Identities = [0-9]+\\/[0-9]+")
        identities <- gsub("Identities = ", "", identities)
        identities <- strsplit(identities, split='/')
        n.matches <- lapply(identities, `[[`, 1)
        n.matches <- as.numeric(n.matches)
        n.total <- lapply(identities, `[[`, 2)
        n.total <- as.numeric(n.total)
        perc.identity <- round(n.matches/n.total*100)
        wh <- which(grepl("Query", block))
        queryLine <- block[wh]
        matchLine <- block[wh+1]
        subjectLine <- block[wh+2]
        # bring it all together
        data.frame(tx.id = tx.id,
                   tx.name = tx.name,
                   gene.id = gene.id,
                   gene.symbol = gene.symbol,
                   biotype = biotype,
                   strand = strand,
                   score = score,
                   evalue = evalue,
                   n.matches = n.matches,
                   n.total = n.total,
                   perc.identity = perc.identity,
                   queryLine = queryLine,
                   matchLine = matchLine,
                   subjectLine = subjectLine)
    }

    .processBlock_ensembl <- function(block){
      desc <- gsub("^>", "", block[1])
        tx.id <- str_extract(desc, "^ENS[A-Z]*T[0-9]+")
        tx.name <- tx.id
        wh <- which(grepl("gene:", block))
        gene.id <- str_extract(block[wh], "ENS[A-Z]*G[0-9]+")
        wh <- which(grepl("gene_symbol:", block))
        if (length(wh)==0){
          gene.symbol <- NA
        } else {
        gene.symbol <- str_extract(block[wh], "gene_symbol:[^\\s]+")
        gene.symbol <- gsub("gene_symbol:", "", gene.symbol)
        }
        wh <- which(grepl("transcript_biotype:", block))
        biotype <- str_extract(block[wh], "transcript_biotype:[^\\s]+")
        biotype <- gsub("transcript_biotype:", "", biotype)
        wh <- which(grepl("Strand=", block))
        strand <- ifelse(block[wh]==" Strand=Plus/Minus", "+-", "++")
        wh <- which(grepl("Score =", block))
        score <- str_extract(block[wh], "Score = [0-9]+\\.{0,1}[0-9]*")
        score <- gsub("Score = ", "", score)
        score <- as.numeric(score)
        evalue <- str_extract(block[wh], "Expect = [0-9]+[e\\.]?-?[0-9]*")
        evalue <- gsub("Expect = ", "", evalue)
        evalue <- as.numeric(evalue)
        wh <- which(grepl("Identities = ", block))
        identities <- str_extract(block[wh], "Identities = [0-9]+\\/[0-9]+")
        identities <- gsub("Identities = ", "", identities)
        identities <- strsplit(identities, split="/")
        n.matches <- lapply(identities, `[[`, 1)
        n.matches <- as.numeric(unlist(n.matches))
        n.total <- lapply(identities, `[[`, 2)
        n.total <- as.numeric(unlist(n.total))
        perc.identity <- round(n.matches/n.total*100)
        wh <- which(grepl("Query", block))
        queryLine <- block[wh]
        matchLine <- block[wh+1]
        subjectLine <- block[wh+2]
        # bring it all together
        data.frame(tx.id = tx.id,
                   tx.name = tx.name,
                   gene.id = gene.id,
                   gene.symbol = gene.symbol,
                   biotype = biotype,
                   strand = strand,
                   score = score,
                   evalue = evalue,
                   n.matches = n.matches,
                   n.total = n.total,
                   perc.identity = perc.identity,
                   queryLine = queryLine,
                   matchLine = matchLine,
                   subjectLine = subjectLine)
    }

    if (isEmpty){
        blocks <- NULL
        return(blocks)
    }
    if (length(aln)>0){
      indices <- lapply(seq_len(length(aln)-1), function(i){
        aln[i]:(aln[i+1]-1)
      })
        blocks <- lapply(indices, function(i) results[i])
        if (style=="gencode"){
            blocks <- lapply(blocks, .processBlock_gencode) 
        } else {
            blocks <- lapply(blocks, .processBlock_ensembl) 
        }
        blocks <- do.call(rbind,blocks)
    } else {
        blocks <- NULL
        return(blocks)
    }
    if (nrow(blocks)>=1 & codingOnly){
        blocks <- blocks[blocks$biotype=="protein_coding",,drop=FALSE]
    }
    if (nrow(blocks)>=1 & oppositeStrandOnly){
        blocks <- blocks[blocks$strand=="+-",,drop=FALSE]
    }
    return(blocks)
}

