library(argparse)
library(MutationTimeR)

readCnTable <- function(cn_path, clonal_freq) {
    cn = read.table(cn_path, header=1)
    bb <- GRanges( # bb file formated for MutationTimeR
      seqnames = Rle( as.character(cn$chromosome) ),
      ranges = IRanges(
        as.numeric(cn$start),
        end = as.numeric(cn$end)
      ),
      major_cn = as.integer(cn$major_1),
      minor_cn = as.integer(cn$minor_1),
      clonal_frequency = rep(clonal_freq, dim(cn)[1]), 
    )
    return(bb)
}

get_args <- function() {
    p <- ArgumentParser(description = "Run MutationTimeR and plot results")

    p$add_argument("vcf", help = "vcf including t_ref_count and t_alt_count")
    p$add_argument("cn", help = "CN tsv file with columns: chromosome, start, end, major_1, minor_1")
    p$add_argument("cf", help = "Clonal frequency")

    p$add_argument("pdf", help = "output plot pdf")
    p$add_argument("rdata", help = "output workspace RData")

    return(p$parse_args())
}


main <- function() {
    argv <- get_args()

    vcf <- readVcf(argv$vcf) # vcf path
    clonal_freq <- as.numeric(argv$cf) # purity
    bb <- readCnTable(argv$cn, clonal_freq) # cn_path

    # run MutationTimeR functions
    mt = mutationTime(vcf, bb, n.boot=10) # TODO: clonality
    vcf <- addMutTime(vcf, mt$V)
    mcols(bb) <- cbind(mcols(bb), mt$T)
    
    # save RData
    save(list = ls(all.names = TRUE), file=argv$rdata, 
         envir = environment())

    # plot output
    pdf(argv$pdf, height=8, width=10, useDingbats = FALSE)
    plotSample(vcf, bb)
    dev.off()
}

main()
