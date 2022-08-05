## Test out transforming cn to bb
#cn_path = '/work/shah/users/chois7/retreat/test/OV-081/results/OV-081.cn.smooth.tsv'
#cn = read.table(cn_path, header=1)
#ov_bb <- GRanges(
#  seqnames = Rle( as.character(cn$chromosome) ),
#  ranges = IRanges(
#    as.numeric(cn$start),
#    end = as.numeric(cn$end)
#  ),
#  major_cn = as.integer(cn$major_1),
#  minor_cn = as.integer(cn$minor_1),
#  clonal_frequency = rep(0.54, dim(cn)[1]), # TODO
#)

#ov_vcf_path = '/work/shah/users/chois7/retreat/test/OV-081/results/OV-081.addi.vcf'
#ov_vcf = readVcf(ov_vcf_path)

#ov_mt = mutationTime(ov_vcf, ov_bb, n.boot=10)

#ov_vcf <- addMutTime(ov_vcf, ov_mt$V)

#mcols(ov_bb) <- cbind(mcols(ov_bb), ov_mt$T)

#pdf('/work/shah/users/chois7/retreat/test/OV-081/results/OV-081.pdf', 
#    width=10, height=8)
#plotSample(ov_vcf, ov_bb)
#dev.off()


library(argparse)

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
      clonal_frequency = rep(0.54, dim(cn)[1]), # TODO
    )
    return(bb)
}

get_args <- function() {
    p <- ArgumentParser(description = "Run MutationTimeR and plot results")

    p$add_argument("vcf", help = "vcf including t_ref_count and t_alt_count")
    p$add_argument("cn", help = "CN tsv file with columns: chromosome, start, end, major_1, minor_1")
    p$add_argument("cf", help = "Clonal frequency")

    p$add_argument("pdf", help = "output plot pdf")

    return(p$parse_args())
}


main <- function() {
    argv <- get_args()

    vcf <- readVcf(argv$vcf) # vcf path
    bb <- readCnTable(argv$cn) # cn_path
    clonal_freq <- as.numeric(argv$cf)

    # run MutationTimeR functions
    library(MutationTimeR)
    mt = mutationTime(vcf, bb, n.boot=10) # TODO: clonality
    vcf <- addMutTime(vcf, mt$V)
    mcols(bb) <- cbind(mcols(bb), mt$T)
    
    # plot output
    pdf(argv$pdf, height=8, width=10, useDingbats = FALSE)
    plotSample(ov_vcf, ov_bb)
    dev.off()
}

main()


