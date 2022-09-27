library(argparse)
library(tidyverse)

parser <- ArgumentParser(description = "Calculate consensus clone copy number")
parser$add_argument('--hscn', type='character', 
                    help="Haplotype-specific copy number object")
parser$add_argument('--cna_bins_consensus', type = 'character',
                    help="Consensus copy number bins")
parser$add_argument('--cna_bins_consensus_mutationtimer', type = 'character',
                    help="Consensus copy number bins")
parser$add_argument('--purity_ploidy', type = 'character',
                    help="Purity and ploidy")
args <- parser$parse_args()

# Read allele-specific copy number
signals_rds_obj <- readr::read_rds(args$hscn)

# Find consensus copy number
cna_bins_consensus <- signals::consensuscopynumber(signals_rds_obj$hscn$data)

# Convert to mutationtimer input
cna_bins_consensus_mutationtimer <- cna_bins_consensus %>%
    select(
        chromosome=chr,
        start=start,
        end=end,
        major_1=Maj,
        minor_1=Min
    )

# Calculate aggregate ploidy across cells
find_purity <- function(cna_bins){
    message('Removing Y chromosome')
    cna_bins <- cna_bins %>%
        dplyr::filter(chr != 'Y')

    # Find the fraction of each cell that is diploid and fraction that == cell ploidy
    cna_bins_ploidy <- cna_bins %>%
        group_by(cell_id) %>%
        summarize(
            frac_nondiploid = sum(state != 2) / n(), 
            frac_cell_ploidy = sum(state == signals:::Mode(state)) / n(), 
        ) %>%
        ungroup

    # Find cells that are probably diploid or that are diploid + misscalled ploidy 
    diploid_cutoff <- 0.05

    message(paste0("Diploid cutoff: ", round(diploid_cutoff, 3)))

    # Compute outlier cutoff
    frac_cell_ploidy_cutoff <- cna_bins_ploidy %>%
        # Remove cells that are probably diploid or that are diploid + misscalled ploidy 
        filter(frac_cell_ploidy < (1 - diploid_cutoff)) %>%
        summarize(cutoff = mean(frac_cell_ploidy) + 3 * sd(frac_cell_ploidy)) %>%
        pull(cutoff)

    message(paste0("Fraction ploidy cutoff: ", round(frac_cell_ploidy_cutoff, 3)))

    cell_type <- cna_bins_ploidy %>% 
        mutate(
            cell_type = ifelse(
                frac_nondiploid > diploid_cutoff & frac_cell_ploidy < frac_cell_ploidy_cutoff,
                "Tumor",
                "Normal"
            )
        ) %>%
        select(cell_id, cell_type)

    purity <- cell_type %>%
        group_by(cell_type) %>%
        summarize(
            num_cells = n_distinct(cell_id)
        ) %>%
        mutate(
            purity = num_cells / sum(num_cells)
        ) %>%
        ungroup %>%
        filter(cell_type == "Tumor") %>%
        pull(purity)

    return(purity)

}

purity <- find_purity(signals_rds_obj$hscn$data)

# Calculate aggregate ploidy across cells
find_ploidy <- function(cna_bins){
    message('Removing Y chromosome')
    cna_bins <- cna_bins %>%
        dplyr::filter(chr != 'Y')

    cell_ploidy <- cna_bins %>%
        group_by(cell_id) %>%
        dplyr::summarise(
            ploidy_mode = signals:::Mode(state),
            ploidy_mean = mean(state, na.rm = TRUE),
            ploidy_median = median(state, na.rm = TRUE)
        ) %>%
        ungroup() %>%
        dplyr::summarise(
            ploidy_mean = mean(ploidy_mean, na.rm = TRUE)
        ) %>%
        pull(ploidy_mean)
}

ploidy <- find_ploidy(signals_rds_obj$hscn$data)

purity_ploidy <- tibble(
    purity = purity,
    ploidy = ploidy
)

# Export copy number consensus
readr::write_tsv(cna_bins_consensus, args$cna_bins_consensus)
readr::write_tsv(cna_bins_consensus_mutationtimer, args$cna_bins_consensus_mutationtimer)

# Export purity and ploidy
readr::write_csv(purity_ploidy, args$purity_ploidy)
