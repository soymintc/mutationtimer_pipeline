library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
rdatafile = args[1]
clone_id_path = args[2]

data = readRDS(rdatafile)
cldf = data$cl$clustering # cluster dataframe
clone_id <- cldf %>% select(cell_id, clone_id)
write.csv(clone_id, clone_id_path, row.names=F, quote=F, sep="\t")
