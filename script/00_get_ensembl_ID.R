# Script to obtain Ensembl annotations for the NCBI gene IDs

library(gprofiler2)
library(dplyr)

df_raw <- readr::read_delim("data/homo_sapiens_gene_id.txt", col_types='ccc')

# extract the ones we expect to be relevant
idx <- nrow(df_raw)
df_subset <- df_raw[1:idx, ] %>% distinct()
nrow(df_subset)  

df_converted <- gconvert(
  query = df_subset$GeneID, organism = "hsapiens",
  numeric_ns = "ENTREZGENE_ACC", target="ENSG", 
  mthreshold = 1, filter_na = FALSE
)

df_out <- df_subset %>% 
  left_join(select(df_converted, GeneID = input, EnsemblID = target)) %>% 
  filter(!is.na(EnsemblID))
nrow(df_out)

df_out %>% select(GeneID, Symbol, EnsemblID) %>% readr::write_delim("data/homo_sapiens_gene_id_mapping.txt", delim='\t')
