############################################################
# Code for obtaining Unique Peptides from Proteome Discoverer
###########################################################
library(dplyr)
library(tidyr)
library(writexl)


file_path <- "your_file.xlsx"
sheet_name <- "PeptideGroups"  
raw <- readxl::read_xlsx(file_path, sheet = sheet_name) %>% as.data.frame()

raw <- as.data.frame(readxl::read_xlsx("25-18_PD3.2_LF_2a.replica_tecnica-DB_Tgallinae-transcript-pepts.xlsx", sheet = "PeptideGroups"))

# Identify abundance columns dynamically
abundance_cols <- grep("Abundance:", colnames(raw), value = TRUE)

# Build peptide dataframe
df_peptides <- raw %>%
  filter(!grepl("Oxidation", Modifications)) %>%     # remove Oxidation-modified peptides
  distinct(Sequence, .keep_all = TRUE) %>%           # keep unique sequences
  select(Sequence, `Master Protein Accessions`, all_of(abundance_cols)) %>%
  pivot_longer(cols = all_of(abundance_cols),
               names_to = "Sample",
               values_to = "Abundance") %>%
  filter(!is.na(Abundance)) %>%
  group_by(`Master Protein Accessions`, Sample) %>%
  summarise(Peptide_Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Peptide_Count, values_fill = 0) %>%
  rename_with(~ paste0("Unique peptides ", .), -`Master Protein Accessions`)

proteins <- as.data.frame(readxl::read_xlsx("25-18_PD3.2_LF_2a.replica_tecnica-DB_Tgallinae-transcript-prots.xlsx", sheet = "Proteins"))

filtered_df <- df_peptides %>%
  inner_join(proteins, by = c("Master Protein Accessions" = "Accession"))

write_xlsx(filtered_df, path = "Protein_groups_Unique_Peptides.xlsx")
