#!/usr/bin/Rscript

#######################
#   DATA CONVERSION   #
# Converts aggregation#
# output to format    #
# that cSCAPE accepts #
#######################

#######################
# Author: Freja Dahl Hede, s174333@student.dtu.dk
# Last edit: 11/24/2021
#######################

## Prework
#-------------------------------------------------------------------------------
library(tidyverse)
library(liftOver)

select <- dplyr::select
filter <- dplyr::filter

## Data load
#-------------------------------------------------------------------------------
cancermutData <- tibble(read.csv("../aggregation/jan2022/metatable.csv", 
                            header = TRUE))


## Data wrangle
#-------------------------------------------------------------------------------
# Select relevant columns and make row for each genomic mutation
data <- cancermutData %>% 
  select(X,
         Position, 
         WT.residue, 
         Mutated.residue,
         Revel.score,
         Genomic.mutation) %>% 
  separate_rows(Genomic.mutation, sep = ", ") # Convert list into sep rows

## Does hgvs DNA format cover only SNVs?
#-------------------------------------------------------------------------------
# Regex: "hg\\d{2}[^:]*:g\\.\\d+[ATCG]>[ATCG]"
# Extraxt all formats not following hgvs DNA format (or empty)
filter(data, 
       !str_detect(Genomic.mutation, "hg\\d{2}[^:]*:g\\.\\d+[ATCG]>[ATCG]") & Genomic.mutation != "") %>% 
  print(n = 100L)

# Conclusion: No, SNVs are covered by the substitution format

## Data wrangle
#-------------------------------------------------------------------------------
data_wrangled <- data %>% filter(str_detect(Genomic.mutation, "hg\\d{2}[^:]*:g\\.\\d+[ATCG]>[ATCG]")) %>% 
  mutate(genome = str_extract(Genomic.mutation, "hg\\d{2}"),
         chromosome = str_sub(str_extract(Genomic.mutation, ",[^:|\\.]*"), start=2),
         position = as.integer(str_sub(str_extract(Genomic.mutation, ":g\\.\\d+"), start=4)),
         ref.base = str_sub(str_extract(Genomic.mutation, "[ATCG]>"), end=-2),
         mut.base = str_sub(str_extract(Genomic.mutation, ">[ATCG]"), start=2),
	 prot.Position = Position,
	 cancermuts.index = X) %>% 
  select(Genomic.mutation, genome, chromosome, position, ref.base, mut.base, WT.residue, Mutated.residue, prot.Position, Revel.score, cancermuts.index) 
         
  
# hg38,17:g.7676587T>C
# Convert from hg39 to hg19 genome assembly
hg38_data <- data_wrangled %>% dplyr::filter(genome == "hg38")
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
grs = GRanges(seqnames=str_c("chr", dplyr::pull(hg38_data, chromosome)), 
              ranges=IRanges(start=dplyr::pull(hg38_data, position), end=dplyr::pull(hg38_data, position)), 
              strand=rep('+', dim(hg38_data)[1]))

seqlevelsStyle(grs) = "UCSC"  # https://bioconductor.org/packages/devel/workflows/vignettes/liftOver/inst/doc/liftov.html
convs = liftOver(grs, ch)

# Add converted chr and position 
convs_19 <- as.data.frame(unlist(convs))
hg38_to_hg19_data <- hg38_data %>% 
  mutate(chromosome = str_replace(pull(convs_19, seqnames), "chr", ""),
         position = pull(convs_19, start), 
         genome = "hg19")

# Print failed conversions if any
message("### The percentage of failed hg38 --> hg19 conversions: ",
        (dim(hg38_data)[1]-dim(convs_19)[1])/dim(hg38_data)[1]*100,
        " %")

# Add back to original data
output_data <- data_wrangled %>% 
  filter(genome != "hg38") %>% 
  bind_rows(hg38_to_hg19_data) %>% 
  select(-Genomic.mutation, -genome) %>%
  unique()

# Write output
write_csv(output_data, "/data/user/shared_projects/p53_jmb_2021/CScape/results/cancermuts_converted_hg19.csv")
