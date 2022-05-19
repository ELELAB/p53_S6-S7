#!/usr/bin/Rscript

######################
#  RESULTS ANALYSES  #
# Summarise and      #
# illustrate CScape  #
# output             #
######################

#######################
# Author: Freja Dahl Hede, s174333@student.dtu.dk
# Last edit: 1/27/2022
#######################

## Prework
#-------------------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(plotly)

## Data load
#-------------------------------------------------------------------------------
CScapeOutput <- tibble(fread("/data/user/shared_projects/p53_jmb_2021/CScape/results/CScape_results.csv", sep = "\t", skip = "Chromosome"))


cancermutData <- tibble(read.csv("/data/user/shared_projects/p53_jmb_2021/CScape/results/cancermuts_converted_hg19.csv", 
                                 header = TRUE, stringsAsFactors = F))

mutation_MCAP <- tibble(read.csv("/data/user/shared_projects/p53_jmb_2021/M-CAP/results/data_MCAP_scores.csv",
                                 header = TRUE))

## Data wrangle
#-------------------------------------------------------------------------------
# Join meta data onto CScape results
finalData <- CScapeOutput %>% 
  rename(Chromosome = `# Chromosome`) %>% 
  mutate(Reference = as.factor(Reference),
         Mutant = as.factor(Mutant)) %>% 
  full_join(cancermutData, by = c("Chromosome" = "chromosome",
                                  "Position" = "position",
                                  "Reference" = "ref.base",
                                  "Mutant" = "mut.base")) %>% 
  unique() %>% 
  left_join(mutation_MCAP, by = c("Chromosome" = "chromosome",
                                  "Position" = "position",
                                  "Reference" = "ref.base",
                                  "Mutant" = "mut.base", 
                                  "WT.residue",
                                  "Mutated.residue",
                                  "prot.Position",
                                  "cancermuts.index")) %>% 
  unique() %>% 
  mutate(MCAP_status_raw = case_when(mcapv1.4 >= 0.025 ~ "pathogenic",
                                mcapv1.4 < 0.025 ~ "possibly benign"),
         MCAP_status_sens = case_when(mcap_sensitivityv1.4 <= 0.95 ~ "pathogenic",
                                 mcap_sensitivityv1.4 > 0.95 ~ "possibly benign"))


## Summarise outputs
#-------------------------------------------------------------------------------
# CScape, nucleotide level, protein level
finalData %>% 
  group_by(Remark) %>% 
  summarise(count = n())

finalData %>% 
  unite(Protein, c(prot.Position, WT.residue, Mutated.residue), sep=",") %>% 
  group_by(Protein) %>% 
  filter(Coding == min(Coding)) %>% 
  select(Protein, Remark) %>% 
  unique() %>% 
  group_by(Remark) %>% 
  summarise(count = n())

# MCAP, nucleotide level, protein level
finalData %>%
  group_by(MCAP_status_raw) %>%
  summarise(count = n())

finalData %>%
  group_by(MCAP_status_sens) %>%
  summarise(count = n())  

finalData %>% 
  unite(Protein, c(prot.Position, WT.residue, Mutated.residue), sep=",") %>% 
  group_by(Protein) %>% 
  filter(mcap_sensitivityv1.4 == min(mcap_sensitivityv1.4)) %>% 
  select(Protein, MCAP_status_sens) %>% 
  unique() %>% 
  group_by(MCAP_status_sens) %>% 
  summarise(count = n())

# Revel: Only protein level
finalData %>%
  unite(Protein, c(prot.Position, WT.residue, Mutated.residue), sep=",") %>% 
  select(Protein, Revel_status) %>% 
  unique() %>% 
  group_by(Revel_status) %>%
  summarise(count = n())  

## Illustrate outputs
#-------------------------------------------------------------------------------

p_CScape <- finalData %>% 
  unite(mutation, c(Chromosome, Position, Reference, Mutant), sep=",") %>%
  unite(Protein, c(prot.Position, WT.residue, Mutated.residue), sep=",") %>% 
  #mutate(mutation = fct_reorder(Protein, desc(Coding))) %>%
  ggplot(mapping = aes(x = Protein, y = Coding, color = Remark)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_manual(values=c("High-confidence" = "red", "Low-confidence" = "blue", "Neutral" = "forestgreen")) +
  #geom_hline(yintercept = MCAP_threshold, linetype = "dashed", color = "red") +
  xlab("Protein mutation") + ylab("Driver score") +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 
  #theme(axis.text.x = element_text(angle = 90, size = 2)) 

p_CScape

p_MCAP <- finalData %>% 
  unite(mutation, c(Chromosome, Position, Reference, Mutant), sep=",") %>%
  unite(Protein, c(prot.Position, WT.residue, Mutated.residue), sep=",") %>% 
  #mutate(mutation = fct_reorder(Protein, desc(Coding))) %>%
  ggplot(mapping = aes(x = Protein, y = mcap_sensitivityv1.4, color = MCAP_status_sens)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_manual(values=c("pathogenic" = "red", "possibly benign" = "forestgreen")) +
  #geom_hline(yintercept = MCAP_threshold, linetype = "dashed", color = "red") +
  xlab("Protein mutation") + ylab("M-CAP sensitivity score") +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 
#theme(axis.text.x = element_text(angle = 90, size = 2)) 

p_MCAP


p_Revel <- finalData %>% 
  unite(mutation, c(Chromosome, Position, Reference, Mutant), sep=",") %>%
  unite(Protein, c(prot.Position, WT.residue, Mutated.residue), sep=",") %>% 
  #mutate(mutation = fct_reorder(Protein, desc(Coding))) %>%
  ggplot(mapping = aes(x = Protein, y = Revel.score, color = Revel_status)) +
  geom_point(size = 0.6, alpha = 0.8) +
  #geom_hline(yintercept = MCAP_threshold, linetype = "dashed", color = "red") +
  scale_color_manual(values=c(pathogenic = "red", neutral = "forestgreen")) +
  xlab("Protein mutation") + ylab("Revel score") +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 
  #theme(axis.text.x = element_text(angle = 90, size = 2)) 

p_Revel 


## Wrangle final data
#-------------------------------------------------------------------------------
# Write output file
cancermutWScores <- finalData %>% 
  group_by(cancermuts.index, prot.Position, WT.residue, Mutated.residue) %>% 
  summarise(CScape.coding.score = paste(sort(unique(Coding)), collapse=", "),
            MCAP.sensitivity.score = paste(sort(unique(mcap_sensitivityv1.4)), collapse=", "),
            Revel.score) %>% unique() %>% arrange(cancermuts.index)

cancermutWScoresDBD <- cancermutWScores %>% 
  filter(prot.Position >= 91 & prot.Position <= 289) 

## Write output file
#-------------------------------------------------------------------------------
write_csv(cancermutWScores, "/data/user/shared_projects/p53_jmb_2021/CScape/results/cancermuts_CScape_MCAP_Revel.csv")
write_csv(cancermutWScoresDBD, "/data/user/shared_projects/p53_jmb_2021/CScape/results/cancermutsDBD_CScape_MCAP_Revel.csv")





