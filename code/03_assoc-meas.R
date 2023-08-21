# Created by Keith Post on 8/20/23
# Association Measures


# Load Packages & DFs===============================================================================
pacman::p_load(here, tidyverse, ade4, adespatial, vegan, gclus, cluster, FD)

source(here("code","01_data-setup.R"))
#once it's determined which is the best data format, then I'll save to RDS and read that in


bio <- df_2010_bio_wide[,-1]

# Q Mode============================================================================================
## Dissimilarity and distance measures for quantitative data--------------------------
### Percentage difference (aka Bray-Curtis) dissimilarity matrix on raw taxonomic data
bio_db <- vegdist(bio)
bio_db[1:6]
#binary formula: (A + B -2J)/(A + B)


### Percentage difference (aka Bray-Curtis) dissimilarity matrix on log-transformed abundances
bio_dbln <- vegdist(log1p(bio))
bio_dbln[1:6]
#binary formula: (A + B -2J)/(A + B)


### Chord distance matrix
bio_dc <- dist.ldc(bio, "chord")
#alternative coding in vegan
bio_dc <- decostand(bio, "nor") %>% dist()
bio_dc[1:6]
#binary formula: sqrt(A + B -2J) where the rows of A and B have been normalized to 1


### Hellinger distnace matrix
bio_dh <- dist.ldc(bio)
#alternative coding in vegan
bio_dh <- decostand(bio, "hel") %>% dist()
#binary formula: sqrt(A + B -2J) where the rows of A and B have been normalized to 1 & the 
  #square root of those values is taken



## Log-chord distance matrix
bio_logchord <- dist.ldc(bio, "log.chord")
#alternative coding in vegan
bio_logchord <- log1p(bio) %>%
  decostand(.,"nor") %>%
  dist()
#binary formula: sqrt(A + B -2J) where sqrt is applied to A & B prior to their rows
  #being normalized to 1


## Dissimilarity measures for binary data--------------------------------------------
### Jaccard dissimilarity matrix
#using vegdist()
bio_dj <- vegdist(bio, "jac", binary=TRUE)

#using dist()
bio_dj2 <- dist(bio, "binary")

#formula: 2B/(1 + B) where B is Bray-Curtis dissimilarity




# R Mode============================================================================================






