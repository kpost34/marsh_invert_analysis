# Created by Keith Post on 8/20/23
# Association Measures


# Load Packages & DFs===============================================================================
pacman::p_load(here, tidyverse, ade4, adespatial, vegan, gclus, cluster, FD, cowplot, GGally)

source(here("code","01_data-setup.R"))
#once it's determined which is the best data format, then I'll save to RDS and read that in
source(here("code", "00-functions.R"))

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
#binary formula: sqrt(A + B -2J) where the rows of A and B have been normalized to 1



## Log-chord distance matrix
bio_logchord <- dist.ldc(bio, "log.chord")
#alternative coding in vegan
bio_logchord <- log1p(bio) %>%
  decostand(.,"nor") %>%
  dist()
#binary formula: sqrt(A + B -2J) where the rows of A and B have been normalized to 1


## Dissimilarity measures for binary data--------------------------------------------
### Jaccard dissimilarity matrix
#using vegdist()
bio_dj <- vegdist(bio, "jac", binary=TRUE)

#using dist()
bio_dj2 <- dist(bio, "binary")

#using dist.binary()
bio_dj3 <- dist.binary(bio, method=1)

#formula: a/(a + b + c) = # of shared 1s/(# of shared 1s, 0-1s, and 1-0s)


## Sorensen dissimilarity matrix
#using dist.ldc()
bio_ds <- dist.ldc(bio, "sorensen")

#using vegdist()
bio_ds2 <- vegdist(bio, method="bray", binary=TRUE)

#using dist.binary()
bio_ds3 <- dist.binary(bio, method=5)

#formula: 2a/(2a + b + c)


## Ochiai dissimilarity matrix
#using dist.ldc()
bio_och <- dist.ldc(bio, "ochiai")

#using dist.binary()
bio_och <- dist.binary(bio, method=7)

#formula: a/sqrt((a + b)(a + c))



## Graphical display of association matrices-----------------------------------------
### Compare dissimilarity and distance matrices obtained from taxonomic data
#### Percentage difference dissimilarity matrix on raw species abundance data
#hard code to re-recreate coldiss()
bio_db %>% 
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var="site_x") %>% 
  pivot_longer(!site_x, names_to="site_y", values_to="dist") %>% 
  mutate(dist=na_if(dist,0),
         dist=dist/max(dist,na.rm=TRUE),
         across(!dist,~factor(.x, levels=1:23)),
         site_y = fct_rev(site_y)) %>% 
  ggplot() +
  geom_tile(aes(x=site_x, y=site_y, fill=dist),
            color = "white") +
  scale_fill_stepsn(breaks = c(.25, .5, .75, 1),
                    colors = c("orchid3", "orchid1", "darkslategray1", "cyan2"),
                    na.value = "white") +
  labs(x="",
       y="") +
  theme(legend.position="none")

#using functions
plot_assoc_matrix(bio_db)
plot_assoc_matrix(bio_db, ordered=TRUE)

#plot together
c(FALSE, TRUE) %>%
  purrr::map(function(x) {
    plot_assoc_matrix(bio_db,ordered = x)
  }) %>%
  plot_grid(plotlist=.)

#using wrapper function
plot_dissim_grid(bio_db)

#### Using log-transformed data
plot_dissim_grid(bio_dbln)
#distances are lessened considerably relative to untransformed data

#### Chord distance matrix
plot_dissim_grid(bio_dc)

#### Hellinger distance matrix
plot_dissim_grid(bio_dh)
#clusters around sites 3-10 & 12-15

#### Log-chord distance matrix
plot_dissim_grid(bio_logchord)
#clusters around 3-10 & 12-20

#### Comparison
#latter two have smaller relative distances than the chord distance matrix & are similar to each
  #other

#### Jaccard dissimilarity matrix
plot_dissim_grid(bio_dj)
#similar to log-transformed Bray-Curtis matrix...which makes sense because log transformation
  #compresses extremes and may resemble a p-a matrix (e.g., Jaccard)
#however, this matrix is different than the others


#### Simple matching dissimilarity (Sokal and Michener index)
#indicates that taxon absences have shared cause(s) (unlike Jaccard)
bio_s1 <- dist.binary(bio, method=2)
plot_dissim_grid(bio_s1^2)

#compare side-by-side
p1 <- plot_assoc_matrix(bio_dj)
p2 <- plot_assoc_matrix(bio_s1^2)

plot_grid(p1, p2, nrow=1)
#similar patterns overall but some differences noted...
  #1-2, 1-9, 1-15, 3-22, 4-22 (as examples) similar in p2 but not p1...
  #but site 11 more similar to other sites (p1) when completely dissimilar in p2...
  #and 5-15 and 7-15 considered similar in p1 but dissimilar in p2


## Dissimilarity measures for quantitative data (excluding taxon abundances)------
env <- df_2010_env_wide[,-1]

### Euclidean distance matrix of the standardized env DF
env_de <- dist(scale(env))
plot_dissim_grid(env_de)
#seeing similar clusters as some of the distance plots for taxa: ~ sites 3-10, 12-15, and
  #18-22

### Compare Euclidean distance matrix of env data with taxonomic data (Hellinger distance matrix)
p1_bio <- plot_assoc_matrix(bio_dh)
p2_env <- plot_assoc_matrix(env_de)

plot_grid (p1_bio, p2_env, nrow=1)
#commonalities include similarities among sites 3-10 and 12-15 in both and general dissimilarities
  #among sites 2-10 with 17-20

### Euclidean distance matrix on spatial coordinates
#[currently can't simply convert lat-lon to x-ys]


## Dissimilarity measures for binary data (excluding taxonomic p-a data)---------------------
#NA--all data are quantitative



## Dissimilarity measures for mixed types (including categorical variables)------------------
#NA--all data are quantitative



# R Mode============================================================================================
#R mode utilizes correlation-type coefficients to compare variables

## Taxonomic abundance data
#in addition to correlations, the following approach (with Chi-square distances) can be computed
  #on transposed matrices in R mode

#transpose matrix of abundances
bio_t <- t(bio)

#Chi-square pre-transformation followed by Euclidean distance
bio_t_chi <- decostand(bio_t, "chi.square")
bio_t_d16 <- dist(bio_t_chi)
plot_dissim_grid(bio_t_d16, fn=plot_assoc_matrix_taxa)
#Mega, Naem, and Tett do not appear correlated with other groups


## Taxonomic presence-absence data
#Jaccard, Sorensen, Ochiai coefficients can also be used in R mode

#apply S7 to presence-absence data after transposing matrix
bio_t_s7 <- vegdist(bio_t, "jaccard", binary=TRUE)
plot_dissim_grid(bio_t_s7, fn=plot_assoc_matrix_taxa)
#clusters include Brac:Gram and Nema:Tetr


## Quantitative and ordinal data (other than taxonomic abundances)
#covariance or Pearson's r correlation can be used for dimensionally homogenous quantitative
  #variables but these are linear indices and may peform poorly with non-inear relationships. In
  #these cases, Spearman's rho or Kenall's tau may be used

### Pearson r linear correlation among env vars
#### Calculate Pearson corrs
env_pearson <- cor(env) %>%
  round(3) 

env_pearson

#### Pull order
env_o <- order.single(env_pearson)


### Sort correlation matrix
env_sorted <- env_pearson[env_o, env_o]


### Plot correlations using ggpairs
#### Create smoother function
lower_smoother <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "black") +
    geom_smooth(method = method, color = "red", se=FALSE, ...) +
    theme_minimal()
  p
}

#### Plot results
#pearson
ggpairs(env,
        upper = list(continuous = wrap("cor", method = "pearson")),
        diag = list(continuous=wrap("barDiag",bins=8,fill="aquamarine2",color="black")),
        lower = list(continuous=wrap(lower_smoother,method="loess")),
        title = "Pearson Correlation Matrix") +
  theme(plot.title = element_text(hjust=0.5, face="bold"))

#spearman
ggpairs(env,
        upper = list(continuous = wrap("cor", method = "spearman")),
        diag = list(continuous=wrap("barDiag",bins=8,fill="aquamarine2",color="black")),
        lower = list(continuous=wrap(lower_smoother,method="loess")),
        title = "Spearman Correlation Matrix") +
  theme(plot.title = element_text(hjust=0.5, face="bold"))

#kendall
ggpairs(env,
        upper = list(continuous = wrap("cor", method = "kendall")),
        diag = list(continuous=wrap("barDiag",bins=8,fill="aquamarine2",color="black")),
        lower = list(continuous=wrap(lower_smoother,method="loess")),
        title = "Kendall Correlation Matrix") +
  theme(plot.title = element_text(hjust=0.5, face="bold"))



## Pre-transformations for species data
#linear methods (e.g., ANOVA, PCA) which utilize Euclidean distance were previously frowned up for 
  #species data because they failed to properly treat cases of 0-0 matches. Later, pre-transformations
  #were found to be appropriate prior to measuring Euclidean distance: relative abundances by site,
  #chord transformation, Hellinger transformation, chi-square double standardization, and log-chord
  #transformation
#note that these five approaches calculate some form of relative abundance per site, removing
  #the effect of total abundance per site
#all transformations can be done using decostand()











