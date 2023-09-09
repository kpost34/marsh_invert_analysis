# Created by Keith Post on 9/1/23
# Cluster Analysis


# Load Packages & DFs===============================================================================
pacman::p_load(here, tidyverse, ade4, adespatial, vegan, gclus, cluster, pvclust, RColorBrewer,
               labdsv, rioja, indicspecies, dendextend, vegclust, 
               colorspace, agricolae, picante, ggdendro, cowplot)

source(here("code","01_data-setup.R"))
#once it's determined which is the best data format, then I'll save to RDS and read that in
source(here("code", "00-functions.R"))

bio <- df_2010_bio_wide[,-1]


# Hierarchical Clustering Based on Links============================================================
## Single linkage agglomerative clustering------------
#objects are grouped by highest pairwise similarity, and subsequent objects (groups) agglomerate
  #to the first group through the closest pair of objects, which results in a chaining pattern of 
  #clusters
#aka "nearest neighbor sorting"


### Compute matrix of chord distance among sites
bio_ch <- decostand(bio, "normalize") %>% vegdist("euc")
#same as bio_dc <- decostand(bio, "nor") %>% dist() from 03

### Compute single linkage agglomerative clustering
bio_ch_single <- hclust(bio_ch, method="single")


### Plot dendogram
#### Using plot()
plot(bio_ch_single,
     labels=rownames(bio),
     main="Chord - Single linkage")


### Using {ggendro}
#### Using ggdendro() with defaults
ggdendrogram(bio_ch_single)


#### With customizations (hard-coded)
#extract dendrogram plot data
bio_ch_single_tree <- bio_ch_single %>% #hclust object
  as.dendrogram() %>% #converts to dendrogram object
  hang.dendrogram() %>% #hangs leaves of dendrogram
  dendro_data(type="rectangle") #extracts plotting data

#determine yend of each leaf and join with site #
site_labs <- segment(bio_ch_single_tree) %>%
  filter(xend %in% 1:23) %>%
  group_by(xend) %>%
  filter(yend==min(yend)) %>% 
  ungroup() %>%
  mutate(yend=yend - 0.01) %>%
  select(ends_with("end")) %>%
  left_join(
    label(bio_ch_single_tree)[,c("x", "label")],
    by=c("xend"="x")
  )


bio_ch_single_tree %>%
  #extracts x-y data from list-object
  segment() %>%
  ggplot() +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data=site_labs,
            aes(x=xend, y=yend, label=label)) +
  coord_cartesian(ylim=c(0, NA)) +
  labs(x="site",
       y="height",
       title="Chord - Single Linkage") +
  theme_void() + #removes axes, axis.titles, and axis.text
  #specify margins as well as center and bold plot title & rotate y-axis title 90o
  theme(plot.title=element_text(hjust=0.5, face="bold", margin=margin(t=10, b=10)),
        axis.line.y = element_line(),
        axis.title.x=element_text(margin=margin(b=5)),
        axis.title.y=element_text(angle=90, margin=margin(r=10, l=5)),
        axis.text.x=element_blank(),
        # axis.text.x=element_text(margin=margin(t=5)),
        axis.text.y=element_text(margin=margin(r=10)))


#### Using function
plot_dendro(bio_ch_single, "Chord - Single Linkage", .01)
#interpretation: 21, 20, and 23 appear to be different from the rest; similar idea with 17 & 18,
  #which appear similar to each other but different from the rest
#three larger clusters seem apparent: 3-1-5-6, 16-9-22, and 19-2-11-15-13-10-14-4-8

p1 <- plot_dendro(bio_ch_single, "Chord - Single Linkage", .01)


## Complete linkage agglomerative clustering-----
#aka "furthest neighbor sorting"; one object/group groups with another group at a similarity/
  #dissimilarity of the most dissimilar objects, ensuring that all members of both groups 
  #are connected

## Compute complete linkage agglomerative clustering
bio_ch_complete <- hclust(bio_ch, method="complete")

### Plot dendrogram
#### Using plot()
plot(bio_ch_complete,
     labels=rownames(bio),
     main="Chord - Complete linkage")


### Using {ggdendro}
#### Using ggdendro() with defaults
ggdendrogram(bio_ch_complete)


#### With customization
plot_dendro(bio_ch_complete, "Chord - Complete Linkage", .01)
#interpretation: unsurprisingly, I see more differentiation than with single linkage clustering.
  #'closely related' clusters are much smaller--typically 2 sites or maybe three--and lots of
  #differentiation (branching) at greater heights

p2 <- plot_dendro(bio_ch_complete, "Chord - Complete Linkage", .02)

#compare the two methods
plot_grid(p1, p2)



## Average agglomerative clustering-----
#four methods based on average dissimilarities. Unweighted pair-group method using arithmetic 
  #averages (UPGMA) is the most known: groups join by mean dissimilarities between an object and 
  #all members of a group. 

## Compute UPGMA agglomerative clustering
bio_ch_UPGMA <- hclust(bio_ch, method="average")

## Plot dendrogram
### Using plot()
plot(bio_ch_UPGMA,
     labels=rownames(bio),
     main="Chord - UPGMA")

### Using {ggdendro}
plot_dendro(bio_ch_UPGMA, "Chord - UPGMA", .015)


# Ward's Minimum Variance Clustering================================================================
#groups are identified by minimizing within-group sum of squares

## Compute Ward's minimum variance clustering
bio_ch_ward <- hclust(bio_ch, method="ward.D2")

## Plot dendrogram
### Using plot()
plot(bio_ch_ward,
     main="Chord - Ward")

### Using {ggdendro}
p3 <- plot_dendro(bio_ch_ward, "Chord - Ward", .02)
#interpretation: a lot more branching occurs towards the bottom of the tree, yielding multiple
  #larger groups--7-1-3-5-6, 16-9-22, 19-15-13-10-14, 2-11-12-4-8--and smaller groups--20-23 and
  #17-18--and site 21 by itself


# Flexible Clustering===============================================================================
#clustering method with four parameters--alpha-h, alpha-i, beta, and gamma--which can be adjusted
  #to generate the aforementioned clustering methods

## Compute beta-flexible clustering using cluster::agnes()
# beta = -0.25
bio_ch_beta2 <- agnes(bio_ch, method="flexible", par.method=0.625)

## Change the class of agnes object
class(bio_ch_beta2) #agnes, twins
bio_ch_beta2 <- as.hclust(bio_ch_beta2) 
class(bio_ch_beta2) #hclust

## Plot dendrogram
### With plot()
plot(bio_ch_beta2,
     labels=rownames(bio),
     main="Chord - Beta-flexible (beta=-0.25)")

### With {ggdendro}
p4 <- plot_dendro(bio_ch_beta2, "Chord - Beta-flexible (beta=-0.25)", .02)


## Compare this method with Ward's Minimum Variance Clustering
plot_grid(p3, p4)
#when considering the possibility of flipping the tree on nodes, these trees are nearly identical:
  #the pairs at the bottom and the larger groups and structure are nearly the same. Site 21 is
  #placed differently in the two trees.


# Interpreting and Comparing Hierarchical Clustering Results========================================
## To extract more information from the clustering result, use summary(cluster.obj)
summary(bio_ch_beta2)
summary(bio_ch_ward)


## Cophenetic correlation-----
#cophenetic distance: distance where two objects become members of the same group. On a dendrogram,
  #this value can be found by identifying two objects, locating the node that separates them,
  #and determining the height (distance) at that node

#cophentic matrix: matrix of cophenetic distances among all pairs of objects

#cophenetic correlation: Pearson's r correlation between the original dissimilarity matrix and the
  #cophenetic matrix; highest cophenetic corr = method that retains most info from dissimilarity
  #matrix

### Compute cophenetic matrix and correlation of four clustering results
#### Single linkage clustering
bio_ch_single_coph <- cophenetic(bio_ch_single)
cor(bio_ch, bio_ch_single_coph) #0.892
#note that significance (if used cor.test()) is uninterpretable
  #because the two sets of data are *not* independent


#### Complete linkage clustering
bio_ch_comp_coph <- cophenetic(bio_ch_complete)
cor(bio_ch, bio_ch_comp_coph) #0.711


#### Average clustering
bio_ch_UPGMA_coph <- cophenetic(bio_ch_UPGMA)
cor(bio_ch, bio_ch_UPGMA_coph) #0.917


#### Ward clustering
bio_ch_ward_coph <- cophenetic(bio_ch_ward)
cor(bio_ch, bio_ch_ward_coph) #0.835

#interpretation: average clustering using UPGMA retains the closest relationship to the chord
  #distance matrix


#### Ward clustering (using Spearman method)
cor(bio_ch, bio_ch_ward_coph, method="spearman") #0.734



### Visualize cophenetic correlations vs dissimilarity matrices
#create vector and list of types and matrices, respectively
types <- c(NA, "single_coph", "complete_coph", "upgma_coph", "ward_coph")

list_mat_dist_coph <- list(bio_ch, bio_ch_single_coph, bio_ch_comp_coph, bio_ch_UPGMA_coph, bio_ch_ward_coph)
list_mat_coph <- list_mat_dist_coph[-1]

#create a crosswalk for titling
method_cwalk <- tibble(method = paste0(types[!is.na(types)],"_dist"),
                       title = c("Single linkage",
                                 "Complete linkage",
                                 "UPGMA",
                                 "Ward"))

### named_vector (for labeller)
method_labeller <- list_mat_coph  %>%
  purrr::map_dbl(cor, bio_ch) %>%
  round(3) %>%
  paste(method_cwalk$title, "\nCophenetic correlation =", .) %>%
  set_names(method_cwalk$method) %>%
  as_labeller()


#feed types & list_matrix_clust into map2 & internal join
map2(.x=list_mat_dist_coph, .y=types, .f=convert_to_long_df) %>%
  reduce(inner_join) %>% 
  pivot_longer(cols=contains("coph"),
               names_to="method",
               values_to="coph_dist") %>%
  left_join(method_cwalk) %>%
  mutate(method=factor(method, levels=method_cwalk$method)) %>%
  ggplot(aes(x=dist, y=coph_dist)) +
  geom_point(shape=1, alpha=0.6) +
  geom_abline(slope=1, intercept=0) +
  geom_smooth(method="loess", color="red", linewidth=0.5, se=FALSE) +
  facet_wrap(~method, labeller=method_labeller) +
  labs(x="Chord distance",
       y="Cophenetic distance") +
  theme_bw() +
  theme(strip.text=element_text(face="bold", size=13),
        strip.background=element_blank())
#interpretation: the smoother that's closest to the 1:1 diagonal yields the greatest correlation,
  #as evidenced by UPGMA's smoother styaing with the line for ~2/3 of the data followed by single
  #linkage which deviates slightly below the line


### Gower distances
#Gower distances, calculated as the sum of squared differences between the original dissimilarities
  #and cophenetic distances, is another method to compare clustering results. Here the smallest
  #Gower distance is evidence of the best clustering result

list_mat_coph %>%
  purrr::map_dbl(function(x) {
    sum((bio_ch - x)^2)
  }) %>%
  set_names(method_cwalk$title)
#interpretation: UPGMA has the smallest Gower distance followed by single linkage (which mirrors
  #the cophenetic correlation result)


## Looking for interpretable clusters-----
#need to decide where to cut a dendrogram to generate clusters. Clusters can 
  #be cut at one or more heights and by using subjectivity or other criteria

### Graph of the fusion levels
#dissimilarity values where a fusion between two branches of a dendrogram occurs. A plot
  #of the fusion levels can help to identify the cutting height

#### Create DF of cluster number vs height
nm_method <- c("Complete", "UPGMA", "Ward", "Beta-flexible")

df_n_clust_height <- list(bio_ch_complete, bio_ch_UPGMA, bio_ch_ward, bio_ch_beta2) %>%
  set_names(nm_method) %>%
  purrr::map_df("height") %>%
  bind_cols(n_clust=23:2,.) %>%
  pivot_longer(!n_clust, names_to="method", values_to="height") %>%
  arrange(n_clust)


#### Plot fusion levels
nm_method_strip <- paste0("Fusion levels - Chord - ", nm_method)


df_n_clust_height %>%
  mutate(method=paste0("Fusion levels - Chord - ", method),
         method=factor(method, levels=nm_method_strip)) %>%
  ggplot(aes(x=height, y=n_clust)) +
  geom_step(color="gray50") +
  geom_text(aes(label=n_clust), size=3, color="red") +
  labs(x="h (node height)",
       y="k (number of clusters)") +
  facet_wrap(~method, scales="free_x") +
  theme_bw() +
  theme(strip.text=element_text(face="bold", size=12),
        strip.background=element_blank())

#interpretation: read from right to left and recognize long horizontal lines which represent cutting
  #areas: complete: 6, UPGMA: 5, Ward: 2 (or 6); and Beta-Flexible: 5
#look at dendros and apply cluster number cut to them

#### Grab heights at these cuts
fuse_h <- df_n_clust_height %>% 
  filter( #grab n_clust from fusion levels plots
    n_clust==5 & method %in% c("Complete", "Beta-flexible")|
      n_clust==6 & method=="UPGMA"|
        n_clust==2 & method=="Ward"
  ) %>%
  pull(height, name="method") %>%
  .[nm_method] #sort


#### Assign dendros names
p_fuse1 <- plot_dendro(bio_ch_complete, "Chord - Complete Linkage", .02)
p_fuse2 <- plot_dendro(bio_ch_UPGMA, "Chord - UPGMA", .015)
p_fuse3 <- plot_dendro(bio_ch_ward, "Chord - Ward", .02)
p_fuse4 <- plot_dendro(bio_ch_beta2, "Chord - Beta-flexible (beta=-0.25)", .02)

list(p_fuse1, p_fuse2, p_fuse3, p_fuse4) %>%
  purrr::map2(.y = fuse_h, function(x, y) {
    x +
      geom_hline(yintercept=y, linetype=2, color="red") +
      geom_text()
  }) %>%
  plot_grid(plotlist=., nrow=2)

  














     

