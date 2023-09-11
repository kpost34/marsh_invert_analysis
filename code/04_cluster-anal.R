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
site_cwalk <- df_2010_wide[,c("site", "loc")]


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


#### Plot fusion levels & dendrograms with cut heights
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

#interpretation: read from right to left and identify h lines preceding large increases which represent 
  #cutting areas: complete: 7, UPGMA: 5, Ward: 6; and Beta-Flexible: 8
#look at dendros and apply cluster number cut to them

#### Grab heights at these cuts
#simple version: height only as vector
fuse_h <- df_n_clust_height %>% 
  filter( #grab n_clust from fusion levels plots
    n_clust==7 & method=="Complete"|
      n_clust==5 & method=="UPGMA"|
        n_clust==6 & method=="Ward"|
          n_clust==8 & method=="Beta-flexible"
  ) %>%
  pull(height, name="method") %>%
  .[nm_method] #sort

#complex version: df of height and n_clust
list_n_clust_h <- df_n_clust_height %>% 
  filter( #grab n_clust from fusion levels plots
    n_clust==7 & method=="Complete"|
      n_clust==5 & method=="UPGMA"|
        n_clust==6 & method=="Ward"|
          n_clust==8 & method=="Beta-flexible"
  ) %>%
  mutate(method=factor(method, levels=nm_method)) %>%
  arrange(method) %>%
  split(1:nrow(.))
  


#### Assign dendros names
p_fuse1 <- plot_dendro(bio_ch_complete, "Chord - Complete Linkage", .02)
p_fuse2 <- plot_dendro(bio_ch_UPGMA, "Chord - UPGMA", .015)
p_fuse3 <- plot_dendro(bio_ch_ward, "Chord - Ward", .02)
p_fuse4 <- plot_dendro(bio_ch_beta2, "Chord - Beta-flexible (beta=-0.25)", .02)


#### Plot dendros with cut lines
#simple version
list(p_fuse1, p_fuse2, p_fuse3, p_fuse4) %>%
  purrr::map2(.y = fuse_h, function(x, y) {
    x +
      geom_hline(yintercept=y, linetype=2, color="red")
  }) %>%
  plot_grid(plotlist=., nrow=2)

#complex version
list(p_fuse1, p_fuse2, p_fuse3, p_fuse4) %>%
  purrr::map2(.y = list_n_clust_h, function(dendro, lab) {
    dendro +
      geom_hline(yintercept=lab[["height"]], linetype=2, color="red") +
      geom_text(aes(x=22, y=lab[["height"]], label=paste("k =", lab["n_clust"])))
  }) %>%
  plot_grid(plotlist=., nrow=2)
#interpretation:
  #complete: k = 7, two groups of 1, two groups of 2, group of 3, group of 4, and one large group;
  #UPGMA: k = 5, group of 1, two groups of 2, group of 5, and one large group;
  #Ward: k = 6, group of 1, two groups of 2, and 3 medium-sized groups;
  #beta-flexible: k = 8, two groups of 1, two groups of 2, group of 3, group of 4, and two 
    #medium sized groups

#interpretation: [add map then describe]
  

### Set common number of groups and compare group contents among dendrograms
#### First choose k (group number) where there is at least a small jump in all four fusion levels
k <- 6 #best value using all four methods

#### Cut dendros and create membership vectors
bioch_single_g <- cutree(bio_ch_single, k=k)
bioch_complete_g <- cutree(bio_ch_complete, k=k)
bioch_UPGMA_g <- cutree(bio_ch_UPGMA, k=k)
bioch_ward_g <- cutree(bio_ch_ward, k=k)
bioch_beta_g <- cutree(bio_ch_beta2, k=k)

#### Compare classification by constructing contingency tables
table(bioch_single_g, bioch_complete_g) #3 classifications differ
table(bioch_single_g, bioch_UPGMA_g) #0 classifications differ
table(bioch_single_g, bioch_ward_g) #4 classifications differ
table(bioch_complete_g, bioch_UPGMA_g) #3 classifications differ
table(bioch_complete_g, bioch_ward_g) #1 classification differs
table(bioch_UPGMA_g, bioch_ward_g) #4 classifications differ
table(bioch_beta_g, bioch_ward_g) #0 classifications differ

#single-UPGMA and beta-Ward are identical and complete-Ward are nearly identical


### Compare two dendrograms to highlight common subtrees
#compare dendros and seek common clusters to help select a partition found by multiple clustering
  #methods

#### Objects of class "hclust" must be first converted into objects of class "dendrogram"
#### Using dendextend::tanglegram()
class(bio_ch_ward)
dend1 <- as.dendrogram(bio_ch_ward)
class(dend1)
dend2 <- as.dendrogram(bio_ch_complete)
dend12 <- dendlist(dend1, dend2)
tanglegram(
  untangle(dend12),
  sort=TRUE,
  common_subtrees_color_branches=TRUE,
  main_left="Ward method",
  main_right="Complete linkage"
)

#interpretation: robust cluters: 1-3-5-6, 2-11, and 4-8-12


#### Using ggendrogram
#complete--l-r
# bio_ch_comp_data <- dendro_data(bio_ch_complete, "rectangle")
# 
# p_comp_lr <- ggplot() +
#   geom_segment(data=segment(bio_ch_comp_data),
#                aes(x=x, y=y, xend=xend, yend=yend)) +
#   geom_text(data=label(bio_ch_ward_data),
#             aes(x=x, y=y, label=label, hjust=0)) +
#   coord_flip() +
#   scale_y_reverse() +
#   theme_void() +
#   theme(axis.text.x=element_text()) 
# 
# 
# #ward--r-l
# # bio_ch_ward_data <- dendro_data(bio_ch_ward, type="rectangle")
# bio_ch_ward_data <- dendro_data(bio_ch_ward, "rectangle")
# 
# p_ward_rl <- ggplot() +
#   geom_segment(data=segment(bio_ch_ward_data),
#                aes(x=x, y=y, xend=xend, yend=yend)) +
#   geom_text(data=label(bio_ch_ward_data),
#             aes(x=x, y=y, label=label, hjust=1)) +
#   coord_flip() +
#   theme_void() +
#   theme(axis.text.x=element_text()) 
# 
# 
# 
# 
# #segments that join the clusters
# comp_ward_link <- tibble(y1=1:nrow(bio_ch_comp_data$labels),
#                          y2=match(bio_ch_comp_data$labels$label, bio_ch_ward_data$labels$label))
# 
# p_m <- ggplot() +
#   geom_segment(data=ward_comp_link, aes(x=0, y=y1, xend=1, yend=y2), 
#                color="grey") +
#   theme_void()
# 
# plot_grid(p_comp_lr, p_m, p_ward_rl, ncol=3, rel_widths=c(1, 0.5, 1))


### Multiscale bootstrap resampling
#assess uncertainty/robustness of a classification (dendrogram) via bootstrap resampling

#### Compute p-values for all clusters (edges) of the dendrogram
bio_norm <- decostand(bio, "normalize") #normalize data (no distances)

bioch_pv <- t(bio_norm) %>% #pvclust takes transposed data
  pvclust(method.hclust="ward.D2", #Ward clustering
          method.dist="euc",
          parallel=TRUE)

#### Plot dendrogram with p-values
plot(bioch_pv)

#### Highlight clusters with high AU p-values
pvrect(bioch_pv, alpha=0.95, pv="au") #au p-values are more accurate than bp p-values
lines(bioch_pv)
pvrect(bioch_pv, alpha=0.91, border=4)
#red boxes + underlined = significant clusters, whereas blue boxes = less significant clusters

#interpretation: most sites are members of significant clusters (p<.05). However, there is one
  #blue box (without underlining) and sites 7 and 21. This supports robust clustering via
  #Ward


### Average silhouette widths
#average measure of degree of membership of an object to its cluster. It uses avg dissimilarity
  #between this object and all objects of the cluster to which it belongs compared to the next
  #closest cluster; range from -1 to 1, with larger values representing stronger clustering

#### Choose and rename the dendrogram ("hclust" object)
hc <- bio_ch_ward

#### Plot average silhouette widths using Ward clustering for all partitions except for trivial
  #partitions 
#book approach
Si <- numeric(nrow(bio))
for(k in 2:(nrow(bio) - 1)) {
  sil <- silhouette(cutree(hc, k=k), bio_ch)
  Si[k] <- summary(sil)$avg.width
}
k.best <- which.max(Si)
plot(
  1:nrow(bio),
  Si,
  type="h",
  main="Silhouette-optimal number of clusters",
  xlab="k (number of clusters)",
  ylab="Average silhouette width"
)
axis(
  1,
  k.best,
  paste("optimum", k.best, sep="\n"),
  col="red",
  font=2,
  col.axis="red"
)
points(k.best,
       max(Si),
       col="red",
       cex=1.5)

#purrr approach
ks <- 2:(nrow(bio)-1)
Si <- rep(0, nrow(bio))

df_sil_width <- ks %>%
  purrr::map_dbl(function(x) {
    silhouette(cutree(hc, k=x), bio_ch) %>%
      summary() %>%
      .[["avg.width"]]
  }) %>%
  tibble(sil_width=.) %>%
  bind_cols(k=2:22, .) %>%
  mutate(k_best=sil_width==max(sil_width))

df_sil_width %>%
  ggplot() +
  geom_linerange(aes(x=k, ymin=0, ymax=sil_width)) +
  geom_point(data=. %>% filter(k_best),
             aes(x=k, y=sil_width),
             color="red",
             size=3) +
  geom_text(data=. %>% filter(k_best),
            aes(x=k, y=sil_width, 
                label=paste0("optimum \n k=",k)),
            color="red",
            nudge_x=2) +
  labs(x="k (number of clusters)",
       y="Average silhouette width",
       title="Silhouette-optimal number of clusters") +
  theme_bw(base_size=13) +
  theme(plot.title=element_text(face="bold", hjust=0.5))





     

