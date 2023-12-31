# Created by Keith Post on 9/1/23
# Cluster Analysis


# Load Packages & DFs===============================================================================
pacman::p_load(here, tidyverse, ade4, adespatial, vegan, gclus, cluster, pvclust, RColorBrewer,
               labdsv, rioja, indicspecies, dendextend, vegclust, 
               colorspace, agricolae, picante, ggdendro, cowplot, sf)

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

#interpretation: robust clusters: 1-3-5-6, 2-11, and 4-8-12


### Multiscale bootstrap resampling
#assess uncertainty/robustness of a classification (dendrogram) via bootstrap resampling

#use UPGMA (had highest cophenetic clustering, cophenetic distance-chord distance plot closest 
  #to diagonal, and smallest Gower's distance)

#### Compute p-values for all clusters (edges) of the dendrogram
bio_norm <- decostand(bio, "normalize") #normalize data (no distances)

#### Compute p-values for all clusters (edges) of the dendrogram
bioch_pv_upgma <- t(bio_norm) %>% #pvclust takes transposed data
  pvclust(method.hclust="average", #UPGMA clustering
          method.dist="euc",
          parallel=TRUE)

#### Plot dendrogram with p-values
plot(bioch_pv_upgma)

#### Highlight clusters with high AU p-values
pvrect(bioch_pv_upgma, alpha=0.95, pv="au") #au p-values are more accurate than bp p-values
lines(bioch_pv_upgma)
pvrect(bioch_pv_upgma, alpha=0.9, border=4) #note: lowered alpha from 0.91 to 0.9
#red boxes + underlined = significant clusters, whereas blue boxes = less significant clusters

#interpretation: this shows 4 clusters (3 clusters of 2+ sites and 1 cluster of 1). Of the three
  #multi-site clusters, two are significant and one is marginally significant. However, note that
  #in comparison to the five-cluster model (from the fusion plots), the two 'subclusters' in the
  #right-most clusters have probabilities (?) of 100 and 98. Note that sites 17 and 18 are adjacent
  #on the map: 17 is in SE NC and 18 is in the peninsula of MD



### Average silhouette widths
#average measure of degree of membership of an object to its cluster. It uses avg dissimilarity
  #between this object and all objects of the cluster to which it belongs compared to the next
  #closest cluster; range from -1 to 1, with larger values representing stronger clustering

#### Choose and rename the dendrogram ("hclust" object)
hc <- bio_ch_UPGMA

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
  bind_cols(k=2:22, .) %>% #ks?
  mutate(k_best=sil_width==max(sil_width))

#hard-coded
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

#by fn
optimize_k_line_graph(df_sil_width, y=sil_width, y_lab="Average silhouette width", 
                      main="Silhouette-optimal number of clusters")

#interpretation: Silhouette plots indicates that k = 5 is the optimal number of clusters. This is
  #supported by the fusion plot (and partially supported by the bootstrapping)


#### Comparison between dissimilarity matrix and binary matrices representing partitions
#compares original dissimilarity matrix with binary matrices of the dendrogram cut at different
  #levels yielding different group memberships. Goal = choose highest correlation


##### Optimal number of clusters according to matrix correlation statistic (Pearson)

grpdist <- function(X) {
  gr <- as.data.frame(as.factor(X))
  distgr <- daisy(gr, "gower")
  distgr
}

#book approach
kt <- data.frame(k=1:nrow(bio), r=0)

for(i in 2:(nrow(bio)-1)) {
  gr <- cutree(hc, i)
  distgr <- grpdist(gr)
  mt <- cor(bio_ch, distgr, method="pearson")
  kt[i, 2] <- mt
}

k.best <- which.max(kt$r)

plot(
  kt$k,
  kt$r,
  type="h",
  main="Matrix correlation-optimal number of clusters",
  xlab="k (number of clusters)",
  ylab="Pearson's correlation"
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
       max(kt$r),
       pch=16,
       col="red",
       cex=1.5)


#purrr approach
df_kt <- 2:22 %>%
  purrr::map_dbl(function(x) {
    cutree(hc, x) %>%
      grpdist() %>%
      cor(bio_ch, method="pearson")
  }) %>%
  c(0, ., 0) %>%
  tibble(k=1:nrow(bio),
         r=.) %>%
  mutate(k_best=r==max(r))

optimize_k_line_graph(df_kt, r, y_lab="Pearson's correlation", 
                      main="Matrix correlation-optimal number of clusters")

#interpretation: 4-6 clusters would yield a high correlation with an optimum at 5 clusters


### Species fidelity analysis
#clusters retained based on presence of "indicator" species, which are species that are more
  #frequent and abundant in sites that make up a cluster
#thus, if clusters are defined in this way, they would maximize 1) sum of indicator values and 2)
  #proportion of clusters with significant indicator species

#book approach
IndVal <- numeric(nrow(bio))
ng <- numeric(nrow(bio))

for(k in 2:(nrow(bio) - 1)) {
  iva <- indval(bio, cutree(hc, k=k), numitr=1000)
  gr <- factor(iva$maxcls[iva$pval <= 0.05])
  ng[k] <- length(levels(gr))/k
  iv <- iva$indcls[iva$pval <= 0.05]
  IndVal[k] <- sum(iv)
}

k.best <- which.max(IndVal)

col3b <- col3a <- rep(1, nrow(bio))
col3a[which.max(IndVal)] <- 3
col3b[which.max(ng)] <- 3


par(mfrow=c(1, 2))
plot(
  1:nrow(bio),
  IndVal,
  type="h",
  main="IndVal-optimal number of clusters",
  xlab="k (number of clusters)",
  ylab="IndVal sum",
  col=col3a
)
axis(
  1,
  k.best,
  paste("optimum", k.best, sep="\n"),
  col="red",
  font=2,
  col.axis="red"
)
points(
  which.max(IndVal),
  max(IndVal),
  pch=16,
  col="red",
  cex=1.5
)
plot(
  1:nrow(bio),
  ng,
  type="h",
  xlab="k (number of groups)",
  ylab="Ratio",
  main="Proportion of clusters with significant indicator species",
  col=col3b
)
points(
  which.max(ng),
  max(ng),
  pch=16,
  col="red",
  cex=1.5
)

#interpretation: no clear best. Optimum number using IndVal sum is 5 clusters but k = 5 does not
  #have ng=1 or the max ng.

#purrr approach
ks <- 2:(nrow(bio)-1)
ng <- numeric(nrow(bio))

df_indval <- ks %>%
  purrr::map_df(function(x) {
    iva <- indval(bio, cutree(hc, k=x), numitr=1000) 
    gr <- factor(iva$maxcls[iva$pval <= 0.05])
    ng[x] <- length(levels(gr))/x
    iv <- iva$indcls[iva$pval <= 0.05]
    IndVal[x] <- sum(iv)
    
    tibble(k=x, ng=ng[x], IndVal=IndVal[x])
  }) %>%
  bind_rows(
    tibble(
      k=c(1, 23), ng=c(0, 0), IndVal=c(0, 0)
    )
  ) %>%
  #k_best = k + 1 where ng==1 (or max ng since none = 1)
  mutate(ng_best=ng==max(ng),
         k_best=IndVal==max(IndVal))

#by fn
#IndVal
p1 <- optimize_k_line_graph(df_indval, 
                      y=IndVal,
                      y_lab="IndVal sum", 
                      main="IndVal-optimal number of clusters",
                      col=TRUE)

#ratio
p2 <- optimize_k_line_graph(df_indval, 
                      y=ng,
                      k_best=ng_best,
                      y_lab="Ratio", 
                      main="Proportion of clusters with significant indicator species",
                      col=TRUE)

plot_grid(p1, p2, labels="auto", label_x=0.9, label_y=0.925)
#interpretation: no clear best. Optimum number using IndVal sum is 5 clusters but k = 5 does not
  #have ng=1 or the max ng. The max ng is 1, meaning that all clusters have indicator species, but
  #this occurs when k=2. However, optimum k = 5 for silhouette widths, matrix correlation, and 
  #fusion plot


### Silhouette plot of the final partition
#go with k=5 and determine if group memberships are appropriate

#book approach
k <- 5

bioch_UPGMA_g <- cutree(bio_ch_UPGMA, k=k)
sil <- silhouette(bioch_UPGMA_g, bio_ch)
rownames(sil) <- row.names(bio)

plot(
  sil,
  main="Silhouette plot - Chord - UPGMA",
  cex.names=0.8, 
  col=2:(k + 1),
  nmax=100
)


#tidyverse approach
#need to recalculate (b/c use k=6 above)
k <- 5
bioch_UPGMA_g <- cutree(bio_ch_ward, k=k)

#plot silhouette profile
#with functions
df_k_sil_upgma <- wrangle_sil_prof_data(data=bio, hclust_obj=bioch_UPGMA_g, dis_mat=bio_ch)
plot_sil_profile(df=df_k_sil_upgma, title="Silhouette plot - Chord - UPGMA")
#interpretation: no misclassifications (all positive), clusters 3 and 4 (smallest clusters aside
  #from solo cluster of site 21) have strongest classification


### Final dendrogram with graphical options
#final tree with custom options
#using dendextend
# Reorder clusters
bio_chUPGMA_o <- reorder.hclust(bio_ch_UPGMA, bio_ch)

# Convert the hclust object into a "dendrogram" object
dend <- as.dendrogram(bio_chUPGMA_o)

# Plot the dendrogram with colored branches
dend %>%
  set("branches_k_color", k=k) %>%
  plot()

# Use standard colors for clusters
clusters <- cutree(dend, k)[order.dendrogram(dend)]

dend %>%
  set("branches_k_color", k=k, value=unique(clusters) + 1) %>%
  plot(main="Final Tree: Chord Distances - UPGMA Method", xlab="Site", ylab="Height")



### Spatial plot of the clustering result
#### Create crosswalk
df_2010_loc_wide %>%
  select(site, lat, lon) %>%
  bind_cols(cluster=as.character(bioch_UPGMA_g)) %>%
  ggplot(aes(x=lon, y=lat)) +
  geom_point(aes(shape=cluster, color=cluster), size=4) +
  geom_text(aes(label=site), nudge_x=0.5, nudge_y=0.25) +
  labs(x="Longitude",
       y="Latitude",
       shape="Cluster") +
  theme_bw()
  
#interpretation: clusters tend to be closer together geographically but there's definitely overlap.
  #Cluster 1 sites occur in the Gulf of Mexico. Cluster 2 sites are in the GoM, se Atlatic coast,
  #and as far north as NJ (site 19). Cluster 3 contains two sites that occur near MD/VA and NC.
  #Cluster 4 comprises two sites that are separated significantly by distance, and Cluster 5
  #has only one site (21) that's on the New England coast.
  

### Heat map and ordered community table
#### Align a heat map with the dendrogram
heatmap(
  1-as.matrix(bio_ch), #am using 1 - association matrix to get a dark diagonal
  Rowv=dend,
  symm=TRUE,
  margin=c(3,3)
)
#interpretation: most closely related sites occur in groups 1-3-6-5-7 (cluster 1), 20-23, 
  #18-17 (cluster 3), and a large group in the middle: 19-10-14-13-15-2-11-4-8-12-16-22 (cluster 2)


#### Reorder and display site-by-species data table
#explore clusters' species content directly by reordering data according to group memberships
#use vegemite() from {vegan} which accepts abundances on 0-9 scale

#reorder data by clusters
#values must be from 0-9 --> use normalized values
#convert to 0-5 classes using Hill scale
bio_hill <- coverscale(bio, scale="Hill") %>%
  mutate(across(everything(), as.numeric))

or <- vegemite(bio_hill, bio_chUPGMA_o)

#apply heatmap
heatmap(
  t(bio_hill[rev(or$species)]),
  Rowv=NA,
  Colv=dend,
  col=c("white", brewer.pal(5, "Greens")),
  scale="none",
  margin=c(4, 4),
  ylab="Taxa (wieghted averages of sites)",
  xlab="Sites"
)

#fyi: similar result can be achieved using tabasco...will look into later but it would be helpful
  #to use final dendrogram (and thus clusters) to color the site-species by those groups and use
  #shading to represent abundance

#interpretation: TBD


# Non-hierarchical Clustering=======================================================================
## k-means partitioning-----
#groups are identified as high-density regions in the data, which minimizes total error sums of
  #squares (i.e., TESS, SSE, sum of within-group SS)
#TESS = summation of the sum of squared Euclidean distances among objects in a group divided by
  #the number of objects among k groups

### k-means with random starts
#inappropriate to use raw species abundance data with lots of 0s because it's a linear method
#one recommendation is to use a non-Euclidean matrix (like Bray-Curtis), sqrt transform it,
  #feed it into a PCoA to get the objects in Euclidean space and provide the PCoA axes into
  #k-means partitioning process. Alternatively, chord transformation followed by computation of
  #Euclidean distances preserves these distances among sites.

#### k-mean partitioning of pre-transformed species data (using same number of groups as UPGMA)
bio_kmeans <- kmeans(bio_norm, centers=5, nstart=100)
#centers = clusters, nstart = # of times analysis is repeated; kmeans() utilizes different
  #random initial configurations

#compare with 5-group results from UPGMA
table(bio_kmeans$cluster, bioch_UPGMA_g)
#given that k-means clustering has a randomness component and that the cluster/group #s are 
  #arbitrarily assigned, sites should not be expected to fall along the diagnonal if they are
  #classified the same using the two methods. Instead, they should be consistently assigned to
  #the same groups, even if those groups differ in number. In other words, there should be only
  #one non-zero value in each row/column, and there are here


#### Choose best number of clusters
#instead of running kmeans() multiple times with different values for centers, vegan's wrapper 
  #cascadeKM() can do this for you
#2-10 groups
bio_KM_cascade <- cascadeKM(
  bio_norm,
  inf.gr=2, #smallest number of groups
  sup.gr=10, #largest number of groups
  iter=100, #number of iterations
  criterion="ssi" #criterion for group selection (= "Simple Structure Index")
)

summary(bio_KM_cascade)
plot(bio_KM_cascade, sortg=TRUE) #sortg=TRUE groups objects according to cluster

#interpretation: per this process using SSI, 9 groups is the best

summary(bio_KM_cascade)
bio_KM_cascade$results
#returns the SSE (lower is better) and ssi (higher is better) for each # of groups
#for a given value of k (within cluster #), SSE used to find the optimal grouping pattern, while
  #ssi (or calinski) are good for determining optimal k value


#### Examine group contents
#reorder sites per k-means result
bio_kmeans_g <- bio_kmeans$cluster #assigns cluster identity for each site
bio[order(bio_kmeans_g), ] #sorts/orders the columns

#reorder sites and species using function vegemite()
ord_KM <- vegemite(bio_hill, bio_kmeans_g) #bio_hill is 'normalized' data (integers from 0-9)
bio[ord_KM$sites, ord_KM$species]


### Use of k-means partitioning to optimize an independently obtained classification
#an advantage of non-agglomerative clustering is using prior information (from clustering) to
  #inform new clustering...in other words, with k-means sites can be seeded with an initial
  #assignment then re-assigned later following new information (with agglomerative clustering,
  #this translocation can't happen)

#### k-means with UPGMA taxa means per group (centroids) as starting points
#seed sites with UPGMA clusters (medoids; best clustering method from hierarchical clustering)
  #and modify with k-means partitioning

#mean taxanomic abundances on UPGMA site clusters
groups <- as.factor(bioch_UPGMA_g) #turn UPGMA cluster ids into a fct
bio_means <- matrix(0, ncol(bio), length(levels(groups))) #create empty matrix populated with 0s
row.names(bio_means) <- colnames(bio) #assign taxonomic names to rows of matrix

#assign normalized values to cells in matrix where columns = clusters
for(i in 1:ncol(bio)) {
  bio_means[i, ] <- tapply(bio_norm[, i], bioch_UPGMA_g, mean)
}

#mean taxonomic abundances as starting points
startpoints <- t(bio_means) #transpose matrix
#k-means on starting points
bio_kmeans2 <- kmeans(bio_norm, centers = startpoints)


#alternatively, return to hierarchical clustering, id the most "typical" object in each group
  #(e.g., silhouette plot), and provide these medoids as starting points
startobjects <- bio_norm[c(6, 13, 18, 20, 21), ]
bio_kmeans3 <- kmeans(bio_norm, centers = startobjects)

#compare with 5-group partition from UPGMA clustering
table(bio_kmeans2$cluster, bioch_UPGMA_g) #seeding with all, normalized values
table(bio_kmeans2$cluster, bioch_UPGMA_g) #seeding with most typical objects
#each kmeans comparison to UPGMA is identical--all sites fall along the diagonal

#compare the two kmeans results after seeding with UPGMA
table(bio_kmeans2$cluster, bio_kmeans3$cluster) #unsurprisingly, they are the same


#Silhouette plot of the final partition
bioch_UPGMA_gk <- bio_kmeans2$cluster
k <- 5
sil <- silhouette(bioch_UPGMA_gk, bio_ch)
rownames(sil) <- row.names(bio)
plot(sil,
     main="Silhouette plot - UPGMA & k-means",
     cex.names=0.8,
     col=2:(k + 1))


#tidyverse approach
bioch_UPGMA_gk <- bio_kmeans2$cluster


#plot silhouette profile
#with functions
df_k_sil_upgma_kmeans <- wrangle_sil_prof_data(data=bio, hclust_obj=bioch_UPGMA_gk, dis_mat=bio_ch)
plot_sil_profile(df=df_k_sil_upgma_kmeans, title="Silhouette plot - Chord - UPGMA & k-means")
#interpretation: clustering is identical between UPGMA and UPGMA-seeded and optimized by k-means



#NOTE: clusters not plotted on map as they would be identical to UPGMA-only plot


## Partitioning around medoids (PAM)----------
#k objects (medoids) are identified that are representative of the structure of the data, and k
  #clusters are constructed by assigning objects to the nearest medioid (goal is to minimize
  #sum of dissimilarities between objects and associated medioids); in comparison, k-means
  #clustering minimizes SS Euclidean distances within groups

#cluster::pam() accepts raw data or dissimilaritiy matrices

### PAM computed from chord distance matrix
#choice of number of clusters--obtain avg silhouette widths (asws) for range of k
asw <- numeric(nrow(bio))

for(k in 2:(nrow(bio) - 1))
  asw[k] <- pam(bio_ch, k, diss=TRUE)$silinfo$avg.width

k.best <- which.max(asw)

plot(
  1:nrow(bio),
  asw,
  type="h",
  main="Choice of the number of clusters",
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
       max(asw),
       pch=16,
       col="red",
       cex=1.5)
#result: k=5 is the best PAM solution, where asw = 0.469. Now let's run a PAM with 5 groups


#tidyverse approach
#purrr approach
ks <- 2:(nrow(bio)-1)

df_sil_width_pam <- ks %>%
  purrr::map_dbl(~pam(bio_ch, .x, diss=TRUE)$silinfo$avg.width) %>%
  tibble(k=ks, asw=.) %>%
  mutate(k_best=asw==max(asw))

optimize_k_line_graph(df_sil_width_pam, y=asw, y_lab="Average silhouette width", 
                      main="PAM Silhouette Plot-optimal number of clusters")


#PAM for k=5
bio_ch_pam <- pam(bio_ch, k=5, diss=TRUE)
summary(bio_ch_pam)

bio_ch_pam_g <- bio_ch_pam$clustering
bio_ch_pam$silinfo$widths

#compare with classification from UPGMA clustering and from k-means
table(bio_ch_pam_g, bioch_UPGMA_g) 
table(bio_ch_pam_g, bio_kmeans_g)

#interpretation: 
  #1) clusters 1 & 2 line up perfectly but 3-5 have 0s on diagonal; however, groups
    #3 & 4 (UPGMA) and 4 & 5 (PAM) align [will need to look into further]; 
  #2) nothing aligns on diagonal...however, a few of the values have no shared values on the same
    #row or column


#silhouette profile for k=5 groups, k-means and PAM
#plot() approach
par(mfrow=c(1, 2))
k <- 5
sil <- silhouette(bio_kmeans_g, bio_ch)
rownames(sil) <- row.names(bio)

#k-means silhouette profile
plot(sil,
     main="Silhouette plot - k-means",
     cex.names=0.8,
     col=2:(k + 1))

#PAM silhouette profile
plot(
  silhouette(bio_ch_pam),
  main="Silhouette plot - PAM",
  cex.names=0.8, 
  col=2:(k + 1)
)

#interpretation: 
  #k-means: 5 cluster: one cluster has one sites, two have two sites, and the other two have at
    #least five sites; asw = 0.5, and they range from 0.45-0.79 (excluding single site cluster)
  #PAM: 5 clusters: two clusters with two sites (same as k-means), one cluster with four sites,
    #one cluster with five sites, and one cluster with 10 sites; asw = 0.47, and they range
    #from 0.14 to 0.79
  #the largest differences reside in clusters 2 and 5


#tidyverse approach
#k-means
df_k_sil_kmeans <- wrangle_sil_prof_data(data=bio, hclust_obj=bio_kmeans_g, dis_mat=bio_ch)
p1 <- plot_sil_profile(df=df_k_sil_kmeans, title="Silhouette plot - k-means")

#PAM
df_k_sil_pam <- wrangle_sil_prof_data(data=bio, hclust_obj=bio_ch_pam)
p2 <- plot_sil_profile(df=df_k_sil_pam, title="Silhouette plot - PAM")

plot_grid(p1, p2)
#largest differences are 1) cluster 2 has 13 sites in k-means and only 10 in PAM; 2) k-means has 
  #a cluster with 1 site (#21) but this site is part of cluster 3 which combines site 21 with the
  #3 sites absent from cluster 2 that are present in the k-means equivalent
  #PAM also has an avg sil that's slightly less than k-means

#NOTE TO SELF: feel free to compare more


# Comparison with Environmental Data================================================================



