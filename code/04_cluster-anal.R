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








     

