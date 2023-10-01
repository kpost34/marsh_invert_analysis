#Created by Keith Post on 8/27/23
#Functions

# Functions for Association Measures================================================================
## Plot dissimilarity matrix
plot_assoc_matrix <- function(dist.obj, ordered = FALSE) {
  #use proportions if values are > 1
  maxD <- max(dist.obj)
  
  dist.obj.ord <- if(ordered & maxD > 1) {
    1-(dist.obj/maxD)
  } else if(ordered){
    1-dist.obj
  } else{dist.obj}
  
  #populate site levels depending on ordered arg
  ord_site <- order.single(dist.obj.ord)
  n_obj <- ncol(dist.obj)
  
  site_levels <- if(ordered) {
    ord_site} else{1:n_obj}

  #plot matrix
  dist.obj %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column(var="site_x") %>% 
    pivot_longer(!site_x, names_to="site_y", values_to="dist") %>% 
    {if(maxD>1) mutate(.,dist=dist/maxD) else .} %>%
    mutate(dist = ifelse(site_x==site_y,
                         NA_real_,
                         dist),
           across(!dist,~factor(.x, levels=site_levels)),
           site_y = fct_rev(site_y)) %>% 
    ggplot() +
    geom_tile(aes(x=site_x, y=site_y, fill=dist),
              color = "white") +
    scale_fill_stepsn(breaks = c(.25, .5, .75, 1),
                      colors = c("orchid3", "orchid1", "darkslategray1", "cyan2"),
                      na.value = "white") +
    labs(x="",
         y="") +
    {if(!ordered) ggtitle("Dissimilarity Matrix") else ggtitle("Ordered Dissimilarity Matrix")} +
    theme(legend.position="none",
          plot.title=element_text(face="bold", hjust=0.5))
}



## Plot dissimilarity matrix for taxa
plot_assoc_matrix_taxa <- function(dist.obj, ordered = FALSE) {


  #use proportions if values are > 1
  maxD <- max(dist.obj)
  
  dist.obj.ord <- if(ordered & maxD > 1) {
    1-(dist.obj/maxD)
  } else if(ordered){
    1-dist.obj
  } else{dist.obj}
  
  
  #populate site levels depending on ordered arg
  #isolate taxon names
  taxon_nm <- names(dist.obj)
  
  ord_taxon <- order.single(dist.obj.ord) %>%
    taxon_nm[.]
    
  
  taxon_levels <- if(ordered) {
    ord_taxon
    } else{taxon_nm}
  
  #plot matrix
  dist.obj %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column(var="taxon_x") %>% 
    pivot_longer(!taxon_x, names_to="taxon_y", values_to="dist") %>% 
    {if(maxD>1) mutate(.,dist=dist/maxD) else .} %>%
    mutate(
           dist = ifelse(taxon_x==taxon_y,
                         NA_real_,
                         dist),
           across(!dist,~factor(.x, levels=taxon_levels)),
           taxon_y = fct_rev(taxon_y)) %>% 
    ggplot() +
    geom_tile(aes(x=taxon_x, y=taxon_y, fill=dist),
              color = "white") +
    scale_fill_stepsn(breaks = c(.25, .5, .75, 1),
                      colors = c("orchid3", "orchid1", "darkslategray1", "cyan2"),
                      na.value = "white") +
    labs(x="",
         y="") +
    {if(!ordered) ggtitle("Dissimilarity Matrix") else ggtitle("Ordered Dissimilarity Matrix")} +
    theme(legend.position="none",
          plot.title=element_text(face="bold", hjust=0.5),
          axis.text.x=element_text(angle=90))
}

  
  
## Create wrapper function
plot_dissim_grid <- function(dist.obj, fn = plot_assoc_matrix){
  c(FALSE, TRUE) %>%
    purrr::map(function(x) {
      fn(dist.obj, ordered = x)
    }) %>%
    plot_grid(plotlist=.)
}


# Functions for Cluster Analysis====================================================================
## Grab site label positions and names
extract_site_tree_info <- function(cluster.tree, nudge=.01) {
  segment(cluster.tree) %>%
    filter(xend %in% 1:23) %>%
    group_by(xend) %>%
    filter(yend==min(yend)) %>% 
    ungroup() %>%
    mutate(yend=yend - nudge) %>%
    select(ends_with("end")) %>%
    left_join(
      label(cluster.tree)[,c("x", "label")],
      by=c("xend"="x")
    )
}


## Plot dendrogram
plot_dendro <- function(cluster.obj, title="", nudge=.01) {
  #convert cluster object to dendro object
  cluster_tree <- cluster.obj %>% #hclust object
    as.dendrogram() %>% #converts to dendrogram object
    hang.dendrogram() %>% #makes leaves hang above x-axis
    dendro_data(type="rectangle") #extracts plotting data

  #extract labels
  site_labs <- extract_site_tree_info(cluster_tree, nudge)
  
  #set ymin
  ymin <- min(site_labs$yend) + nudge
  ymin <- if(ymin < 0) {ymin} else{0}
  
  #plot dendrogram
  cluster_tree %>%
    #extracts x-y data from list-object
    segment() %>%
    ggplot() +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data=site_labs,
              aes(x=xend, y=yend, label=label)) +
    coord_cartesian(ylim=c(ymin-0.02, NA)) +
    labs(x="site",
         y="height",
         title=title) +
    theme_void() + #removes axes, axis.titles, and axis.text
    #specify margins as well as center and bold plot title & rotate y-axis title 90o
    theme(plot.title=element_text(hjust=0.5, face="bold", margin=margin(t=10, b=10)),
          axis.line.y=element_line(),
          axis.title.x=element_text(size=14, margin=margin(b=5)),
          axis.title.y=element_text(size=14, angle=90, margin=margin(r=10, l=5)),
          axis.text.x=element_blank(),
          axis.text.y=element_text(margin=margin(r=10)))
}


## Convert distances from wide (matrix) to long form
convert_to_long_df <- function(dist.obj, type=NA) {
  dist.obj %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column(var="site_x") %>% 
    pivot_longer(!site_x, names_to="site_y", 
                 values_to=ifelse(!is.na(type),
                                  paste(type,"dist",sep="_"),
                                  "dist"))
}


## Make line graphs to determine optimal k (e.g., silhouette plots, IndVal)
optimize_k_line_graph <- function(data, y, k_best=k_best, y_lab, main, col=FALSE) {
  p <- data %>%
    ggplot() +
    geom_linerange(aes(x=k, ymin=0, ymax={{y}})) +
    geom_point(data=. %>% filter({{k_best}}),
               aes(x=k, y={{y}}),
               color="red",
               size=3) +
    geom_text(data=. %>% filter({{k_best}}),
              aes(x=k, y={{y}}, 
                  label=paste0("optimum \n k=",k)),
              color="red",
              nudge_x=2) +
    labs(x="k (number of clusters)",
         y=y_lab,
         title=main) +
    theme_bw(base_size=13) +
    theme(plot.title=element_text(face="bold", hjust=0.5))
  
  if(col){
    p +
      geom_linerange(data=. %>% filter({{k_best}}),
                     aes(x=k, ymin=0, ymax={{y}}),
                     color="green")
  } else{p}
}


## Silhouette plots----------
### Wrangle data
wrangle_sil_prof_data <- function(data, hclust_obj, dis_mat=NA) {
  
  n <- nrow(data)
  
  sil <- if(!is.na(dis_mat)) {
    silhouette(hclust_obj, dis_mat)} else{
      silhouette(hclust_obj)
    }
  
  # df_k_sil <- tibble(site=seq_len(n)) %>%
  #   bind_cols(sil) %>%
    
  df_k_sil <- sil %>%
    as.data.frame() %>%
    rownames_to_column(var="site") %>%
    group_by(cluster) %>%
    #'switch' for assigning clusters depending on whether dis_mat is NA or not
    {if(!is.na(dis_mat))
      mutate(., k_size=n(),
             sil_width_avg=mean(sil_width) %>%
               round(., 2),
             text_x=case_when(
               cluster=="1" ~ 0.93*n,
               cluster=="2" ~ 0.5*n,
               cluster=="3" ~ 0.2*n,
               cluster=="4" ~ 0.1*n,
               cluster=="5" ~ NA_real_)
             ) else mutate(., k_size=n(),
             sil_width_avg=mean(sil_width) %>%
               round(., 2),
             text_x=case_when(
               cluster=="1" ~ 0.93*n,
               cluster=="2" ~ 0.6*n,
               cluster=="3" ~ 0.3*n,
               cluster=="4" ~ 0.15*n,
               cluster=="5" ~ 0.07*n)
             )} %>%
    ungroup() %>%
    mutate(sil_width_global_avg=mean(sil_width) %>%
             round(., 2)) %>%
    arrange(desc(cluster), sil_width) %>%
    mutate(cluster=as.factor(cluster),
           site=as.character(site),
           site=fct_inorder(site))
}


### Graph silhouette profile
plot_sil_profile <- function(df, title) {
  
  #text labels
  df_k_sil_lab <- df %>%
    select(cluster, k_size, sil_width_avg, text_x) %>%
    distinct()
  
  #colors (fills)
  k_cols <- c("1" = "red", "2" = "green", "3" = "blue", "4" = "aquamarine", "5" = "purple")
  
  n <- nrow(df)
  
  df %>%
    ggplot() +
    geom_col(aes(x=site, y=sil_width, fill=cluster)) +
    geom_text(data=df_k_sil_lab,
              aes(x=text_x, y=0.925,
                  label=paste0(cluster, ": ", k_size, " | ", sil_width_avg))) +
    annotate("text", x= n, y=0.925,
             label=paste("Avg si:", unique(df[["sil_width_global_avg"]]))) +
    scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
    scale_fill_manual(values=k_cols, guide="none") +
    coord_flip() +
    labs(title=title,
         subtitle=paste0("n = ", n),
         y="Silhouette width si") +
    theme_void() +
    theme(axis.text=element_text(),
          axis.line.x=element_line(),
          axis.title.x=element_text(margin=margin(t=5, b=10)),
          plot.title=element_text(face="bold"))
}








