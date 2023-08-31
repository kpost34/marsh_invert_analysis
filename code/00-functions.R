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
    mutate(dist = ifelse(taxon_x==taxon_y,
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













