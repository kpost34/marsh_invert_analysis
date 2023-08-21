#Created by Keith Post on 8/13/23
#Exploratory Data Analysis

# Load Packages & DFs===============================================================================
pacman::p_load(here,tidyverse,janitor,sf, googleVis, cowplot, vegan, GGally, car)

source(here("code","01_data-setup.R"))
#once it's determined which is the best data format, then I'll save to RDS and read that in


# Taxonomic Data====================================================================================
## Basic exploration----------
bio <- df_2010_bio_wide[,-1]

head(bio)
tail(bio)
nrow(bio) #23 sites
ncol(bio) #30 taxa
dim(bio)

colnames(bio) #site + taxon names
summary(bio,-1)


## Overall distribution of abundances----------
#range of ns over all taxa
range(bio) #0-1408

#min-max for each taxon
apply(bio, 2, range) #0-1--0-1408

#barplot of distribution
#create species classes: 0: 0, 1:1-10, 2:11-100, 3:101-1000, and 4:1001+
bio %>%
  mutate(across(everything(),
                ~case_when(
                  .x==0                ~ 0,
                  between(.x,1,10)     ~ 1,
                  between(.x,11,100)   ~ 2,
                  between(.x,101,1000) ~ 3,
                  .x > 1000            ~ 4))) -> df_bio_n_class

df_bio_n_class %>%
  unlist() %>%
  tabyl() %>%
  rename(class=".") %>%
  mutate(class=as.character(class)) -> df_bio_n_class_freq

#graph classes
fills <- gray(5:1/6) %>%
  set_names(as.character(0:4))

ggplot(df_bio_n_class_freq) +
  geom_col(aes(x=class,y=n,fill=class)) +
  scale_fill_manual(values=fills,guide="none") +
  labs(x="Abundance class",
       y="Frequency") +
  theme_classic()

#number of absences
sum(bio==0) #279

#proportion of absences
sum(bio==0)/(nrow(bio) * ncol(bio)) #0.404


# Spatial Data======================================================================================
spa <- df_2010_loc_wide %>%
  select(site, longitude="lon",latitude="lat")

## Map of the locations of sites
#assume datum WGS84 (default setting for device)

### Create sf object
sf_spa <- st_as_sf(x=spa,
                   coords=c("longitude","latitude"),
                   crs=4326)

### Plot sites onto lat-long coordinates
ggplot(sf_spa) +
  # geom_sf() + 
  geom_sf_label(aes(label=site),color="red") +
  labs(title="Site Locations",
       x="Longitude",
       y="Latitude") +
  theme_bw(base_size=13) +
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        panel.grid=element_blank())


### Plot onto Google Maps
mymap1 <- spa %>%
  mutate(latlong2 = paste(latitude,longitude,sep=":"),
         .keep="unused") %>% 
  gvisMap(locationvar = "latlong2",
          tipvar = "site",
          options = list(showTip = TRUE)
  )

plot(mymap1)


### Plot data onto maps
#population classes by species
df_bio_n_class %>% 
  select(Chlo, Ceci, Ot_sa, Coll) %>%
  mutate(site=row_number()) %>%
  pivot_longer(cols=!site,
               names_to="taxon",
               values_to="n_class") %>%
  mutate(n_class=as.integer(n_class))-> df_bio_n_class_sub

sf_spa %>%
  left_join(df_bio_n_class_sub) %>%
  left_join(df_bio_cwalk) %>%
  ggplot() +
  geom_sf(aes(size=n_class),shape=1,color="brown") +
  facet_wrap(~taxon_name) +
  labs(x="Longitude",
       y="Latitude") +
  theme_bw(base_size=14) +
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        panel.grid=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(face="bold"),
        legend.position="none")


# Non-spatial Data==================================================================================
## Plot taxonomic occurrences and relative frequencies
### Generate sorted vector of taxonomic presence data
bio_pres <- apply(bio > 0, 2, sum) %>%
  sort()

### Compute rounded percentage frequencies
bio_relf <- (100 * bio_pres/nrow(bio)) %>%
  {round(.,1)} 

### Produce plots
p_freq_occur <- bio_pres %>%
  tibble(n=.) %>%
  ggplot() +
  geom_histogram(aes(n),bins=5,binwidth=5,fill="bisque",color="black",boundary=0) +
  scale_y_continuous(breaks=1:10) +
  labs(x="Number of occurrences",
       y="Number of taxonomic groups",
       title="Taxa Occurrences") +
  theme_classic(base_size=14) +
  theme(plot.title=element_text(hjust=0.5, face="bold"))
#understanding: last bar = there's 7 tax that were found in 20-25 sites


p_freq_rel <- bio_relf %>%
  tibble(freq = .) %>%
  ggplot() +
  geom_histogram(aes(freq),bins=7,fill="bisque",color="black",boundary=0) +
  labs(x="Frequency of occurrences (%)",
       y="Number of taxonomic groups",
       title="Taxa Relative Frequencies") +
  theme_classic(base_size=14) +
  theme(plot.title=element_text(hjust=0.5, face="bold"))
#understanding: number of taxa that occur in #% of sites (i.e., (sites present/total sites)*100)

p_freq <- plot_grid(p_freq_occur, p_freq_rel)


# Data Transformations=============================================================================
## Transform to presence-absence data--------------
bio_pa <- decostand(bio, method="pa")



## Standardize by columns (taxa)------------------
### Scale by taxon max
#divide abundances by max value of each taxon
bio_scal <- decostand(bio, "max")

#display max in each transformed col
apply(bio_scal, 2, max) #all 1s, as expected

### Scale by taxon totals
bio_relsp <- decostand(bio, "total", MARGIN = 2)

#display sum by column
colSums(bio_relsp) #again, all 1s as expected


## Standardize by rows (sites)--------------------
### Scale by site totals
#divide abundances by total individuals per site
bio_rel <- decostand(bio, "total")

#check scaling
rowSums(bio_rel) #yep, all 1s

### Chord transformation
#length (norm) of 1 to each row vector
bio_norm <- decostand(bio, "normalize") #note: default MARGIN = 1

#verify transformation
vec_norm <- function(x) sqrt(sum(x^2)) #computes norm
apply(bio_norm, 1, vec_norm) #apply function to normalized data; all 1s


### Hellinger transformation
#compute square root of relative abundances per site
bio_hel <- decostand(bio, "hellinger")

#verify transformation
apply(bio_hel, 1, vec_norm) #all 1s


## Double standardization by columns and rows------
### Chi-square transformation
bio_chi <- decostand(bio, "chi.square")


### Wisconsin standardization
#range abundances by taxon maxima then site totals
bio_wis <- wisconsin(bio)


## Boxplots of a common taxon-----------------------
bio %>%
  mutate(site=row_number()) %>%
  pivot_longer(cols=!site, names_to = "taxon", values_to = "n") %>%
  left_join(df_bio_cwalk) %>%
  relocate(taxon_name, .after="taxon") %>%
  mutate(across(n))

#create dfs
df_chlo_simple <- bio["Chlo"] %>%
  mutate(sqrt=sqrt(Chlo),
         log=log1p(Chlo)) %>%
  rename(raw ="Chlo") %>%
  pivot_longer(cols=everything(),
              names_to="transformation",
              values_to="values") %>%
  mutate(transformation=fct_relevel(transformation,c("raw","sqrt","log")))

df_chlo_stanspp <- tibble(max=bio_scal$Chlo,
                          total=bio_relsp$Chlo) %>%
  pivot_longer(cols=everything(),
              names_to="transformation",
              values_to="values")
  
df_chlo_stansite <- tibble(Hellinger=bio_hel$Chlo,
                           total_sitetot=bio_rel$Chlo,
                           norm=bio_norm$Chlo) %>%
  pivot_longer(cols=everything(),
          names_to="transformation",
          values_to="values")
  
df_chlo_doubstan <- tibble(`Chi-square`=bio_chi$Chlo,
                           Wisconsin=bio_wis$Chlo) %>%
  pivot_longer(cols=everything(),
          names_to="transformation",
          values_to="values")
  

#develop function
make_transformed_boxplot <- function(dat, fill_color, plot_title) {
  dat %>%
    ggplot(aes(x=transformation,y=values)) +
    #adds horizontal line at end of whisker (and when before geom_boxplot, the vertical line in each
      #box is not displayed)
    stat_boxplot(geom="errorbar",width=0.5) +
    geom_boxplot(fill=fill_color,fatten=3) +
    labs(x="",
         y="",
         title=plot_title) +
    theme_bw(base_size=13) +
    theme(axis.text.x=element_text(size=12),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold",
                                  hjust=0.5))
}

#test function
make_transformed_boxplot(df_chlo_simple,"bisque","Simple transformation")

#create list for pmap
dat_list <- list(df_chlo_simple, df_chlo_stanspp, df_chlo_stansite, df_chlo_doubstan)
col_vec <- c("bisque","lightgreen","lightblue","orange")
title_vec <- c("Simple transformations","Standardizations by species",
               "Standardizations by sites","Double standardizations")

#iterate function over DFs
purrr::pmap(list(dat_list,col_vec,title_vec),make_transformed_boxplot) %>%
  set_names(title_vec) %>%
  plot_grid(plotlist=.)


# Environmental Data================================================================================
env <- df_2010_env_wide[,-1]

# Basic exploration-----------------------------------------
env                       # Display the whole data frame in the 
                          # console
                          # Not recommended for large datasets!
env[1:5, 1:10]            # Display only 5 lines and 10 columns
head(env)                 # Display only the first 6 lines
tail(env)                 # Display only the last 6 rows
nrow(env)                 # Number of rows (sites)
ncol(env)                 # Number of columns (species)
dim(env)                  # Dimensions of the data frame (rows, 
                          # columns)
colnames(env)             # Column labels (descriptors = species)
names(env)                # Alternative
rownames(env)             # Row labels (objects = sites)
summary(env)              # Descriptive statistics for columns


# Bubble Maps of Some Environmental Variables---------------
#isolate env variables & names
env_vars <- c("soc","swc", "ln", "lint")
env_names <- c("Soil Organic Content", "Soil Water Content", "Leaf Nitrogen", "Light Intercepted")

#set fills
env_fills <- c("Soil Organic Content"="red",
               "Soil Water Content"="blue",
               "Leaf Nitrogen"="brown",
               "Light Intercepted"="green3")

#create crosswalk
df_env_crosswalk <- tibble(
  env_code = env_vars,
  env_name = env_names
)

#create long DF
env_sub_long <- env[,env_vars] %>%
  rownames_to_column(var="site") %>%
  mutate(across(!site,~.x/max(.x))) %>%
  pivot_longer(cols=!site,names_to="var",values_to="value") 


#join env data to spatial data
df_xy_env <- sf_spa %>%
  mutate(site=as.character(site)) %>%
  left_join(env_sub_long) %>%
  #pull in names
  left_join(df_env_crosswalk,by=c("var"="env_code")) %>%
  relocate(name="env_name",.after="var") %>%
  mutate(name=fct_relevel(name,df_env_crosswalk$env_name))


#plot data
df_xy_env %>%
  ggplot() +
  geom_sf(aes(fill=name, size=value),shape=21) +
  scale_fill_manual(values=env_fills) +
  labs(x="longitude",
       y="latitude") +
  facet_wrap(~name) +
  theme_bw(base_size = 14) +
  theme(legend.position="none",
        strip.background=element_blank(),
        strip.text=element_text(face="bold"),
        panel.grid=element_blank())



## Scatter plots for all pairs of environmental variables
upperFn <- function(data, mapping, method = "loess", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "black") +
    geom_smooth(method = method, color = "red", se=FALSE, ...)
  p
}

lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "black") +
    geom_smooth(method = method, color = "red", se=FALSE, ...)
  p
}


ggpairs(env,
        upper = list(continuous=wrap(upperFn,method="loess")),
        diag = list(continuous=wrap("barDiag",bins=8,fill="aquamarine2",color="black")),
        lower = list(continuous=wrap(lowerFn,method="lm")),
        title = "Bivariate Plots with Histograms and Smooth Curves") +
  theme(plot.title = element_text(hjust=0.5, face="bold"))


## Simple transformation of an environmental variable
### First example:logit transformation (because variable is measured as %) to improve normality
env_theme <- theme_bw(base_size=12) +
  theme(
    panel.border=element_blank(),
    axis.line=element_line(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold",
                                  hjust=0.5))
p1_env <- env %>%
  ggplot() +
  geom_histogram(aes(x=swc),boundary=0,breaks=seq(0,100,20),binwidth=10,fill="bisque",color="black") +
  labs(x="env$swc",
       y="Frequency",
       title="Histogram of env$swc") +
  scale_x_continuous(breaks=seq(0,100,20)) +
  scale_y_continuous(breaks=seq(0,8,2)) + 
  env_theme

p2_env <- env %>%
  mutate(swc=car::logit(swc)) %>%
  ggplot() +
  geom_histogram(aes(x=swc),boundary=0,breaks=seq(-2,2,.5),binwidth=1,fill="lightgreen",color="black") +
  scale_x_continuous(breaks=seq(-2,2,.5)) +
  scale_y_continuous(breaks=seq(0,10,2)) +
  labs(x="log(env$lint)",
       y="Frequency",
       title="Histogram of logit(env$swc)") +
  env_theme

p3_env <- env %>%
  ggplot(aes(y=swc)) +
  stat_boxplot(geom="errorbar",width=0.25) +
  geom_boxplot(fill="bisque",fatten=2,outlier.shape=1) +
  scale_x_continuous(breaks=NULL) +
  labs(y="env$lint",
       title="Boxplot of env$swc") +
  env_theme

p4_env <- env %>%
  mutate(swc=car::logit(swc)) %>%
  ggplot(aes(y=swc)) +
  stat_boxplot(geom="errorbar",width=0.5) +
  geom_boxplot(fill="light green",fatten=2,outlier.shape=1) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=-2:4) +
  labs(y="log(env$swc)",
       title="Boxplot of logit(env$swc)") +
  env_theme
  
plot_grid(p1_env,p2_env,p3_env,p4_env,ncol=2)

#test normality statistically
shapiro.test(env$swc)
shapiro.test(car::logit(env$swc))



### Another example: non-normal -> normal using log transform
env %>%
  mutate(site=row_number()) %>%
  pivot_longer(!site, names_to="var", values_to="value") %>%
  ggplot() +
  geom_histogram(aes(x=value)) +
  facet_wrap(~var,scales="free_x")


p1_env <- env %>%
  ggplot() +
  geom_histogram(aes(x=pbc),boundary=0,breaks=seq(0,2000,200),binwidth=10,fill="bisque",color="black") +
  labs(x="env$swc",
       y="Frequency",
       title="Histogram of env$pbc") +
  scale_x_continuous(breaks=seq(0,2000,400)) +
  scale_y_continuous(breaks=seq(0,15,3)) +
  env_theme

p2_env <- env %>%
  mutate(pbc=log(pbc)) %>%
  ggplot() +
  geom_histogram(aes(x=pbc),boundary=0,breaks=0:10,binwidth=1,fill="lightgreen",color="black") +
  scale_x_continuous(breaks=0:10) +
  scale_y_continuous(breaks=seq(0,8,2)) +
  labs(x="log(env$pbc)",
       y="Frequency",
       title="Histogram of log(env$pbc)") +
  env_theme

p3_env <- env %>%
  ggplot(aes(y=pbc)) +
  stat_boxplot(geom="errorbar",width=0.25) +
  geom_boxplot(fill="bisque",fatten=2,outlier.shape=1) +
  scale_x_continuous(breaks=NULL) +
  labs(y="env$pbc",
       title="Boxplot of env$pbc") +
  env_theme

p4_env <- env %>%
  mutate(pbc=log(pbc)) %>%
  ggplot(aes(y=pbc)) +
  stat_boxplot(geom="errorbar",width=0.5) +
  geom_boxplot(fill="light green",fatten=2,outlier.shape=1) +
  scale_x_continuous(breaks=NULL) +
  # scale_y_continuous(breaks=-2:4) +
  labs(y="log(env$swc)",
       title="Boxplot of log(env$pbc)") +
  env_theme
  
plot_grid(p1_env,p2_env,p3_env,p4_env,ncol=2)


#test normality statistically
shapiro.test(env$pbc) #non-normal
shapiro.test(log(env$pbc)) #normal




## Standardization of all environmental variables 

# Center and scale = standardize the variables (z-scores)
env.z <- decostand(env, "standardize")
apply(env.z, 2, mean)	# means = 0
apply(env.z, 2, sd)	# standard deviations = 1

# Same standardization using the scale() function (which returns a matrix)
env.z <- as.data.frame(scale(env))






