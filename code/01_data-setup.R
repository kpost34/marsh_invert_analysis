#Created by Keith Post on 8/6/23
#Objectives: 1)import data, 2) clean data, 3) assess missingness, 4) truncate and impute, 
  #5) spin data to different formats

# Load Packages=====================================================================================
pacman::p_load(here, tidyverse, janitor, skimr, mice, naniar)


# Data Import and Initial Cleaning==================================================================
## Set path
data_fp <- here("data", "raw_data", "MLT-GCET-1106a_1_1.CSV")

## Read in data
### Extra units
extra_units <- "none|^d$|^m$|^yyyy$"

### Read data and headers
df_input <- read_csv(data_fp, skip=5, col_names=FALSE, guess_max=50)
header_input <- read_csv(data_fp, skip=2, n_max=2, col_names=FALSE)

### Wrangle header rows
df_header <- header_input %>%
  mutate(across(everything(), make_clean_names),
         across(everything(), ~ifelse(row_number()==2,
                                      str_remove(.x,extra_units) %>%
                                        str_replace("e_m_2_s","ue_m_2_s"),
                                      .x)),
         across(everything(), ~ifelse(nchar(.x[row_number()==2])>0,
                                      paste(.x[row_number()==1],
                                            .x[row_number()==2],
                                            sep="__"),
                                      .x))) %>% 
  .[1,] %>%
  as.character()
  

### Append wrangled header to data, remove oil-exposed samples, and create date field
df_full_clean <- df_input %>%
  set_names(df_header) %>% 
  filter(oil_presence==0) %>%
  mutate(date=make_date(year, month, day),
         .after="longitude__degrees")



# Partition Data====================================================================================
## By sampling period to compare
df_full_clean %>%
  filter(year==2010) -> df_2010

df_full_clean %>%
  filter(year==2009,
         month %in% 4:5) -> df_2009_1

df_full_clean %>%
  filter(year==2009,
         month %in% 7:8) -> df_2009_2

df_list <- list(df_2010,
                df_2009_1,
                df_2009_2) %>%
  set_names(c("subset2010","subset2009_spring","subset2009_summer"))


## By loc/env, non-dvac spp, and d-vac spp
df_env_list <- df_list %>%
  purrr::map(select,1:19)

df_sp_ndvac_list <- df_list %>%
  purrr::map(select,20:41)

df_sp_dvac_list <- df_list %>%
  purrr::map(select,contains("d_vac"))



# Compare Partitions===============================================================================
## Dimensions
df_list %>%
  purrr::map(dim)


## Look at summaries of location and env fields
df_env_list  %>%
  purrr::map(skim)
#missingness:
#2010: 5 (all 1): soil water, soil salinity, par above canopy, par below canopy, light intercepted
#2009 spring: 10 (4 = all missing)
#2009 summer: 5 (4 = all missing, including both par and light intercepted)


## Total missingness
#overall
df_list %>%
  purrr::map(n_miss)
#2010: 61
#2009 spring: 850
#2009 summer: 705

#env vars
df_env_list %>%
  purrr::map(n_miss)
#5, 106, 89

#non-dvac vars
df_sp_ndvac_list %>%
  purrr::map(n_miss)
#56, 480, 440


#dvac vars
df_sp_dvac_list %>%
  purrr::map(n_miss)
#0, 264, 176


## Missing case summaries
#env vars
df_env_list %>%
  purrr::map(miss_case_summary)
#2010: only 2 cases with miss values
#2009 spring: at least 10 cases have at least 4 missing values
#2009 summer: same as 2009 spring

#non-dvac vars
df_sp_ndvac_list %>%
  purrr::map(miss_case_summary)
#2010: at least 10 cases with at least 9% missing values; at least 90% missing
  #values for at least 10 cases for both 2009 datasets

#dvac vars
df_sp_dvac_list %>%
  purrr::map(miss_case_summary)
#2010: none, both 2009: at least 10 cases with at least 25% missingness


## Missing variables
#env vars
df_env_list %>%
  purrr::map(miss_var_table)
#2010: 14 vars complete; 5 vars with 1 missing
#2009 spring: 9 vars complete; 2 vars with 1 missing; 4 vars with 2 missing; 4 vars with 24 missing
#2009 summer: 14 vars complete; 1 var with 1 missing; 4 vars with 22 missing

#non-dvac vars
df_sp_ndvac_list %>%
  purrr::map(miss_var_table)
#2010: 17 vars with 0 or 2 missing; 5 vars with 4-8 missing
#2009 spring: heavy missingness (e.g., 20 vars with 24 missing)
#2009 summer: heavy missingness (e.g., 20 vars with 22 missing)

#dvac vars
df_sp_dvac_list %>%
  purrr::map(miss_var_table)
#2010: all vars complete
#2009 spring: 24 vars with 3 missing; 8 vars with 24 missing
#2009 summer: 24 vars complete; 8 vars with 22 missing


#Conclusion: given all of this, the 2010 sampling period has far less missing data than the
  #two 2009 sampling periods, so this sampling period will be used in the analysis


# Closer Look at 2010 Data==========================================================================
## Summary
df_2010 %>%
  skim()
#as mentioned above, there's missing loc/env data in SWC, soil salinity, PAR above, PAR below,
  #and light intercepted ratio (all 1 per var for 5 total)
#non-dvac sampling: many cases with 2 missing values, two with 4 missing values, one with 5,
  #one with 7, and one with 8
#dvac complete

## Investigate data types more deeply
#env data
df_2010 %>%
  select(location:light_intercepted__ratio) %>% 
  miss_case_summary()
#case 3 has 3 missing values and case 19 has 2 missing values

#non-dvac org data
df_2010 %>%
  select(littoraria_density__count_m_2:unidentified_infauna_count_sediment__count) %>%
  miss_case_table()
#12 complete cases
#9 cases with 2-4 NAs
#2 cases with 14 NAs


#Conclusion: 
#Given that there's only 5 missing values (from 2 cases) for the loc/env data and all dvac
  #sampling data are complete, these fields will be retained. The non-dvac sampling will be 
  #removed (this will leave the same sampling method for the dependent response matrix). Further,
  #given that there are few missing env data in a small sample size, imputation is the best
  #approach.


# Data Imputation===================================================================================
## Generate a subset of the data
df_2010 %>% 
  #loc & env vars
  select(location:year,
         soil_organic_content__percent:light_intercepted__ratio, #omit oil_presence
         contains("d_vac")) %>% #only dvac dependent vars
  remove_constant(quiet=FALSE) %>%
  #bring back year
  mutate(year = year(date), .after="month") -> df_2010_sub

#remove location & date info
df_2010_sub %>%
  select(soil_organic_content__percent:light_intercepted__ratio) -> df_2010_env

df_2010_sub %>%
  select(contains("d_vac")) -> df_2010_bio

## Look at missingness more deeply
gg_miss_upset(df_2010_env)
gg_miss_var(df_2010_env)
gg_miss_case(df_2010_env)
gg_miss_which(df_2010_env)


## Test for whether data are MCAR [only independent variables]
naniar::mcar_test(df_2010_env)
#rationale: missingness lies only in independent variables


## Impute missing values
df_2010_env %>% 
  mice(method="cart") -> df_imp_env

df_2010_full_env <- complete(df_imp_env) %>%
  as_tibble()


## Develop colnames
#loc
loc_nm <- c("loc","lat","lon","date","day","mon","yr")

#env
env_nm <- c("soc__percent","swc__percent","ss__psu","sbit__percent","ln__percent","spc__percent",
            "dpc__percent","sph__cm","pac__ue_m2_s","pbc__ue_m2_s","lint__ratio")
env_abbr <- str_remove(env_nm,"__.+$")

#taxonomic
bio_nm <- names(df_2010_bio) %>%
  str_remove("_count_d_vac__count$") %>%
  ifelse(str_detect(.,"_"),
         paste0(str_sub(.,1,2),str_extract(.,"_[a-z]{2}")),
         .) %>%
  ifelse(!str_detect(.,"_"),
         str_sub(.,1,4),
         .) %>%
  str_to_title()

#develop crosswalk
df_bio_cwalk <- names(df_2010_bio) %>%
  str_remove("_count_d_vac__count$") %>%
  str_to_title() %>%
  tibble(taxon_name=.) %>%
  bind_cols(
    taxon = bio_nm
  )

## Combine DFs
df_2010_wide <- df_2010_full_env %>%
  bind_cols(
    df_2010_sub %>%
      select(location:year),
    .,
    df_2010_bio
  ) %>% 
  set_names(c(loc_nm, env_nm, bio_nm)) %>%
  #add a simple site number
  mutate(site=row_number(),.before="loc")
                                    
  
  

# Data Partitioning=================================================================================
## Split env and bio data
df_2010_wide %>%
  select(site:yr) -> df_2010_loc_wide #site and loc-date

df_2010_wide %>% 
  select(site,soc__percent:lint__ratio) %>%
  set_names(c("site",env_abbr)) -> df_2010_env_wide #site & env only
  
df_2010_wide %>%
  select(site,Tr_uh:last_col()) -> df_2010_bio_wide #site & taxonomic data

df_2010_wide %>%
  set_names(c("site",loc_nm,env_abbr,bio_nm)) -> df_2010_wide_abbr #full DF (with abbrv)


## Pivot to long formats
### env data
df_2010_wide %>%
  select(site,all_of(env_nm)) %>%
  pivot_longer(cols=!site, 
               names_to=c("var", "units"),
               names_sep="__", 
               values_to="meas_value") -> df_2010_env_long #site & env

### bio data
df_2010_bio_wide %>%
  pivot_longer(cols=!site, 
               names_to="taxon",
               values_to="count") -> df_2010_bio_long #site & bio
  
            
### all data
#site, env, and taxonomic data
df_2010_env_long %>%
  left_join(df_2010_bio_long) -> df_2010_long

#site, loc-date, env, and taxonomic data
df_2010_loc_wide %>%
  left_join(df_2010_env_long) %>%
  left_join(df_2010_bio_long) -> df_2010_full_long


  


## Remove all DFs but
rm(list=setdiff(ls(),ls(pattern="wide|long|cwalk")))




