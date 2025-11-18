01-growth_analysis
================
Compiled at 2025-11-18 17:29:20 UTC

``` r
here::i_am(paste0(params$name, ".Rmd"), uuid = "898af63f-053d-4795-bae4-02a7cb54a3b2")
```

The purpose of this document is to analyze the growth curves

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(readxl)
library(foreach)
```

    ## 
    ## Attaching package: 'foreach'
    ## 
    ## The following objects are masked from 'package:purrr':
    ## 
    ##     accumulate, when

``` r
library(DGrowthR)
```

``` r
# create or *empty* the target directory, used to write this file's data: 
projthis::proj_create_dir_target(params$name, clean = TRUE)

# function to get path to target directory: path_target("sample.csv")
path_target <- projthis::proj_path_target(params$name)

# function to get path to previous data: path_source("00-import", "sample.csv")
path_source <- projthis::proj_path_source(params$name)
```

## Read all experiment IDs

Here we assess all experiment IDs to get an overview of what experiments
have been performed.

``` r
read_expids <- function(file_name, input_directory){
  
  
  # Read file
  file_path <- file.path(input_directory, file_name)
  timeseries_df <- read.table(file_path, sep='\t', header = TRUE) %>% 
    
    pivot_longer(cols = -c(Time), names_to = "well", values_to = "od") %>% 
    
    select(well) %>% 
    distinct()
  
  
  # Add information
  assay.info <- str_split(file_name, "/")
  batch_name <- assay.info[[1]][1]
  
  plate_id <- assay.info[[1]][2]
  plate_id <- str_remove(plate_id, "_OD.tsv.gz")
  
  
  timeseries_df <- timeseries_df %>% 
    
    mutate(file_origin = file_name,
           batch = batch_name,
           plate_id = plate_id) %>% 
    
    separate(plate_id, sep = "_", into = c("promoter", "libplate", "replicate"))
  
  
  return(timeseries_df)
}


# List files
od_files <- list.files(path_source("00-import/Salmonella/od_data"), recursive = TRUE)

# Do not consider validation screen files yet. 
od_files <- od_files[!str_starts(od_files, "validataion_screen/")] 

registerDoSEQ()
complete_metadata <- foreach(fname = od_files, .combine = "rbind") %dopar% read_expids(fname, path_source("00-import/Salmonella/od_data"))
```

``` r
complete_metadata <- complete_metadata %>% 
  unite("experiment_id", c(promoter, libplate, well, replicate), sep="_", remove=FALSE) %>% 
  unite("Plate_id", c(promoter, libplate, replicate), sep="_", remove=FALSE)
```

## Review data presence.

``` r
file_database <- complete_metadata %>%
  select(promoter, libplate, replicate) %>% 
  distinct() %>% 
  mutate(present = "yes") %>% 
  
  group_by(libplate) %>% 
  pivot_wider(id_cols = c(promoter, libplate), names_from = replicate, values_from = present, values_fill = "no") %>% 
  pivot_longer(cols = -c(promoter, libplate), names_to = "replicate", values_to = "present") %>% 
  
  ungroup()


experiment_presence.g <- file_database %>% 
  #filter(!promoter %in% c("Pzor", "Pvip", "PCj0569", "PCj0164c")) %>% 
  ggplot(aes(x=replicate, y=promoter, fill=present)) +
  geom_tile(color="black") +
  facet_wrap(~libplate, ncol=7) +
  
  theme_minimal()

experiment_presence.g
```

![](01-growth_analysis_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggsave(path_target("experiment_presence.png"), plot=experiment_presence.g, dpi=300)
```

    ## Saving 7 x 5 in image

``` r
complete_metadata %>% select(experiment_id) %>% n_distinct()
```

    ## [1] 161280

## Read optical density.

``` r
update_curve_ids <- function(cid_original){
  
  if(cid_original == "Time"){
    return(cid_original)
  }
  
  # Remove batch
  cid_nobatch <- str_extract(cid_original, ".*?:(P.*)", group=1)
  
  # Obtain parts
  cid_split <- str_split(cid_nobatch, "_", simplify = TRUE)
  
  # Create new CID
  cid_new <- paste(cid_split[1,1], # Promoter
                   cid_split[1,2], #Libplate
                   cid_split[1,4], # Well
                   cid_split[1,3], # Replicate
                   sep = "_")
  
  return(cid_new)
  
}

read_individual_plate_srn <- function(plate_file, input_dir, time_vector){
  
  plate_df <- read.delim(file.path(input_dir, plate_file), sep='\t')
  plate_df$Time <- time_vector
  
  
  # Update the well identifiers to include the plate name
  plate_name <- tools::file_path_sans_ext(plate_file, compression = TRUE)
  plate_name <- str_replace(plate_name, "/", ":")
  plate_name <- str_remove(plate_name, "_OD")
  
  plate_df <- plate_df %>% 
    rename_with(~paste0(plate_name, "_", .x, recycle0 = TRUE), -Time)
  
}

read_multiple_plates_srn <- function(plate_files,
                                     target_path, 
                                     common_time_vector){
  
  
  # List files in path_target
  #plate_files <- list.files(target_path, recursive = TRUE)
  #plate_files <- plate_files[!str_starts(plate_files, "validataion_screen")]
  
  #plate_files <- plate_files
  
  
  #message(paste("The following OD data was found:", paste0(plate_files, collapse = ", ")))
  
  # Read the data from the other lists
  # TODO: can be optimized so as to not re-read the first file
  plates_read_list <- lapply(plate_files, read_individual_plate_srn, input_dir=target_path, time_vector=common_time_vector)
  
  # Combine all dataframes into one big dataframe
  complete_od <- reduce(plates_read_list, full_join, by="Time")
  
  # Create a mapping dataframe
  curve_ids <- colnames(complete_od)[colnames(complete_od) != "Time"]
  mapper_df <- data.frame(curve_info = curve_ids) %>% 
    separate(curve_info, into = c("batch", "plate_well"), sep=":") %>% 
    separate(plate_well, into = c("promoter", "libplate", "replicate", "well"), sep="_") %>% 
    unite("curve_id", c(promoter, libplate, well, replicate), sep="_", remove = FALSE) %>% 
    unite("srn_code", c(libplate, well), sep="_", remove=FALSE)
  
  # Update curve ids
  colnames(complete_od) <- sapply(colnames(complete_od), update_curve_ids)
  
    
  # Return OD and metadata df
  
  return(list("od_data" = complete_od,
              "metadata" = mapper_df))
                                     }
```

Here we must consider whether the measurement was made later or not.

``` r
plate_database <- read_tsv(path_source("00-import/Salmonella", "PlateDatabase.tsv.gz"), show_col_types = FALSE)

plates_lastread.2_1 <- plate_database %>% 
  filter(Notes == "Robot problem. Last read time was later",
         Batch == "2_1") %>% 
  pull(Plate_id)

plates_lastread.3_1 <- plate_database %>% 
  filter(Notes == "Robot problem. Last read time was later",
         Batch == "3_1") %>% 
  pull(Plate_id)


plates_lastread <- c(plates_lastread.2_1, plates_lastread.3_1)
```

We can now consider those experiments separate from the rest

``` r
regular_measurement_files <- complete_metadata %>% filter(!Plate_id %in% plates_lastread) %>% pull(file_origin) %>% unique()


salm_common_time_vector <- c(0, seq(from=0.72, length.out=17, by=0.72))
salm_od.regular <- read_multiple_plates_srn(regular_measurement_files,
                                            path_source("00-import/Salmonella/od_data"),
                                            salm_common_time_vector)
```

Late measurement batch 3_1

``` r
irregular_measurement_files <- complete_metadata %>% filter(Plate_id %in% plates_lastread.3_1) %>% pull(file_origin) %>% unique()


salm_common_time_vector <- c(0, seq(from=0.72, length.out=16, by=0.72), 18.77)
salm_od.irregular.3_1 <- read_multiple_plates_srn(irregular_measurement_files,
                                            path_source("00-import/Salmonella/od_data"),
                                            salm_common_time_vector)
```

Late measurement batch 2_1

``` r
irregular_measurement_files <- complete_metadata %>% filter(Plate_id %in% plates_lastread.2_1) %>% pull(file_origin) %>% unique()


salm_common_time_vector <- c(0, seq(from=0.72, length.out=16, by=0.72), 13.62)
salm_od.irregular.2_1 <- read_multiple_plates_srn(irregular_measurement_files,
                                            path_source("00-import/Salmonella/od_data"),
                                            salm_common_time_vector)
```

``` r
saveRDS(salm_od.regular, path_target("salm_regular_list.rds"))
saveRDS(salm_od.irregular.2_1, path_target("salm_irregular_21_list.rds"))
saveRDS(salm_od.irregular.3_1, path_target("salm_irregular_31_list.rds"))
```

## Update metadata

Join metadata files

``` r
salm_od.regular[["metadata"]] <- salm_od.regular[["metadata"]] %>% 
  mutate(notes = "None")


salm_od.irregular.3_1[["metadata"]] <- salm_od.irregular.3_1[["metadata"]] %>% 
  mutate(notes = "robot_error_3_1")


salm_od.irregular.2_1[["metadata"]] <- salm_od.irregular.2_1[["metadata"]] %>% 
  mutate(notes = "robot_error_2_1")


salm.metadata <- bind_rows(salm_od.regular[["metadata"]], salm_od.irregular.2_1[["metadata"]], salm_od.irregular.3_1[["metadata"]])
```

## Combine OD.

``` r
salm_od.irregular.2_1.OD_DATA <- salm_od.irregular.2_1[["od_data"]] %>% 
  mutate(timepoint_n = 1:n()) %>% 
    pivot_longer(cols = -c(Time, timepoint_n), names_to = "curve_id", values_to = "od") %>% 
    rename("timepoint" = "Time")

salm_od.irregular.3_1.OD_DATA <- salm_od.irregular.3_1[["od_data"]] %>% 
  mutate(timepoint_n = 1:n()) %>% 
    pivot_longer(cols = -c(Time, timepoint_n), names_to = "curve_id", values_to = "od") %>% 
    rename("timepoint" = "Time")

salm_od.regular.OD_DATA <- salm_od.regular[["od_data"]] %>% 
   mutate(timepoint_n = 1:n()) %>% 
    pivot_longer(cols = -c(Time, timepoint_n), names_to = "curve_id", values_to = "od") %>% 
    rename("timepoint" = "Time")


salm_od.long <- bind_rows(salm_od.regular.OD_DATA, salm_od.irregular.2_1.OD_DATA, salm_od.irregular.3_1.OD_DATA)
```

## Create DGrowthR object.

``` r
salm.dgobj <- new("DGrowthR", 
                od_data = salm_od.long, 
                metadata = salm.metadata,
                raw_od = salm_od.long,
                verbose=TRUE)
```

## Pre-process samples.

``` r
salm.dgobj.linear <- preprocess_data(salm.dgobj, 
                                     baseline = c(3,4),
                                     log_transform = FALSE,
                                     skip_first_n_timepoints = 2)
```

    ## Ignoring first 2 for analysis.

    ## Using linear growth measurements...

    ## WARNING: Not log-transforming growth mesaurements might result in inaccurate estimation of growth parameters...

    ## Using minimum OD from timepoints 3, 4 as baseline

``` r
salm.dgobj.log <- preprocess_data(salm.dgobj, 
                                     baseline = 4,
                                     log_transform = TRUE,
                                     skip_first_n_timepoints = 3)
```

    ## Ignoring first 3 for analysis.

    ## Log-transforming the growth measurements...

    ## Using OD from timepoint 4 as baseline.

## Write files.

Since these steps already take a long time, we can just set these as a
separate step.

``` r
saveRDS(salm.dgobj, path_target("salm_dgobj.rds"))
saveRDS(salm.dgobj.linear, path_target("salm_dgobj_linear.rds"))
saveRDS(salm.dgobj.log, path_target("salm_dgobj_log.rds"))

saveRDS(salm_od.long, path_target("complete_od.rds"))
saveRDS(salm.metadata, path_target("complete_metadata.rds"))
```

## Files written

These files have been written to the target directory,
`data/01-growth_analysis`:

``` r
projthis::proj_dir_info(path_target())
```

    ## # A tibble: 9 × 4
    ##   path                       type         size modification_time  
    ##   <fs::path>                 <fct> <fs::bytes> <dttm>             
    ## 1 complete_metadata.rds      file      656.13K 2025-11-18 17:33:54
    ## 2 complete_od.rds            file       12.27M 2025-11-18 17:33:54
    ## 3 experiment_presence.png    file      161.69K 2025-11-18 17:29:35
    ## 4 salm_dgobj.rds             file       25.18M 2025-11-18 17:33:39
    ## 5 salm_dgobj_linear.rds      file       24.19M 2025-11-18 17:33:45
    ## 6 salm_dgobj_log.rds         file       31.59M 2025-11-18 17:33:51
    ## 7 salm_irregular_21_list.rds file      193.35K 2025-11-18 17:32:26
    ## 8 salm_irregular_31_list.rds file      194.72K 2025-11-18 17:32:26
    ## 9 salm_regular_list.rds      file        6.14M 2025-11-18 17:32:26
