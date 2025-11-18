00-import
================
Compiled at 2025-11-18 10:36:17 UTC

``` r
here::i_am(paste0(params$name, ".Rmd"), uuid = "c614387e-b8b5-49bd-a8ab-99ba8e026b90")
```

The purpose of this document is read in the data that we generated from
our experiments, as well as some pre-processed information about our
chemical library

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.2     ✔ tibble    3.3.0
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.4     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(readxl)
```

``` r
# create or *empty* the target directory, used to write this file's data: 
projthis::proj_create_dir_target(params$name, clean = TRUE)

# function to get path to target directory: path_target("sample.csv")
path_target <- projthis::proj_path_target(params$name)

# function to get path to previous data: path_source("00-import", "sample.csv")
path_source <- projthis::proj_path_source(params$name)
```

## Read SAMI data.

Here we read in the data first we prepare some functions for this

``` r
write_family <- function(fam_id, complete_df, out_directory, suffix){
  
  # Check if directory exists
  if(!dir.exists(out_directory)){
    
    # If it does not exist, create it
    dir.create(out_directory, recursive = TRUE)
    
  }
  
  
  # Gather relevant data
  family_df <- complete_df %>% 
    filter(Family == fam_id) %>% 
    
    # Order timepoint
    mutate(time_posit = strptime(Time, format = "%m/%d/%Y %H:%M:%S", tz = "CET")) %>% 
    arrange(time_posit)
  
  
  # Gather Plate ID
  plate_id <- family_df %>% 
    select(Plate) %>% 
    unlist() %>% 
    unique()
  
  # Get rid of uneeded data
  family_df <- family_df %>% 
    select(-c(Family, Plate, time_posit))
  
  # Write file
  filename = paste0(plate_id, "_", suffix, ".tsv.gz")
  
  complete.outfile <- file.path(out_directory, filename)
  write_tsv(family_df, complete.outfile)
  
  message(paste("Wrote", filename, "to", out_directory))
  stopifnot(file.exists(complete.outfile))
  
}

process_logfile <- function(file_path, output_dir, timezone="CET", invert=T, Lux=T, batchname = NULL){
  
  # Gather batchname and create subdirectory
  if(is.null(batchname)){
    batch_name <- str_split(file_path, "/")[[1]][4]
  } else{
    batch_name <- batchname %>% 
      filter(`SAMI file` == basename(file_path)) %>%
      pull(batch) %>% 
      unique()
  }
  
  if(is.na(batch_name)){
    print(basename(file.path))
  }
  
  # Read the raw logfile
  dataset <- read.table(file_path, header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE, skip = 4)
  
  # Remove empty columns and remove duplicated rows
  dataset <- dataset %>% 
    select(where(function(x) any(!is.na(x)))) %>% 
    
    distinct()
  
  # Assign new column names depending on inversion or no inversion
  if(invert){
    well_coords <- rev(paste0(rep(LETTERS[1:16], each=24),rep(seq(1,24), 16)))
  } else{
    
    well_coords <- paste0(rep(LETTERS[1:16], each=24),rep(seq(1,24), 16))
    
  }
  
  header <- c("Family","Plate","Time", well_coords)
  stopifnot(length(header) == ncol(dataset))
  colnames(dataset) <- header
  
  
  # Separate Lux and OD data
  if(Lux){
    od_data <- dataset[seq(1, nrow(dataset), by=2),]
    lux_data <- dataset[seq(2, nrow(dataset), by=2),]
    
    stopifnot(nrow(od_data) == nrow(lux_data))
    
    # Write separate files for od and lux per family
    sapply(unique(od_data$Family), write_family, complete_df=od_data, out_directory=file.path(output_dir, "od_data", batch_name), suffix="OD")
    sapply(unique(lux_data$Family), write_family, complete_df=lux_data, out_directory=file.path(output_dir, "lux_data", batch_name), suffix="LUX")
    
    
    #return(list("od_data" = od_data, "lux_data" = lux_data))
  }else{
    
    return(dataset)
  }

}
```

## Copy Plate Database and Experimental Metadata.

``` r
dir.create(path_target("Salmonella/od_data"), recursive = TRUE)
dir.create(path_target("Salmonella/lux_data"), recursive = TRUE)
```

Copy Plate Database

``` r
plate_db <- read_excel("../RAW_DATA/Salmonella/PlateDatabase_OCT2023.xlsx")
write_tsv(file = path_target("Salmonella/PlateDatabase.tsv.gz"), x=plate_db)
```

Copy library map

``` r
libmap_df <- read_tsv("../RAW_DATA/Salmonella/LibMap.txt", show_col_types = FALSE)
write_tsv(x=libmap_df, file = path_target("Salmonella/LibMap.tsv.gz"))
```

## Process SAMI logfiles

Next we can read in all of the SAMI files for each batch

Iterate over sami files

``` r
sami_files <- list.files("../RAW_DATA/Salmonella/", pattern = ".+log", recursive = TRUE)
sami_paths <- file.path("../RAW_DATA/Salmonella", sami_files)

process.status <- sapply(sami_paths, process_logfile, output_dir = path_target("Salmonella"))
```

## Chemical information.

``` r
dir.create(path_target("chemical_library"))
```

``` r
file.copy("../RAW_DATA/ChemicalLibrary/chemical_information_final.tsv.gz",
          path_target("chemical_library/chemical_library.tsv.gz"))
```

    ## [1] TRUE

## External information.

#### Maier 2018

``` r
dir.create(path_target("Maier_2018"))
```

``` r
MAIER_INPUT <- "../RAW_DATA/External/Maier_2018/"

maier_files <- file.path(MAIER_INPUT, list.files(MAIER_INPUT))

file.copy(maier_files, path_target("Maier_2018"))
```

    ## [1] TRUE TRUE TRUE TRUE TRUE

## Files written

These files have been written to the target directory, `data/00-import`:

``` r
projthis::proj_dir_info(path_target())
```

    ## # A tibble: 3 × 4
    ##   path             type             size modification_time  
    ##   <fs::path>       <fct>     <fs::bytes> <dttm>             
    ## 1 Maier_2018       directory         224 2025-11-18 10:37:05
    ## 2 Salmonella       directory         192 2025-11-18 10:36:19
    ## 3 chemical_library directory          96 2025-11-18 10:37:05
