02-growth_filtering
================
Compiled at 2025-11-18 19:30:08 UTC

``` r
here::i_am(paste0(params$name, ".Rmd"), uuid = "cddd9485-e271-4bd8-ba61-563b2e390f26")
```

The purpose of this document is …

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

## Read DGrowthR object.

``` r
salm.dgobj.linear <- readRDS(path_source("01-growth_analysis", "salm_dgobj_linear.rds"))
```

## Determine flat curves

``` r
max_growth <- salm.dgobj.linear@od_data %>% 
  
  #filter(curve_id %in% pilot_curves) %>% 
  filter(curve_id != "PsroC_lp2_A2_r2", # Potential overgrowth
         timepoint < 13) %>% # To account for late measurements in the july experiments
  
  filter(od == max(od)) %>% 
  pull(od)

minimum_growth_filter <- max_growth * 0.2
```

``` r
flat.curves.df <- salm.dgobj.linear@od_data %>% 
  group_by(curve_id) %>% 
  filter(timepoint_n < 18) %>% # To account for the july experiments. Also, this means at least two timepoints have to pass the threshold
  summarise(flat_curve = if_else(max(od) < minimum_growth_filter, TRUE, FALSE))

flat.curves.df %>% count(flat_curve) %>% knitr::kable()
```

    ## Warning: 'xfun::attr()' is deprecated.
    ## Use 'xfun::attr2()' instead.
    ## See help("Deprecated")
    ## Warning: 'xfun::attr()' is deprecated.
    ## Use 'xfun::attr2()' instead.
    ## See help("Deprecated")

| flat_curve |      n |
|:-----------|-------:|
| FALSE      | 155580 |
| TRUE       |   5700 |

``` r
max_od.df <- salm.dgobj.linear@od_data %>% 
  filter(timepoint_n < 18) %>% 
  group_by(curve_id) %>% 
  summarise(max.od = max(od)) %>% 
  ungroup()

salm.filter.hist <- max_od.df %>% 
  filter(curve_id != "PsroC_lp2_A2_r2") %>% 
  #mutate(max_growth = if_else(is.na(max_growth), 0, max_growth)) %>%
  
  
  ggplot(aes(x=max.od)) +
  geom_histogram(binwidth=0.01) +
  
  annotate("text", x= 0.09, y=6000, label="Failed filter\n5,700", size=3) +
  annotate("text", x= 0.55, y=6000, label="Passed filter\n155,580", size=3) +
  
  geom_vline(xintercept = minimum_growth_filter, color="red") +
  
  theme_bw() +
  
  labs(title = latex2exp::TeX("\\textit{S.} Typhimurium"),
       x = "Maximum OD",
       y = "Count") 

salm.filter.hist
```

![](02-growth_filtering_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
ggsave(path_target("filter_histogram.png"), plot=salm.filter.hist, dpi=300)
```

    ## Saving 7 x 5 in image

## Determining growth inhibitors.

``` r
salm_killers.df <- flat.curves.df %>% 
  separate(curve_id, into=c("promoter", "libplate", "well", "replicate"), remove = FALSE) %>% 
  unite("srn_code", c(libplate, well), remove=FALSE) %>% 
  
  group_by(srn_code) %>% 
  summarise(n_tests = n(),
            n_dead = sum(flat_curve)) %>% 
  
  mutate(prop.killed = n_dead/n_tests) 


salm_killers <- salm_killers.df %>% filter(prop.killed >= 0.5) %>% pull(srn_code)

length(salm_killers)
```

    ## [1] 93

We can begin tracking which experiments should be removed and why.

``` r
filter.metadata <- salm.dgobj.linear@metadata %>% 
  mutate(filter.status = if_else(srn_code %in% salm_killers, "remove", "keep"),
         
         filter.reason = if_else(srn_code %in% salm_killers, "killer.compound", "None"))
```

## Remove low growers.

Now we can remove compounds that have low growth despite not being a
killer compounds

``` r
flat_curves <- flat.curves.df %>% filter(flat_curve) %>% pull(curve_id)
```

``` r
filter.metadata <- filter.metadata %>% 
  mutate(filter.status = case_when(filter.status != "remove" & curve_id %in% flat_curves ~ "remove",
                                   TRUE ~ filter.status),
         
         filter.reason = case_when(filter.reason == "None" &  curve_id %in% flat_curves ~ "low.growth",
                                   TRUE ~ filter.reason))
```

Now we make sure that at two replicates pass.

``` r
nrep_after_filter <- filter.metadata %>% 
  unite("promoter_compound", c(promoter, libplate, well), remove = FALSE) %>% 
  
  filter(filter.status == "keep") %>% 
  count(promoter_compound) %>% 
  
  rename("nrep_after_filter" = "n")

filter.metadata <- filter.metadata %>% 
  unite("promoter_compound", c(promoter, libplate, well), remove = FALSE) %>% 
  left_join(nrep_after_filter, by="promoter_compound") %>% 
  
 mutate(filter.status = case_when(filter.status != "remove" & nrep_after_filter < 2 ~ "remove",
                                   TRUE ~ filter.status),
         
         filter.reason = case_when(filter.reason == "None" & nrep_after_filter < 2 ~ "nrep_after_growth_filters",
                                   TRUE ~ filter.reason))


filter.metadata %>% count(filter.reason)
```

    ##               filter.reason      n
    ## 1                      None 155326
    ## 2           killer.compound   5580
    ## 3                low.growth    268
    ## 4 nrep_after_growth_filters    106

## Removing outlier curves.

Some growth curves have unsual shapes we can remove them.

``` r
odd.curves.p <- salm.dgobj.linear@od_data %>% 
  filter(curve_id %in% c("PgcvB_lp7_D5_r1", "PgcvB_lp7_D5_r2",
                         "PsroC_lp2_A2_r2", "PsroC_lp2_A2_r1")) %>% 
  
  left_join(salm.dgobj.linear@metadata %>% select(curve_id, promoter, srn_code, replicate)) %>% 
  unite("promoter_compound", c(promoter, srn_code)) %>% 
  
  
  ggplot(aes(x=timepoint, y=od, group=curve_id, color=replicate)) +
  geom_line() +
  
  facet_wrap(~promoter_compound) +
  
  theme_bw() +
  
  labs(x="Time (h)",
       y="OD")
```

    ## Joining with `by = join_by(curve_id)`

``` r
odd.curves.p
```

![](02-growth_filtering_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
ggsave(path_target("odd_curves.png"), plot=odd.curves.p, dpi=300)
```

    ## Saving 7 x 5 in image

``` r
odd_curves <- c("PgcvB_lp7_D5_r1", "PgcvB_lp7_D5_r2", "PsroC_lp2_A2_r2", "PsroC_lp2_A2_r1")

filter.metadata <- filter.metadata %>% 
  mutate(filter.status = case_when(filter.status != "remove" & curve_id %in% odd_curves ~ "remove",
                                   TRUE ~ filter.status),
         
         filter.reason = case_when(filter.reason == "None" &  curve_id %in% odd_curves ~ "odd.growth",
                                   TRUE ~ filter.reason))
```

``` r
filter.metadata %>% count(filter.reason)
```

    ##               filter.reason      n
    ## 1                      None 155322
    ## 2           killer.compound   5580
    ## 3                low.growth    268
    ## 4 nrep_after_growth_filters    106
    ## 5                odd.growth      4

## Dissimilar growth.

Here we remove experiments where there is a large difference in OD.

``` r
# OD at 10 hours
od_at_10h.complete <- salm.dgobj.linear@od_data %>% 
  mutate(diff_10h = abs(timepoint - 10)) %>% 
  
  group_by(curve_id) %>% 
  filter(diff_10h == min(diff_10h)) %>% 
  ungroup() %>% 
  
  select(curve_id, timepoint, od) %>% 
  rename("od_10h" = "od",
         "timepoint_10h" = "timepoint")
```

``` r
kept.experiments <- filter.metadata %>% filter(filter.status == "keep") %>% pull(curve_id)

od_10h.wide <- od_at_10h.complete %>% 
  filter(curve_id %in% kept.experiments) %>% 
  
  select(curve_id, od_10h) %>% 
  
  separate(curve_id, c("promoter", "libplate", "well", "replicate")) %>% 
  unite("srn_code", c(libplate, well)) %>% 
  
  pivot_wider(id_cols = c(promoter, srn_code), names_from = replicate, values_from = od_10h) %>% 
  mutate(od_diff = abs(r1-r2))
```

``` r
od.comparison.p <- ggplot(od_10h.wide, aes(x=r1, y=r2, color=od_diff)) +
  facet_wrap(~promoter, ncol=6) +
  geom_point(size=1) +
  
  geom_abline(color="orange", linetype="longdash") +
  geom_smooth(formula = y~x, method="lm", color="red") +
  
  theme_bw() +
  ggtitle("OD at 10h comparison") +
  
  labs(x="Replicate 1",
       y="Replicate 2")

od.comparison.p
```

![](02-growth_filtering_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
ggsave(path_target("od_comparison.png"), plot=od.comparison.p, dpi=300)
```

    ## Saving 7 x 5 in image

We can eliminate those at an od diff of 0.5

``` r
od_10h.wide <- od_10h.wide %>% 
  mutate(removal_candidate = if_else(od_diff > 0.5, "Remove", "Remain"))


od.removal.p <- ggplot(od_10h.wide, aes(x=r1, y=r2, color=removal_candidate))+
  
  facet_wrap(~promoter, ncol=6) +
  geom_point(size=1) +
  scale_color_manual(values = c(alpha("#C5C5C5", 0.3), "#d95f02")) +
  geom_abline(color="red") +
  theme_light() +
  ggtitle("OD at 10h comparison")

od.removal.p
```

![](02-growth_filtering_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
ggsave(path_target("od_removal.png"), plot=od.removal.p, dpi=300)
```

    ## Saving 7 x 5 in image

``` r
example.diffOD.p <- salm.dgobj.linear@od_data %>% 
  filter(curve_id %in% c("PEVC_lp1_E2_r1", "PEVC_lp1_E2_r2")) %>% 
  left_join(salm.dgobj.linear@metadata %>% select(curve_id, promoter, srn_code, replicate)) %>% 
  
  unite("promoter_compound", c(promoter, srn_code)) %>% 
  
  
  ggplot(aes(x=timepoint, y=od, group=curve_id, color=replicate)) +
  geom_line() +
  
  facet_wrap(~promoter_compound) +
  
  theme_bw() +
  
  labs(x="Time (h)",
       y="OD")
```

    ## Joining with `by = join_by(curve_id)`

``` r
example.diffOD.p
```

![](02-growth_filtering_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
ggsave(path_target("example_diffOD.png"), plot=example.diffOD.p, dpi=300)
```

    ## Saving 7 x 5 in image

``` r
remove.promchem <- od_10h.wide %>% 
  unite("promoter_compound", c(promoter, srn_code)) %>% 
  filter(removal_candidate == "Remove") %>% 
  pull(promoter_compound)


filter.metadata <- filter.metadata %>% 
  unite("promoter_compound", c(promoter, srn_code), remove=FALSE) %>% 
  
  mutate(filter.status = case_when(filter.status != "remove" & promoter_compound %in% remove.promchem ~ "remove",
                                   TRUE ~ filter.status),
         
         filter.reason = case_when(filter.reason == "None" &  promoter_compound %in% remove.promchem ~ "large.od.diff",
                                   TRUE ~ filter.reason))

filter.metadata %>% count(filter.reason)
```

    ##               filter.reason      n
    ## 1                      None 155274
    ## 2           killer.compound   5580
    ## 3             large.od.diff     48
    ## 4                low.growth    268
    ## 5 nrep_after_growth_filters    106
    ## 6                odd.growth      4

``` r
filter.metadata %>% dim()
```

    ## [1] 161280     12

## Write files.

``` r
write_tsv(filter.metadata, path_target("filter_metadata.tsv.gz"))
write_tsv(od_at_10h.complete, path_target("od_10hours.tsv.gz"))
```

## Files written

These files have been written to the target directory,
`data/02-growth_filtering`:

``` r
projthis::proj_dir_info(path_target())
```

    ## # A tibble: 7 × 4
    ##   path                   type         size modification_time  
    ##   <fs::path>             <fct> <fs::bytes> <dttm>             
    ## 1 example_diffOD.png     file       83.61K 2025-11-18 19:30:59
    ## 2 filter_histogram.png   file        49.5K 2025-11-18 19:30:38
    ## 3 filter_metadata.tsv.gz file        1.41M 2025-11-18 19:31:03
    ## 4 od_10hours.tsv.gz      file      729.21K 2025-11-18 19:31:07
    ## 5 od_comparison.png      file      570.44K 2025-11-18 19:30:51
    ## 6 od_removal.png         file      479.56K 2025-11-18 19:30:58
    ## 7 odd_curves.png         file      114.04K 2025-11-18 19:30:42
