04-hit_determination
================
Compiled at 2025-11-18 22:16:34 UTC

``` r
here::i_am(paste0(params$name, ".Rmd"), uuid = "7f90b0e4-3dbc-4a8f-a2d8-1ac4b37dd6c3")
```

The purpose of this document is …

``` r
library("tidyverse")
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
# create or *empty* the target directory, used to write this file's data: 
projthis::proj_create_dir_target(params$name, clean = TRUE)

# function to get path to target directory: path_target("sample.csv")
path_target <- projthis::proj_path_target(params$name)

# function to get path to previous data: path_source("00-import", "sample.csv")
path_source <- projthis::proj_path_source(params$name)
```

## Read data.

``` r
expression_data <- read_tsv(path_source("03-estimate_lux", "lux_auc_filtered_median.tsv.gz"), show_col_types = FALSE) %>% 
  
  mutate(expression_value = log2.lux.normed.centered)

libplate_map <- read_tsv(path_source("00-import/Salmonella", "LibMap.tsv.gz"), show_col_types = FALSE) %>% 
  mutate(libplate = str_replace(`Library plate`, "LibPlate", "lp"),
         srn_code = paste(libplate, `New well`, sep="_")) %>% 
  select(srn_code, `Catalog Number`)
```

## Remove DMSO noisy experiments

``` r
## Determine DMSO noisy curves
dmso_noisy.srn <- libplate_map %>% 
  filter(`Catalog Number` == "DMSO noisy") %>% 
  select(srn_code) %>% 
  unlist()
  


expression_data <- expression_data %>% 
  filter(!srn_code %in% dmso_noisy.srn)
```

## Calculate change w.r.t. EVC

``` r
evc.fc.10h <- expression_data %>% 
  
  unite("promoter_replicate", c(promoter, replicate), sep="_") %>% 
  
  pivot_wider(id_cols=srn_code, names_from = promoter_replicate, values_from = expression_value) %>% 
  
  # Remove cases where EVC3 is missing. This is 6 compounds
  filter(!is.na(PEVC3_r1) | !is.na(PEVC3_r2)) %>% 
  
 # Get an average EVC expression
  mutate(PEVC3_avg = rowMeans(select(., PEVC3_r1, PEVC3_r2))) %>% 
  
  # Prepare dataframe for normalisation
  select(-c(PEVC3_r1, PEVC3_r2)) %>% 
  column_to_rownames("srn_code") %>%
    
  # Normalise with respect to EVC average
  mutate(across(.cols = everything(), .fns = ~. - PEVC3_avg)) %>% 
    
  # Convert back to long format
  select(-PEVC3_avg) %>% 
  rownames_to_column("srn_code") %>% 
  
  
  pivot_longer(cols=-srn_code, names_to = "promoter_replicate", values_to = "log.evcfc") %>% 
  filter(!is.na(log.evcfc))
```

## Test for normality.

``` r
promoter.normtest <- function(promrep_str, express.val.str, auc_df){
  
  # Gather the relevant information
  promrep.aucs <- auc_df %>% 
    filter(promoter_replicate == promrep_str) %>% 
    select(all_of(express.val.str)) %>%
    unlist()
  
  # Prepare output
  r1.df <- data.frame(promoter_replicate=promrep_str,
                      shapiro.pval=shapiro.test(promrep.aucs)$p.value,
                      lillie.pval=nortest::lillie.test(promrep.aucs)$p.value)
  return(r1.df)
  }
```

Gather DMSO values.

``` r
dmso.evcfc <- evc.fc.10h %>% 
  left_join(libplate_map %>% select(srn_code, `Catalog Number`)) %>% 
  filter(`Catalog Number` == "DMSO")
```

    ## Joining with `by = join_by(srn_code)`

``` r
dmso.evcfc %>% 
  group_by(promoter_replicate) %>% 
  summarise(n_dmso_wells = n())
```

    ## # A tibble: 58 × 2
    ##    promoter_replicate n_dmso_wells
    ##    <chr>                     <int>
    ##  1 PEVC2_r1                    126
    ##  2 PEVC2_r2                    126
    ##  3 PEVC_r1                     126
    ##  4 PEVC_r2                     126
    ##  5 PacrAB_r1                   126
    ##  6 PacrAB_r2                   126
    ##  7 ParcZ_r1                    126
    ##  8 ParcZ_r2                    126
    ##  9 PcpxP_r1                    126
    ## 10 PcpxP_r2                    126
    ## # ℹ 48 more rows

``` r
salm.normtest.pvals <- bind_rows(lapply(unique(dmso.evcfc$promoter_replicate), promoter.normtest, express.val.str="log.evcfc", auc_df=dmso.evcfc)) %>% 
  
  mutate(shapiro.pval.adj = p.adjust(shapiro.pval, method = "BH"),
         lillie.pval.adj = p.adjust(lillie.pval, method = "BH")) %>% 
  
  separate(promoter_replicate, into=c("promoter", "replicate"), sep="_") %>%
  mutate(replicate = str_split(replicate, "r", simplify = TRUE)[,2],
         replicate = as.numeric(replicate),
         replicate = factor(replicate, levels=unique(sort(replicate))))



salm.normtest.pvals.long <- salm.normtest.pvals %>% 
  pivot_wider(names_from = replicate, values_from = c(shapiro.pval.adj, lillie.pval.adj), id_cols = promoter) %>% 
  pivot_longer(cols=-promoter, names_to="test_replicate", values_to="adj.pval") %>% 
  mutate(decision = if_else(adj.pval < 0.05, "Not normal", "Normal"),
         decision = if_else(is.na(adj.pval), "NA", decision),
         decision = factor(decision, levels=c("Normal", "Not normal", "NA")),
         
         pval_label = if_else(is.na(adj.pval), "NA", as.character(round(adj.pval, 2)))) %>% 
  
  separate(test_replicate, into=c("test", "replicate"), sep="_") %>%
  mutate(replicate = factor(replicate, levels=unique(sort(replicate))))
```

``` r
lillie.salm.g <- salm.normtest.pvals.long %>% 
  mutate(promoter = case_when(promoter == "PEVC" ~ "EVC1",
                              promoter == "PEVC2" ~ "EVC2",
                              promoter == "PEVC3" ~ "EVC3",
                              promoter == "PmarR" ~ "PmarRAB",
                              TRUE ~ promoter),
         promoter = factor(promoter, levels = rev(sort(unique(promoter))))) %>%
  filter(test == "lillie.pval.adj") %>% 
  
   ggplot(aes(x=replicate, y=promoter, fill = decision)) +
  geom_tile(color="black") +
  geom_text(aes(label = pval_label), size=3) +
  
  scale_fill_manual(values=c(alpha("#CC79A7", 0.95), alpha("#b3b3b3", 0.65)))+

  theme_minimal() +
  labs(title = "Lilliefors test: adj. p-values",
       fill="Decision",
       x="Replicate",
       y="Promoter") +
  theme(legend.position = "right",
        
        legend.key.size = unit(0.3, "cm"))

lillie.salm.g
```

![](04-hit_determination_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
ggsave(path_target("lillie_dmso_tests.png"), plot=lillie.salm.g, dpi=300, width = 5)
```

    ## Saving 5 x 5 in image

``` r
dnorm_promrep <- function(promrep_str, express.val.str, auc_df){
  
  # Gather the relevant information
  promrep.aucs <- auc_df %>% 
    filter(promoter_replicate == promrep_str) %>% 
    select(all_of(express.val.str)) %>%
    unlist()
  
  # Prepare h
  h <- hist(promrep.aucs, breaks = 30, density = 10)
  
  xfit <- seq(min(promrep.aucs), max(promrep.aucs), length = 40) 
  
  yfit <- dnorm(xfit, mean = mean(promrep.aucs), sd = sd(promrep.aucs)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(promrep.aucs) 
  
  # Prepare output
  df.out <- data.frame(promoter_replicate=promrep_str,
                      x=xfit,
                      y=yfit)
  
  return(df.out)
}


salm.normtest.dnorm <- bind_rows(lapply(unique(dmso.evcfc$promoter_replicate), dnorm_promrep, express.val.str="log.evcfc", auc_df=dmso.evcfc)) %>% 
  
  separate(promoter_replicate, into=c("promoter", "replicate"), sep="_") %>%
  mutate(replicate = str_split(replicate, "r", simplify = TRUE)[,2],
         replicate = as.numeric(replicate),
         replicate = factor(replicate, levels=unique(sort(replicate))))
```

![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-4.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-5.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-6.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-7.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-8.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-9.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-10.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-11.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-12.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-13.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-14.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-15.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-16.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-17.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-18.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-19.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-20.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-21.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-22.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-23.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-24.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-25.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-26.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-27.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-28.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-29.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-30.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-31.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-32.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-33.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-34.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-35.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-36.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-37.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-38.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-39.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-40.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-41.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-42.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-43.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-44.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-45.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-46.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-47.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-48.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-49.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-50.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-51.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-52.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-53.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-54.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-55.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-56.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-57.png)<!-- -->![](04-hit_determination_files/figure-gfm/unnamed-chunk-9-58.png)<!-- -->

``` r
df_ <- dmso.evcfc %>% 
  separate(promoter_replicate, into=c("promoter", "replicate"), sep="_") %>% 
  mutate(replicate = str_split(replicate, "r", simplify = TRUE)[,2],
         replicate = as.numeric(replicate),
         replicate = factor(replicate, levels=unique(sort(replicate))))
  

distrib.dens.salm.g <- ggplot(df_) +
  geom_histogram(aes(x=log.evcfc, fill=replicate), alpha=0.2, bins = 50, color="black",
                 position = "identity", linewidth=0.01) +
  geom_line(data=salm.normtest.dnorm, aes(x, y, color=replicate), linewidth=0.5) +
  
  facet_wrap(~promoter, scales="free_y") +
  theme_bw() +
  
  labs(x=latex2exp::TeX("$log_2$ Norm. Lux $- log_2$ average EVC"),
       y = "Count \\ Density") +
  theme(text=element_text(size=7),
        axis.text = element_text(size=5),
        axis.text.x = element_text(size=5, angle = 30, hjust = 1),
        legend.position = "bottom",
        legend.key.size = unit(0.1, "cm")) +
  labs(color="Replicate",
       fill = "Replicate",
       title = "DMSO experiments")
distrib.dens.salm.g
```

![](04-hit_determination_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
ggsave(plot = distrib.dens.salm.g, filename = path_target("dmso_distribution.pdf"), dpi=300,
       width = 20, height = 12, units="cm")
```

## Pvalue calculation

``` r
dmso.params <- dmso.evcfc %>% 
  group_by(promoter_replicate) %>% 
  summarise(dmso.mean = mean(log.evcfc),
            dmso.stdv = sd(log.evcfc))

pvalues.df <- evc.fc.10h %>% 
  left_join(dmso.params) %>% 
  
  mutate(zscore = (log.evcfc-dmso.mean)/dmso.stdv,
         pvalue = 2 * pnorm(q=abs(zscore), mean=0, sd=1, lower.tail=FALSE)) %>% 
  separate(promoter_replicate, c("promoter", "replicate"), remove = FALSE)
```

    ## Joining with `by = join_by(promoter_replicate)`

``` r
pval.histogram <- pvalues.df %>% 
  
   mutate(promoter = case_when(promoter == "PEVC" ~ "EVC1",
                              promoter == "PEVC2" ~ "EVC2",
                              promoter == "PEVC3" ~ "EVC3",
                              promoter == "PmarR" ~ "PmarRAB",
                              TRUE ~ promoter)) %>% 
  
  ggplot(aes(x=pvalue, fill=replicate)) +
  geom_histogram(binwidth = 0.01, position="identity", alpha=0.5) +
  facet_wrap(~promoter, scales = "free_y") +
  
    theme_bw() +
  
  theme(text = element_text(size=6),
        legend.position = "right") +

  labs(x="p-value",
       y="Count",
       fill="Replicate")

pval.histogram
```

![](04-hit_determination_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
ggsave(plot = pval.histogram, filename = path_target("pvalue_histogram.pdf"), dpi=300,
       width = 18, height = 12, units="cm")
```

## Hit determination.

``` r
pvalues_adj.df <- pvalues.df %>% 
  
  group_by(promoter, srn_code) %>% 
  slice_max(pvalue, n=1) %>% 
  ungroup() %>% 
  
  group_by(promoter) %>% 
  mutate(pvalue.adj = p.adjust(pvalue, method = "BH")) %>% 
  ungroup()
```

``` r
hit_table <- pvalues_adj.df %>% 
   mutate(hit = case_when(pvalue.adj < 0.01 & zscore > 0 ~ "Upregulated",
                                pvalue.adj < 0.01 & zscore < 0 ~ "Downregulated",
                                TRUE ~ "Not DE"))
```

``` r
hit_table %>% 
  filter(hit != "Not DE") %>% 

  count(hit) %>% 
  knitr::kable()
```

    ## Warning: 'xfun::attr()' is deprecated.
    ## Use 'xfun::attr2()' instead.
    ## See help("Deprecated")
    ## Warning: 'xfun::attr()' is deprecated.
    ## Use 'xfun::attr2()' instead.
    ## See help("Deprecated")

| hit           |   n |
|:--------------|----:|
| Downregulated |  96 |
| Upregulated   | 520 |

``` r
hit_table %>% 
  filter(hit != "Not DE") %>% 
  group_by(promoter) %>% 
  
  count(hit) %>% 
  knitr::kable()
```

    ## Warning: 'xfun::attr()' is deprecated.
    ## Use 'xfun::attr2()' instead.
    ## See help("Deprecated")
    ## Warning: 'xfun::attr()' is deprecated.
    ## Use 'xfun::attr2()' instead.
    ## See help("Deprecated")

| promoter | hit           |   n |
|:---------|:--------------|----:|
| PEVC     | Upregulated   |   1 |
| PEVC2    | Upregulated   |   1 |
| PacrAB   | Downregulated |   1 |
| PacrAB   | Upregulated   |  55 |
| PcpxP    | Downregulated |   6 |
| PcpxP    | Upregulated   |  14 |
| PdsrA    | Upregulated   |   6 |
| PfeoB    | Downregulated |   1 |
| PfeoB    | Upregulated   |  11 |
| PgcvB    | Downregulated |   6 |
| PgcvB    | Upregulated   |   3 |
| PglmY    | Upregulated   |   4 |
| PiscR    | Upregulated   |  12 |
| PisrK    | Downregulated |   1 |
| PisrK    | Upregulated   |   4 |
| PkatE    | Downregulated |   6 |
| PkatE    | Upregulated   |   1 |
| PmarR    | Downregulated |   1 |
| PmarR    | Upregulated   |  99 |
| PmicA    | Downregulated |  16 |
| PmicA    | Upregulated   |  11 |
| PmicF    | Downregulated |   1 |
| PmicF    | Upregulated   |  45 |
| Pnc1390  | Downregulated |   1 |
| Pnc1390  | Upregulated   |   5 |
| Pnc700   | Downregulated |   1 |
| Pnc700   | Upregulated   |   5 |
| PomrB    | Upregulated   |   8 |
| PoxyS    | Upregulated   |   9 |
| PpinT    | Upregulated   |  23 |
| PramA    | Downregulated |   3 |
| PramA    | Upregulated   | 137 |
| Prob     | Upregulated   |   2 |
| PryhB1   | Upregulated   |  28 |
| PsdsR    | Downregulated |   9 |
| PsdsR    | Upregulated   |   2 |
| PsoxS    | Downregulated |   1 |
| PsoxS    | Upregulated   |  12 |
| PsroC    | Downregulated |   4 |
| PsroC    | Upregulated   |   1 |
| PvirK    | Downregulated |  38 |
| PvirK    | Upregulated   |  21 |

## Additional data.

``` r
chemical_information <- read_tsv(path_source("00-import", "chemical_library/chemical_library.tsv.gz"), show_col_types = FALSE)
```

    ## New names:
    ## • `` -> `...1`

``` r
chemical_information <- chemical_information %>% 
  select(srn_code, ProductName, antibiotic,chemtype_mce,ATC_level1,ATC_level2)


hit_table.metadata <- hit_table %>% 
  select(-c(promoter_replicate, replicate)) %>% 
  
  left_join(chemical_information, by="srn_code")
```

``` r
hit_table.metadata %>% 
  filter(hit != "Not DE",
         promoter == "PvirK") %>% 
  arrange(hit)
```

    ## # A tibble: 59 × 14
    ##    srn_code promoter log.evcfc dmso.mean dmso.stdv zscore   pvalue pvalue.adj
    ##    <chr>    <chr>        <dbl>     <dbl>     <dbl>  <dbl>    <dbl>      <dbl>
    ##  1 lp1_H12  PvirK       -0.550   0.0116      0.141  -4.00 6.44e- 5   3.02e- 3
    ##  2 lp1_I4   PvirK       -0.572   0.0116      0.141  -4.15 3.28e- 5   1.67e- 3
    ##  3 lp1_O19  PvirK       -0.565  -0.00771     0.135  -4.13 3.65e- 5   1.82e- 3
    ##  4 lp1_O22  PvirK       -1.15    0.0116      0.141  -8.28 1.27e-16   1.73e-14
    ##  5 lp1_P2   PvirK       -1.27    0.0116      0.141  -9.09 9.73e-20   1.70e-17
    ##  6 lp2_G15  PvirK       -0.671  -0.00771     0.135  -4.91 8.98e- 7   6.26e- 5
    ##  7 lp2_N2   PvirK       -0.561   0.0116      0.141  -4.08 4.58e- 5   2.20e- 3
    ##  8 lp2_N22  PvirK       -0.549   0.0116      0.141  -3.99 6.69e- 5   3.08e- 3
    ##  9 lp3_G19  PvirK       -0.968   0.0116      0.141  -6.97 3.25e-12   3.30e-10
    ## 10 lp3_G21  PvirK       -0.542  -0.00771     0.135  -3.96 7.55e- 5   3.41e- 3
    ## # ℹ 49 more rows
    ## # ℹ 6 more variables: hit <chr>, ProductName <chr>, antibiotic <lgl>,
    ## #   chemtype_mce <chr>, ATC_level1 <chr>, ATC_level2 <chr>

## Volcano plot.

``` r
ggplot(hit_table.metadata, aes(x=zscore, y=-log10(pvalue.adj), color=zscore)) +
  
  geom_point(size=1) +
  
  theme_bw()
```

![](04-hit_determination_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

## Write files.

``` r
write_tsv(hit_table.metadata, path_target("hit_table.tsv.gz"))
```

## Files written

These files have been written to the target directory,
`data/04-hit_determination`:

``` r
projthis::proj_dir_info(path_target())
```

    ## # A tibble: 4 × 4
    ##   path                  type         size modification_time  
    ##   <fs::path>            <fct> <fs::bytes> <dttm>             
    ## 1 dmso_distribution.pdf file        33.9K 2025-11-18 22:16:44
    ## 2 hit_table.tsv.gz      file        3.47M 2025-11-18 22:17:05
    ## 3 lillie_dmso_tests.png file      209.76K 2025-11-18 22:16:37
    ## 4 pvalue_histogram.pdf  file       39.92K 2025-11-18 22:16:49
