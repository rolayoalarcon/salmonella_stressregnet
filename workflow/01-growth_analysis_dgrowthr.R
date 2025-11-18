library(tidyverse)
library(DGrowthR)

closeAllConnections()
gc()

salm.dgobj <- readRDS("data/01-growth_analysis/salm_dgobj_log.rds")


salm.dgobj <- estimate_growth_parameters(salm.dgobj,
                                         od_auc_at_t = 10,
                                         save_gp_data = TRUE,
                                         n_cores = 80)


saveRDS(salm.dgobj, "tmp/salm_dgobj_LOG_gparams_18112025.rds")

############ Linear

salm.dgobj <- readRDS("data/01-growth_analysis/salm_dgobj_linear.rds")


salm.dgobj <- estimate_growth_parameters(salm.dgobj,
                                         od_auc_at_t = 10,
                                         save_gp_data = TRUE,
                                         n_cores = 80)


saveRDS(salm.dgobj, "tmp/salm_dgobj_LIN_gparams_18112025.rds")
