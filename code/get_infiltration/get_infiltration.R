suppressMessages(library(CIBERSORT))
suppressMessages(library(tidyverse))


###TCGA
cancer_type <- list.files("./TLS_project/data/TCGA_expr") %>% strsplit("\\.") %>% sapply("[", 1)
for(i in cancer_type){
  data <- read.table(paste0("./TLS_project/data/TCGA_expr/", i, ".txt"), header = T, sep = "\t")
  sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
  results <- CIBERSORT::cibersort(sig_matrix, 
                                  paste0("./TLS_project/data/TCGA_expr/", i, ".txt"),
                                  perm = 1000,
                                  QN = T)
  write.table(results, file = paste0("./TLS_project/data/TCGA_infiltration/", i, "_result.txt"), quote = FALSE, row.names = TRUE, sep = "\t")
}

###ICB
cancer_type <- list.files("./TLS_project/data/ICB_expr") %>% strsplit("\\.") %>% sapply("[", 1)

for(i in cancer_type){
  data <- read.table(paste0("./TLS_project/data/TCGA_expr/", i, ".txt"), header = T, sep = "\t")
  sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
  results <- CIBERSORT::cibersort(sig_matrix, paste0("./TLS_project/data/ICB_expr/", i, ".txt"),
                                  perm = 1000,
                                  QN = T)
  write.table(results, file = paste0("./zcfang/TLS_project/data/ICB_infiltration/", i, "_results.txt"), quote = FALSE, row.names = TRUE, sep = "\t")
}
