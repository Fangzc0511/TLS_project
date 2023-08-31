
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(stringr))

cancer_type <- list.files("./TLS_project/data/ICB_pheno") %>% strsplit("\\.") %>% sapply("[", 1)
ICB_pheno <- list()

cancer_type <- cancer_type[-which(cancer_type=="Mariathasan_UC")]
pheno <- read.table("./TLS_project/data/ICB_pheno/Mariathasan_UC.tsv", header = T, sep = "\t")
pheno_filter <- subset(pheno, select = c("ID", "OS", "status", "response"))
#pheno_filter <- cbind(rownames(pheno), pheno_filter)
colnames(pheno_filter) <- c("sample_id", "os", "status", "response")
pheno_filter$os <- as.numeric(pheno_filter$os)
pheno_filter$status <- as.numeric(pheno_filter$status)
pheno_filter$response <- as.numeric(pheno_filter$response)
ICB_pheno[["Mariathasan_UC"]] <- pheno_filter

cancer_type <- cancer_type[-which(cancer_type=="miao_1_CRCC_1")]
pheno <- read.table("./TLS_project/data/ICB_pheno/miao_1_CRCC_1.tsv", header = T, sep = "\t")
pheno_filter <- subset(pheno, select = c("patient_id", "os_days", "os_censor", "response"))
colnames(pheno_filter) <- c("sample_id", "os", "status", "response")
pheno_filter$os <- as.numeric(pheno_filter$os)
pheno_filter$status <- as.numeric(pheno_filter$status)
pheno_filter$response <- as.numeric(pheno_filter$response)
ICB_pheno[["miao_1_CRCC_1"]] <- pheno_filter

cancer_type <- cancer_type[-which(cancer_type=="Puch_melanoma")]
pheno <- read.table("./TLS_project/data/ICB_pheno/Puch_melanoma.tsv", header = T, sep = "\t")
pheno_filter <- subset(pheno, select = c("X", "OS", "status", "response"))
colnames(pheno_filter) <- c("sample_id", "os", "status", "response")
pheno_filter[which(pheno_filter$response==-1), "response"] <- 1
pheno_filter$os <- as.numeric(pheno_filter$os)
pheno_filter$status <- as.numeric(pheno_filter$status)
pheno_filter$response <- as.numeric(pheno_filter$response)
pheno_filter$sample_id <- chartr("-",".", pheno_filter$sample_id)
pheno_filter$sample_id <- paste0("X", pheno_filter$sample_id)
ICB_pheno[["Puch_melanoma"]] <- pheno_filter

cancer_type <- cancer_type[-which(cancer_type=="Snyder_UC")]
pheno <- read.table("./TLS_project/data/ICB_pheno/Snyder_UC.tsv", header = T, sep = "\t")
pheno_filter <- subset(pheno, select = c("ID", "OS", "status", "response"))
colnames(pheno_filter) <- c("sample_id", "os", "status", "response")
pheno_filter$os <- as.numeric(pheno_filter$os)
pheno_filter$status <- as.numeric(pheno_filter$status)
pheno_filter$response <- as.numeric(pheno_filter$response)
ICB_pheno[["Snyder_UC"]] <- pheno_filter

for (i in 1:length(cancer_type)){
    pheno <- read.table(paste0("./TLS_project/data/ICB_pheno/", cancer_type[i], ".tsv"), header = T, sep = "\t")
    #pheno_filter <- subset(pheno, select = c("overall.survival..days.", "vital.status", "response_NR"))
    pheno_filter <- data.frame(sample_id=pheno$sample_id, os=pheno$overall.survival..days., status=pheno$vital.status, response=pheno$response_NR)
    pheno_filter$status <- lapply(pheno_filter$status, function(x){ifelse(x == "Dead", 1, 0)})
    pheno_filter$response <- lapply(pheno_filter$response, function(x){ifelse(x == "R", 1, ifelse(x == "N", 0, NA))})
    #print(head(pheno_filter, n=5))
    pheno_filter$os <- as.numeric(pheno_filter$os)
    pheno_filter$status <- as.numeric(pheno_filter$status)
    pheno_filter$response <- as.numeric(pheno_filter$response)
    if (cancer_type[i]=="Nathanson_Melanoma_2017"){pheno_filter$sample_id <- paste0("X", pheno_filter$sample_id)}
    ICB_pheno[[cancer_type[i]]] <- pheno_filter
    
}

save(ICB_pheno, file = "./TLS_project/data/ICB_pheno.Rdata")

cancer_type <- list.files("./TLS_project/data/ICB_expr") %>% strsplit("\\.") %>% sapply("[", 1)

ICB_expr <- list()

for (i in 1:length(cancer_type)){
    expr <- read.table(paste0("./TLS_project/data/ICB_expr/", cancer_type[i], ".tsv"), header = T, row.names = 1, sep = "\t") 
    ICB_expr[[cancer_type[i]]] <- expr
}

save(ICB_expr, file = "./TLS_project/data/ICB_expr.Rdata")

cancer_type <- list.files("./TLS_project/data/TCGA_expr") %>% strsplit("\\.") %>% sapply("[", 1)
TCGA_pheno <- list()
TCGA_expr <- list()

for (i in 1:length(cancer_type)){
    expr <- read.table(paste0("./TLS_project/data/TCGA_expr/", cancer_type[i], ".txt"), header = T, row.names = 1, sep = "\t") 
    colnames(expr) <- chartr(".","-", colnames(expr))
    TCGA_expr[[cancer_type[i]]] <- expr
    
    pheno <- read.table(paste0("./TLS_project/data/TCGA_pheno/survival%2F", cancer_type[i], "_survival.txt"), header = T, row.names = 1, sep = "\t") 
    pheno_filter <- data.frame(sample_id=rownames(pheno), os=pheno$OS.time, os.status=pheno$OS, pfi=pheno$PFI.time, pfi.status=pheno$PFI)
    TCGA_pheno[[cancer_type[i]]] <- pheno_filter
}

save(TCGA_expr, file = "./TLS_project/data/TCGA_expr.Rdata")
save(TCGA_pheno, file = "./TLS_project/data/TCGA_pheno.Rdata")

cancer_type <- list.files("./TLS_project/data/TCGA_expr") %>% strsplit("\\.") %>% sapply("[", 1)
TCGA_TMB <- list()

for (i in 1:length(cancer_type)){
    data <- fread(paste0("./TLS_project/data/TCGA_TMB/", cancer_type[i], "_mc3.txt"), header = T, stringsAsFactors = F)
    data_filter <- data[which(data$effect %in% c("Missense_Mutation","Frame_Shift_Del","Nonstop_Mutation")),]
    sample_id <- unique(data_filter$sample)
    TMB <- c()
    for (j in 1:length(sample_id)){
        data_per_sample <- data_filter[which(data_filter$sample==sample_id[j]),]
        TMB <- c(TMB, nrow(data_per_sample)/38)
    }
    TMB_data <- data.frame(sample_id=sample_id, TMB=TMB)
    TCGA_TMB[[cancer_type[i]]] <- TMB_data  
}

save(TCGA_TMB, file = "./TLS_project/data/TCGA_TMB.Rdata")

gene_signatures <- list()
gene_signatures_names <- list.files("./TLS_project/data/gene_signatures") %>% strsplit("\\.") %>% sapply("[", 1)
for (i in 1:length(gene_signatures_names)){
    gene_list <- read.table(paste0("./TLS_project/data/gene_signatures/", gene_signatures_names[i], ".txt"))
    gene_signatures[[gene_signatures_names[i]]] <- gene_list  
}

save(gene_signatures, file = "./TLS_project/data/gene_signatures.Rdata")



source("./TLS_project/code/get_infiltration/get_infiltration.R")

cancer_type <- list.files("./TLS_project/data/ICB_infiltration") %>% strsplit("_results") %>% sapply("[", 1)
ICB_infiltration <- list()

for (i in 1:length(cancer_type)){
    infiltration_data <- read.table(paste0("./TLS_project/data/ICB_infiltration/", cancer_type[i], "_results.txt"), header = T, sep = "\t")
    infiltration_data <- subset(infiltration_data, select = c("B.cells.naive", "B.cells.memory",
                                                            "T.cells.CD8", "T.cells.CD4.naive",
                                                            "T.cells.CD4.memory.resting", "T.cells.CD4.memory.activated",
                                                            "T.cells.follicular.helper", "T.cells.regulatory..Tregs.",
                                                            "T.cells.gamma.delta"))
    average <- apply(infiltration_data, 1, sum)
    if(cancer_type[i]=="Puch_melanoma"){
        names(average) <- chartr("-",".", names(average))
        names(average) <- paste0("X", names(average))
    }
    ICB_infiltration[[cancer_type[i]]] <- data.frame(sample_id=names(average), infiltration=average)
    
}

save(ICB_infiltration, file = "./TLS_project/data/ICB_infiltration.Rdata")

cancer_type <- list.files("./TLS_project/data/TCGA_infiltration") %>% strsplit("_results") %>% sapply("[", 1)
TCGA_infiltration <- list()

for (i in 1:length(cancer_type)){
    infiltration_data <- read.table(paste0("./TLS_project/data/TCGA_infiltration/", cancer_type[i], "_results.txt"), header = T, sep = "\t")
    infiltration_data <- subset(infiltration_data, select = c("B.cells.naive", "B.cells.memory",
                                                            "T.cells.CD8", "T.cells.CD4.naive",
                                                            "T.cells.CD4.memory.resting", "T.cells.CD4.memory.activated",
                                                            "T.cells.follicular.helper", "T.cells.regulatory..Tregs.",
                                                            "T.cells.gamma.delta"))
    average <- apply(infiltration_data, 1, sum)
    TCGA_infiltration[[cancer_type[i]]] <- data.frame(sample_id=rownames(infiltration_data), infiltration=average)
}

save(TCGA_infiltration, file = "./TLS_project/data/TCGA_infiltration.Rdata")
