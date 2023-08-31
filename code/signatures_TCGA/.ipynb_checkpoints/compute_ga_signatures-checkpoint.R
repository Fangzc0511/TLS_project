{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "suppressMessages(library(GSVA))\n",
    "suppressMessages(library(survcomp))\n",
    "suppressMessages(library(tidyverse))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [],
   "source": [
    "create_directory <- function(signature){\n",
    "    dir <- paste0(\"./TLS_project/results/siganatures_TCGA/\", signature)\n",
    "    if (file.exists(dir)){unlink(dir, recursive = TRUE)}\n",
    "    dir.create(dir, recursive = TRUE)\n",
    "    dir.create(paste0(dir, \"/KMPlot\"))\n",
    "    dir.create(paste0(dir, \"/KMPlot/OS\"))\n",
    "    dir.create(paste0(dir, \"/KMPlot/PFI\"))\n",
    "    dir.create(paste0(dir, \"/Overall\"))\n",
    "    dir.create(paste0(dir, \"/Overall/OS\"))\n",
    "    dir.create(paste0(dir, \"/Overall/PFI\"))\n",
    "}"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {},
   "outputs": [],
   "source": [
    "source(\"./TLS_project/code/summary_figures/Get_KMplot.R\")\n",
    "source(\"./TLS_project/code/meta_analysis/Get_Association.R\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "metadata": {},
   "outputs": [],
   "source": [
    "load(\"./TLS_project/data/TCGA_expr.Rdata\")\n",
    "load(\"./TLS_project/data/TCGA_pheno.Rdata\")\n",
    "load(\"./TLS_project/data/gene_signatures.Rdata\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "scrolled": true
   },
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Warning message in .filterFeatures(expr, method):\n",
      "“1 genes with constant expression values throuhgout the samples.”\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Estimating ssGSEA scores for 1 gene sets.\n",
      "[1] \"Calculating ranks...\"\n",
      "[1] \"Calculating absolute values from ranks...\"\n",
      "  |======================================================================| 100%\n",
      "\n",
      "[1] \"Normalizing...\"\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Warning message in .filterFeatures(expr, method):\n",
      "“1 genes with constant expression values throuhgout the samples.”\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Estimating ssGSEA scores for 1 gene sets.\n",
      "[1] \"Calculating ranks...\"\n",
      "[1] \"Calculating absolute values from ranks...\"\n",
      "  |======================================================================| 100%\n",
      "\n",
      "[1] \"Normalizing...\"\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Warning message in .filterFeatures(expr, method):\n",
      "“1 genes with constant expression values throuhgout the samples.”\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Estimating ssGSEA scores for 1 gene sets.\n",
      "[1] \"Calculating ranks...\"\n",
      "[1] \"Calculating absolute values from ranks...\"\n",
      "  |======================================================================| 100%\n",
      "\n",
      "[1] \"Normalizing...\"\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Warning message in .filterFeatures(expr, method):\n",
      "“1 genes with constant expression values throuhgout the samples.”\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Estimating ssGSEA scores for 1 gene sets.\n",
      "[1] \"Calculating ranks...\"\n",
      "[1] \"Calculating absolute values from ranks...\"\n",
      "  |==                                                                    |   4%"
     ]
    }
   ],
   "source": [
    "for (i in 1:length(names(gene_signatures))){\n",
    "    create_directory(names(gene_signatures)[i])\n",
    "    gene_list <- gene_signatures[[names(gene_signatures)[i]]]\n",
    "    geneset=data.frame(gene_list)\n",
    "    cox_os <- NULL\n",
    "    cox_pfi <- NULL\n",
    "    for (j in 1:length(names(TCGA_expr))){\n",
    "        expr <- TCGA_expr[[names(TCGA_expr)[j]]]\n",
    "        pheno <- TCGA_pheno[[names(TCGA_expr)[j]]]\n",
    "        \n",
    "        ssgsea <- GSVA::gsva(unique(as.matrix(expr)),\n",
    "                         geneset,\n",
    "                         method='ssgsea',\n",
    "                         mx.diff=TRUE,\n",
    "                         kcdf='Gaussian',\n",
    "                         parallel.sz=1)\n",
    "        df_ssgsea=as.data.frame(t(ssgsea))\n",
    "        colnames(df_ssgsea)=\"ssgsea\"\n",
    "        \n",
    "        pheno$score <- NA\n",
    "        for (sample in intersect(pheno$sample_id, rownames(df_ssgsea))){\n",
    "            pheno[which(pheno$sample_id==sample), \"score\"] <- as.numeric(df_ssgsea[sample,])\n",
    "        }\n",
    "        pheno <- pheno[which(!is.na(pheno$score)),]\n",
    "        \n",
    "        OS_KMplot_dir <- paste0(\"./TLS_project/results/siganatures_TCGA/\", names(gene_signatures)[i], \"/KMPlot/OS/\", names(TCGA_expr)[j], \".pdf\")\n",
    "        PFI_KMplot_dir <- paste0(\"./TLS_project/results/siganatures_TCGA/\", names(gene_signatures)[i], \"/KMPlot/PFI/\", names(TCGA_expr)[j], \".pdf\")\n",
    "        \n",
    "        Get_KMplot(cancer_type = names(TCGA_expr)[j], status = pheno$os.status, time = pheno$os, score = pheno$score, data = pheno, dir = OS_KMplot_dir)\n",
    "        Get_KMplot(cancer_type = names(TCGA_expr)[j], status = pheno$pfi.status, time = pheno$pfi, score = pheno$score, data = pheno, dir = PFI_KMplot_dir)\n",
    "        \n",
    "        cox_os <- rbind(cox_os, Get_HR_continous(status = pheno$os.status, time = pheno$os, score = pheno$score, data = pheno))\n",
    "        cox_pfi <- rbind(cox_pfi, Get_HR_continous(status = pheno$pfi.status, time = pheno$pfi, score = pheno$score, data = pheno))\n",
    "\n",
    "    }\n",
    "    cox_os <- cbind(names(TCGA_expr), cox_os)\n",
    "    cox_pfi <- cbind(names(TCGA_expr), cox_pfi)\n",
    "    colnames(cox_os) <- c(\"study\", \"HR\", \"SE\", \"95di_low\", \"95di_high\", \"Pval\")\n",
    "    colnames(cox_pfi) <- c(\"study\", \"HR\", \"SE\", \"95di_low\", \"95di_high\", \"Pval\")\n",
    "    save(cox_os, file = paste0(\"./TLS_project/results/siganatures_TCGA/\", names(gene_signatures)[i], \"/COX_OS.Rdata\"))\n",
    "    save(cox_pfi, file = paste0(\"./TLS_project/results/siganatures_TCGA/\", names(gene_signatures)[i], \"/COX_PFI.Rdata\"))\n",
    "}"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "R_4.2.3",
   "language": "R",
   "name": "r"
  },
  "language_info": {
   "codemirror_mode": "r",
   "file_extension": ".r",
   "mimetype": "text/x-r-source",
   "name": "R",
   "pygments_lexer": "r",
   "version": "4.2.3"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
