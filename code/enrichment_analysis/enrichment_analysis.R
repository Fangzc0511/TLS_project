
suppressMessages(library(GSEABase))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(enrichplot))

file = "./TLS_project/results/summary_figures/enrichment_analysis/"
if( dir.exists( file ) ) {
   unlink( file , recursive = TRUE )
}
dir.create( file )

gene_list <- read.table("./TLS_project/results/denovo_Single_Gene/response_gene.txt",sep="\t",header = F)

gene_symbols <- gene_list[,1]

gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

enrich_res_GO <- enrichGO(gene_entrez$ENTREZID, 
                       OrgDb = org.Hs.eg.db,
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.2, 
                       minGSSize = 10, 
                       maxGSSize = 5000, 
                       universe = NULL, 
                       keyType = "ENTREZID", 
                       pAdjustMethod = "BH")

pdf(paste0(file, "/gene_enrichment_analysis_GO.pdf"),width = 10, height = 5, family="Times")
print(barplot(enrich_res_GO,
              title = "GO Enrichment Analysis",
              xlab = "Gene Ontology Term",
              ylab = "Adjusted p-value",
              showCategory = 10))
dev.off()

enrich_res_KEGG <- enrichKEGG(gene = gene_entrez$ENTREZID,
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.2,
                              minGSSize = 10,
                              maxGSSize = 5000)

pdf(paste0(file, "/gene_enrichment_analysis_KEGG.pdf"),width = 10, height = 5, family="Times")
print(barplot(enrich_res_KEGG,
              title = "KEGG Enrichment Analysis",
              xlab = "Gene Ontology Term",
              ylab = "Adjusted p-value",
              showCategory = 10))
dev.off()
