# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"))

library(TCGAbiolinks)
library(SummarizedExperiment)

# pasta onde os dados serão salvos
setwd("~/Documents/workspace/miRNA-tcga-brca")
dir.create("data", showWarnings = FALSE)
OUTPUT <- "data/raw"

# como descobrir o que existe para baixar
projetos <- getGDCprojects()
head(projetos[, c("project_id", "name")])
projetos[grepl("^TCGA", projetos$project_id), c("project_id", "name")]

# como eu quero o brca, vou verificar quais categorias de dados existem e quantos casos
resumo <- getProjectSummary("TCGA-BRCA")
resumo$data_categories

# tive 10 resultados como categorias, como clinical, copy number variation, DNA Methylation, Transcriptome Profiling,...
#   como eu quero estudar miRNA irei pegar o 'Transcriptome Profiling'

# dentro de uma categoria, quais data.type existem?
#   pode demorar um pouco pois lista TODOS os arquivos de transcriptoma
q_tudo <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling")
res_tud <- getResults(q_tudo)
table(res_tudo$data_type)

# para miRNA interessam:
#   miRNA Expression Quantification -> contagem por miRNA maduro/precursor
#                                        (1 linha por miRNA; o mais usado)
#   Isoform Expression Qunatification -> contagem por isoforma (isomiR):
#                                        mais detalhado e maior (~417 MB vs ~60 MB)
# (os outros, Gene Expression Quantification etc, sao de mRNA)



# -------------------------------------------------------------
#        DOWNLOAD DOS DADOS
# -------------------------------------------------------------

# Montando a consultado (query) de miRNA:
# Principais argumentos do GDCquery():
#   project             "TCGA-BRCA"
#   data.category       "Transcriptome Profiling"
#   data.type           "miRNA Expression Quantification"
#   workflow.type       pipeline que gerou o dado para miRNA RCGSC mRNA Profiling
#   sample.type         filtra tipo de amostra: "Primary Tumor", "Solid Tissue Normal", "Metastatic"
#   barcode             vetor de barcodes para baixar só alguns pacientes
#   acess               "open" ou "controlled" (miRNA é open, não precisa token)

query_mirna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification"
)

# Inspecionando o resultado antes de baixar
res <- getResults(query_mirna)
dim(res)                    # tem 1207 arquivos
colnames(res)               # todas as informações de cada arquivo
table(res$sample_type)      # Primary Tumor: 1096 | Solid Tumor: 104 | Metastatic: 7
sum(res$file_size) / 1e6    # Tamanho total em MB (~60 MB, é pequeno)
head(res$cases)

# Para baixar:
GDCdownload(
    query  = query_mirna,
    method = "api",
    directory = "GDCdata",
    files.per.chunk = 100
)

# Juntando tudo em uma matriz
# cada arquivo .txt representa um paciente. dentro de cada arquivo txt tem um df contendo:
#   miRNA_ID:  com o id do miRNA
#   read_count
#   reads_per_million_miRNA_mapped: RPM já normalizado
#   cross-mapped: N or Y


mirna <- GDCprepare(
    query = query_mirna,
    directory = "GDCdata",
    save = TRUE,
    save.filename = paste0(OUTPUT, "/TCGA_BRCA_miRNA_raw.rda")
)

dim(mirna) #1881 x 3622
mirna[1:5, 1:4]

# separando contagens brutas (para DESeq2)
rownames(mirna) <- mirna$miRNA_ID # transforma o nome das linhas para a coluna de ID
counts <- mirna[, grep("^read_count_", colnames(mirna))]
colnames(counts) <- sub("^read_count_", "", colnames(counts))

rpm <- mirna[, grep("^reads_per_million_miRNA_mapped_", colnames(mirna))]
colnames(rpm) <- sub("^reads_per_million_miRNA_mapped_", "", colnames(rpm))

# Salva o dado de contagem bruta e normalizado
# irei salver em .CSV (se quiser pode salvar em rds(que é mais leve) é só trocar o write.csv por saveRDS)
write.csv(counts, paste0(OUTPUT, "/TCGA_BRCA_miRNA_counts.csv"))
write.csv(rpm, paste0(OUTPUT, "/TCGA_BRCA_miRNA_rpm.csv"))

# Download de dados clinicos
clinical <- GDCquery_clinic(
  project = "TCGA-BRCA",
  type = "clinical"
)

colnames(clinical)
write.csv(clinical, paste0(OUTPUT, "/TCGA_BRCA_clinical.csv"))

# os subtipos moleculares não estao incluidos no clinical data, para isso iremos importar
# Subtipos moleculares PAM50 (LumA, LumB, Her2, Basal, Normal):
subtipos <- TCGAquery_subtype(tumor = "BRCA")

# distribuicao do numero de amostras por subtipo
table(subtipos$BRCA_Subtype_PAM50)
write.csv(subtipos, paste0(OUTPUT, "/TCGA_BRCA_subtypes.csv"))