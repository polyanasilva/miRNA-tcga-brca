# Desenvolvimento de modelo prognóstico baseado em miRNAs e aprendizado de máquina no câncer de mama

O objetivo desse projeto é iniciar o desenvolvimento de uma abordagem interpretável de aprendizado de máquina baseada em dados de expressão de miRNAs para a predição do prognóstico de pacientes com câncer de mama.

Para isso irei começar fazendo o download dos dados de miRNA do banco de dados do TCGA do projeto TCGA-BRCA, via TCGAbiolinks.

O plano inicial é explorar a expressão dos miRNA e buscar padrões entre os subtipos classificados pelo próprio TCGA.

Até o momento foi utilizado o script "download_miRNA_BRCA.R" para baixar o "miRNA Expression Quantification" e o dado clinico do projeto "TCGA-BRCA". Foi salvo em data/ o dado de contagem bruta, o rpm e o dado clinico. 