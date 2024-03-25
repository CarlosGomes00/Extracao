# packages que podemos precisar (os que estão como comentário são aqueles que não tenho a certeza)
library(readr)
library(readxl)
library(data.table)
library(tibble)
#library(tidyverse)
library(arules)
library(dplyr)
#library(summarytools)
#library(ggpubr)
library(limma)
library(Glimma)
library(edgeR)
library(org.Hs.eg.db)
library(GO.db)
#library(car)
#library(caret)
library(RColorBrewer)
library(gplots)
#library(FactoMineR)
#library(factoextra)

# Carregamento dos dados através do url: https://raw.githubusercontent.com/lais-carvalho/Trabalho-pratico/main/Dados/

# Primeiro associamos o url a uma variável
RNA_rpkm_zscore_url = "https://raw.githubusercontent.com/lais-carvalho/Trabalho-pratico/main/Dados/data_mrna_seq_rpkm_zscores_ref_all_samples.txt"
data_patient_url = "https://raw.githubusercontent.com/lais-carvalho/Trabalho-pratico/main/Dados/data_clinical_patient.txt"
data_sample_url = "https://raw.githubusercontent.com/lais-carvalho/Trabalho-pratico/main/Dados/data_clinical_sample.txt"

# Depois fazemos o download do ficheiro no url e guardamos uma cópia local 
download.file(RNA_rpkm_zscore_url, destfile="RNA_seq_rpkm_zscore.txt")
download.file(data_patient_url, destfile="data_patient.txt") 
download.file(data_sample_url, destfile ="data_sample.txt")

# Leitura do ficheiro 
rna_rpkm_zscore = read.table("RNA_seq_rpkm_zscore.txt", header = T, sep = '\t')
data_patient = read.table("data_patient.txt", header = T, sep = '\t')
data_sample = read_tsv("data_sample.txt", na = "")

# Verificar a estrutura dos dados
class(rna_rpkm_zscore)
## É um dataframe, o que é o suposto

# Verificação da classe de cada coluna 
unlist(lapply(rna_rpkm_zscore,class))  
## todas têm a classe suposto (ou seja, o nome dos genes é visto como string e os z-cores sao numéricos)

# Dimensão do dataframe
dim(rna_rpkm_zscore)
# temos 22843 linhas e 453 colunas
# dado que os genes estão representados nas linhas e as amostras estão representadas nas colunas, isto significa que temos 22842 genes nesta análise e (dado que há um coluna vazia) temos 251 amostras tal como esperado
# Verificar se as colunas são as corretas
rownames(rna_rpkm_zscore) #como esperado o nome dos genes não está associado ao nome das linhas porque não o conseguimos fazer aquando da função read.table()
# há nomes duplicados por isso temos que encontrar os que estão duplicados e adicionar um '_duplicado' para sabermos que o são e conseguirmos associar o nome dos genes às linhas

which(duplicated(rna_rpkm_zscore[,1])) # há 6 genes que têm os nomes repetidos
rna_rpkm_zscore[which(duplicated(rna_rpkm_zscore[,1])),1]
which(rna_rpkm_zscore == 'BTBD8', arr.ind = T)
which(rna_rpkm_zscore == 'COMMD9', arr.ind = T)
#depois de verificação manual de dois exemplos destes nomes de genes duplicados, verificámos que de facto o nome é o mesmo porém os valores de Z-cores são diferentes
#isto faz com que tenhamos de os distinguir adicionando no final do duplicado '_duplicado'
rna_rpkm_zscore[which(duplicated(rna_rpkm_zscore[,1])),1] <- paste(rna_rpkm_zscore[which(duplicated(rna_rpkm_zscore[,1])),1], "_duplicado", sep = "")
which(duplicated(rna_rpkm_zscore[,1]))
rna_rpkm_zscore[13532, 1]
#os nomes duplicados foram tratados por isso já podemos associar o nome dos genes às linhas
rownames(rna_rpkm_zscore) = rna_rpkm_zscore[, 1]
rownames(rna_rpkm_zscore)

# Verificar se existem missing values
any(is.na(rna_rpkm_zscore))
# Suspeitamos que os valores da coluna 'Entrez_gene_id' são os únicos valores únicos por isso vamos eliminar esta coluna do dataframe
rna_rpkm_zscore$Entrez_Gene_Id <- NULL
# Verificar outra vez se ainda há missing values
sum(is.na(rna_rpkm_zscore)) # ainda há 8118 valores omissos no nosso dataset
# Os valores omissos (NA) não nos vão permitir realizar análises posteriores com estes dados logo temos que tratar estes NA de alguma forma
# Provavelmente a melhor opção é substituir os NA pela média de cada coluna dado que uma coluna é uma amostra e uma linha um gene 
#( acho que o que faz mais sentido é substituir pelo valor média da amostra e não do gene porque amostras diferentes podem comportar-se de maneira diferente consoante o gene)
for (i in 2:ncol(rna_rpkm_zscore)){
  m = mean(rna_rpkm_zscore[, i], na.rm = T)
  rna_rpkm_zscore[is.na(rna_rpkm_zscore[,i]), i] = m
}

sum(is.na(rna_rpkm_zscore))
# Já não existem NA no nosso dataset


############### Metadados #########################
# falta fazer esta verificação inicial; ver se são dataframe (aquele que não é, precisamos de transformar em dataframe)
#ver se o nº de linhas corrresponde ao nº de amostras descritas
#ver se o nome das variaveis e amostras esão associados às colunas e linhas, respetivamente
#eu acho que é preciso ver se há NAs porém não os devemos substituir ou eliminar porque maior parte deles são caracteristicas das amostras logo não faz sentido eliminarou substituir pela média
# quanto ao juntar os datasets estive a ver trabalhos de outros anos e eles criam subsets dos metadados que têm as amostras presentes nos dados
#provavelmente é isto que temos de fazer

