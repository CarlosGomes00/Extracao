# packages que podemos precisar (os que estão como comentário são aqueles que não tenho a certeza)
library(readr)
library(readxl)
library(data.table)
library(tibble)
library(genefilter)
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
library(stringr)
library(sparseMatrixStats)
library(DelayedMatrixStats)

# Carregamento dos dados através do url: https://raw.githubusercontent.com/lais-carvalho/Trabalho-pratico/main/Dados/

# Primeiro associamos o url a uma variável
RNA_cpm_url = "https://raw.githubusercontent.com/lais-carvalho/Trabalho-pratico/main/Dados/data_mrna_seq_cpm.txt"
data_patient_url = "https://raw.githubusercontent.com/lais-carvalho/Trabalho-pratico/main/Dados/data_clinical_patient.txt"
data_sample_url = "https://raw.githubusercontent.com/lais-carvalho/Trabalho-pratico/main/Dados/data_clinical_sample.txt"
data_mutations_url = "https://raw.githubusercontent.com/lais-carvalho/Trabalho-pratico/main/Dados/data_mutations.txt"
# Depois fazemos o download do ficheiro no url e guardamos uma cópia local 
download.file(RNA_cpm_url, destfile="RNA_seq_cpm.txt")
download.file(data_patient_url, destfile="data_patient.txt") 
download.file(data_sample_url, destfile ="data_sample.txt")
download.file(data_mutations_url, destfile="data_mutations.txt")

# Leitura do ficheiro 
RNA_cpm = read.table("data_mrna_seq.txt", header = T, sep = '\t')
data_patient = read.table("data_patient.txt", header = T, sep = '\t')
data_sample = read_tsv("data_sample.txt", na = "")
data_mutations = read.table("data_mutations.txt", header = T, sep = '\t')

## Verificar a estrutura dos dados ##
class(RNA_cpm)
# É um dataframe, o que é o suposto

## Verificação da classe de cada coluna ##
unlist(lapply(RNA_cpm,class))  
# todas têm a classe suposto (ou seja, o nome dos genes é visto como string e nas colunas relativas às amostras sao numéricos)
#a unica que está "estranha" é a 2ª coluna

## Dimensão do dataframe ##
dim(RNA_cpm)
# temos 22843 linhas e 453 colunas
# dado que os genes estão representados nas linhas e as amostras estão representadas nas colunas, isto significa que temos, no total, 22843 genes nesta análise e (dado que há um coluna vazia e a 1ª coluna são o nome dos genes) temos 251 amostras tal como esperado

## Verificar se o nome das linhas está associado aos genes e as colunas às amostras ##
colnames(RNA_cpm)  #já sabiamos isto devido ao termos verificado a class de cada coluna
#os nomes das amostras parecem bons
# vamos que ter cuidado porque nos metadados o nome das amostras não é exatamente igual a este por isso teremos que mudar para corresponder

head(rownames(RNA_cpm)) #como esperado o nome dos genes não está associado ao nome das linhas porque não o conseguimos fazer aquando da função read.table()
# há nomes duplicados por isso temos que encontrar os que estão duplicados e adicionar um '_duplicado' para sabermos que o são e conseguirmos associar o nome dos genes às linhas

which(duplicated(RNA_cpm[,1])) # há 6 genes que têm os nomes repetidos
RNA_cpm[which(duplicated(RNA_cpm[,1])),1]  #para obter o nome dos genes duplicados
which(RNA_cpm == 'BTBD8', arr.ind = T)     #para verificar que de facto o nome está duplicado
#depois de verificação manual de dois exemplos destes nomes de genes duplicados, verificámos que de facto o nome é o mesmo porém os valores de Z-cores são diferentes
#isto faz com que tenhamos de os distinguir adicionando no final do duplicado '_duplicado'
RNA_cpm[which(duplicated(RNA_cpm[,1])),1] = paste(RNA_cpm[which(duplicated(RNA_cpm[,1])),1], "_duplicado", sep = "")
which(duplicated(RNA_cpm[,1]))
RNA_cpm[13532, 1]
#os nomes duplicados foram tratados por isso já podemos associar o nome dos genes às linhas
rownames(RNA_cpm) = RNA_cpm[, 1]
rownames(RNA_cpm)

## Verificar se existem missing values ##
any(is.na(RNA_cpm))
# Suspeitamos que os valores da coluna 'Entrez_gene_id' são os únicos valores únicos por isso vamos eliminar esta coluna do dataframe
RNA_cpm$Entrez_Gene_Id <- NULL
# Verificar outra vez se ainda há missing values
any(is.na(RNA_cpm))
# Já não existem NA no nosso dataset

## Pré-processamento dos dados feito anteriormente ##
# O nosso dataset consiste em valores CPM (counts per million), uma medida usada para normalizar a expressão génica e facilitar a comparação entre amostras com nº total de reads diferentes
# Valor CPM = Raw count para 1 gene / 1 milhão / número total de reads mapeadas na amostra
# Ou seja, quanto maior for o CPM, maior o na expressão de determinado gene numa determinada amostra (relativamente a outros genes na mesma amostra)
# antes desta conversão para CPM, os dados foram normalizados usando o 'Trimmed Mean of M caling' que serve para corrigir diferenças sistemáticas na composição das bibliotecas de RNA entre amostras (não é uma normalização dos dados em si, mas sim do tamanho da biblioteca)
# esta normalização TMM já ajuda a reduzir a variabilidade técnica entre amostras, mas os dados ainda podem ter uma distribuição não normal por isso podemos aplicar uma transformação logaritmica

## ver se há outliers ??
## aplicar filtros de flat pattern para retirar genes com baixa expressão e variabilidade não significativa
verificar_zero_frequente <- function(RNA_cpm) {
  # Inicializa um vetor para armazenar os nomes dos genes com expressão zero frequente
  genes_zero_frequente <- character(0)
  total_amostras <- ncol(RNA_cpm) - 1 
  # Itera sobre as linhas do dataframe
  for (i in 1:nrow(RNA_cpm)) {
    # Calcula a proporção de valores zero na linha atual
    proporcao_zero <- sum(RNA_cpm[i, -1] == 0) / total_amostras
    # Se a proporção de valores zero for maior que 90%, adiciona o nome do gene à lista
    if (proporcao_zero > 0.9) {
      genes_zero_frequente <- c(genes_zero_frequente, RNA_cpm[i, 1])
    }
  }
  # Retorna os nomes dos genes com expressão zero frequente
  return(genes_zero_frequente)
}

verificar_zero_frequente(RNA_cpm)  # não temos dados que tenham contagens com mais de 90% de 0

## supostamente o edgeR faz o TMM por isso teoricamente poderiamos utilizar esse package para a DE?

############### Metadados #########################
#Verificar a classe do metadado:
class(data_patient) # É um data frame, como queremos.
dim(data_patient)
unlist(lapply(data_patient,class))
str(data_patient)
row.names(data_patient)
row.names(data_patient) <- data_patient$PATIENT_ID
row.names(data_patient)
data_patient <- data_patient[,-1]

#Não sei se era apropriado, fiz de outra forma
#data_patient[, -1] <- lapply(data_patient[, -1], as.factor) # não sei se é correto aplicar a todas as colunas
# Decidi não selecionar nada, acho que talvez selelcionar a medida que analisamos:
## Transformação das variáveis:
str(data_patient)
data_patient[, c(1,2,4,5,7, 8, 9, 10, 11, 12, 13, 14, 18, 20, 21)] <- lapply(data_patient[, c(1,2,4,5,7, 8, 9, 10, 11, 12, 13, 14, 18, 20, 21)], as.factor)

sum(is.na(data_patient)) # há NA mas é normal em metadados e não iremos fazer nada
str(data_patient)

#Data_sample:
class(data_sample) # Há vários formatos, um deles é de dataframe, como queremos
data_sample = as.data.frame(data_sample) #Se quiser estabelecer em um formato só.
class(data_sample)
str(data_sample)

colnames(data_sample) # para facilitar a nossa análise vamos pôr o nome das colunas a corresponder ao nome que está nos outros datasets
colnames(data_sample) = data_sample[4,]
colnames(data_sample)
data_sample = data_sample[-c(1:4),] # Retirando as duas primeiras linhas que são apenas explicações das variáveis 
dim(data_sample) # aproximadamente 73 variáveis e 672 amostras como esperado

row.names(data_sample)  # as linhas não estão associadas a nada, como isto se trata de metadados referentes às amostras vamos associar cada linha ao ID de uma amostra (já que tambem são as amostras que estão nos nossos dados)
row.names(data_sample) = data_sample$SAMPLE_ID
row.names(data_sample)  # ainda é necessário certificarmo-nos que o ID da amostra corresponde ao ID que estão nos dados de RNA seq por isso vamos transformar o '-' em '.'
row.names(data_sample) = gsub("-", ".", row.names(data_sample)) 
row.names(data_sample)
data_sample = data_sample[,-2]

#data_sample[, -c(1, 2)] <- lapply(data_sample[, -c(1, 2)], as.factor)  #não sei se isto é o mais correto
# talvez aqui seja melhor descrever todas as váriaveis (colunas) e já identificar quais pomos como fatores e que ireos usar a seguir
str(data_sample)
colunas_fator <- c(2, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 21, 23, 24, 26, 28, 34, 45, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 70, 72)
data_sample[,colunas_fator] <- lapply(data_sample[,colunas_fator], as.factor)
colunas_numericas <- c(25, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 46, 47, 49, 50, 52, 54, 55, 69, 74)
data_sample[,colunas_numericas] <- lapply(data_sample[,colunas_numericas], as.numeric)
str(data_sample)

# Os dados estão em relação as amostras, então o ID dos pacientes podem aparecer duplicados
any(duplicated(data_sample$PATIENT_ID)) #Id de paciente duplicados, mas é normal dado que do mesmo paciente foram tiradas mais do que uma amostra
# Não é suposto ter ID de amostra duplicado
any(duplicated(data_sample$SAMPLE_ID)) #confirmado que não há

###### Dados sobre mutações nos genes que serão usados posteriormente à análise da expressão diferencial ####
class(data_mutations) # É um data frame
dim(data_mutations)
unlist(lapply(data_mutations,class))
row.names(data_mutations)
row.names(data_mutations) = data_mutations$Hugo_Symbol  # não dá para fazer porque temos genes com o mesmo nome, mas aqui não podemos por 'duplicado' no final então vamos deixar assim
sum(is.na(data_mutations))
## Preparação dos dados:
str(data_mutations)
colunas_fator <- c(1, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,29,30,31,32, 33, 45, 47, 48, 51, 57, 60, 64)
colunas_numericas <- c(34,35,36,37,43, 54, 58, 63)
colunas_char <- c(2,6,7)
data_mutations[, colunas_fator]<- lapply(data_mutations[, colunas_fator], as.factor)
data_mutations[, colunas_numericas]<- lapply(data_mutations[, colunas_numericas], as.numeric)
data_mutations[, colunas_char]<- lapply(data_mutations[, colunas_char], as.character)
str(data_mutations)

#Nem todas as amostras foram utilizadas para o processo de análise de RNAseq
# Há duas variáveis nos metadados: RNA sequenced e RNA seq analysis
# Através da tabela de dados de RNA-seq sabemos que 451 amostras foram analisadas
table(data_sample$RNA_SEQ_ANALYSIS) #corresponde ao esperado por isso vamos criar um subset de metadados apenas com estas amostras para facilitar análises posteriores
pie(table(data_sample$RNA_SEQ_ANALYSIS), main = "RNA-Seq analysis?")
sample_data_rna = data_sample[data_sample$RNA_SEQ_ANALYSIS == 'Yes',]
dim(sample_data_rna) # ficámos com 451 linhas, ou seja, metadados relativos a 451 amostras cujo RNA foi sequenciado e analisado
sum(row.names(sample_data_rna) %in% colnames(RNA_cpm)) # filtragem bem sucedida pois temos 451 correspondências

#fazer a mesma filtragem (criação do subset) para os metadados dos pacientes
patient_data_rna = data_patient %>% filter(data_patient$PATIENT_ID %in% sample_data_rna$PATIENT_ID)

#para análises posterior vamos querer ver se há alguma assoicação entre mutações nos genes e a sua expressão daí fazermos a filtragem 
mutations_subset = data_mutations %>% filter(data_mutations$Hugo_Symbol %in% RNA_cpm$Hugo_Symbol)

############# Sumarização dos dados ##############
library(summarytools)
## data patient:
print(dfSummary(data_patient, style = 'grid', graph.magnif = 1, valid.col = FALSE,
                max.distinct.values = 5, col.widths = c(100, 200, 200, 350, 500, 250),
                dfSummary.silent  = TRUE, headings = FALSE, justify = 'l')
      , method = 'render', max.tbl.height = 500)

# Considerando variáveis relevantes:
colunas_relevantes = c("SEX", "ETHNICITY", "AGE_AT_DIAGNOSIS", "PRIOR_CANCER", "TREATMENT_TYPE", "OS_STATUS", "OS_MONTHS", "CAUSE_OF_DEATH_SOURCE", "DIAGNOSIS")
relev_data_patient <- data_patient[,colunas_relevantes]
boxplot(AGE_AT_DIAGNOSIS ~ OS_STATUS, data = relev_data_patient, main = "Idade do diagnostico vs status", xlab = "Status", ylab = "Idade em anos", col = c("blue", "green", "red"), names = c("NA", "0:LIVING", "1:DECEASED"))

valores = table(relev_data_patient$ETHNICITY)
nome = c("AdmixedAsian", "AdmixedBlack", "AdmixedWhite", "Asian", "Black", "HispNative", "White")

pie(valores, labels = rep("", length(valores)), col = c("black",3,4,5,6,7,2), main = "Etinias") 
legend("right", legend = nome, bty = "n", cex = 0.8, fill = c("black",3,4,5,6,7,2))

boxplot(AGE_AT_DIAGNOSIS ~ SEX, data = relev_data_patient, main = "Idade do diagnostico vs gênero", xlab = "Gênero", ylab = "Idade em anos", col = c(2,5))

d = table(relev_data_patient$DIAGNOSIS)
d = as.data.frame(d)
indice_maximo <- which.max(d$Freq)
nome_maximo <- d$Var1[indice_maximo]
nome_maximo #Diagnóstico mais recorrente

## data_sample
print(dfSummary(data_sample, style = 'grid', graph.magnif = 1, valid.col = FALSE,
                max.distinct.values = 5, col.widths = c(100, 200, 200, 350, 500, 250),
                dfSummary.silent  = TRUE, headings = FALSE, justify = 'l')
      , method = 'render', max.tbl.height = 500)

colunas_relevantes = c("PATIENT_ID", "GROUP", "SAMPLE_SITE", "PLATFORM", "PB_LYMPHOCYTES_PERCENTAGE", "PB_MONOCYTES_PERCENTAGE")
relev_data_sample <- data_sample[,colunas_relevantes]
pie(table(relev_data_sample$SAMPLE_SITE), main = "Local de recolha das amostras")

x = table(relev_data_sample$GROUP)
x = as.data.frame(x)
indice_maximo <- which.max(x$Freq)
nome_maximo <- x$Var1[indice_maximo]
nome_maximo ## grupo do diagnóstico mais analisado.

## Data_mutations
print(dfSummary(data_mutations, style = 'grid', graph.magnif = 1, valid.col = FALSE,
                max.distinct.values = 5, col.widths = c(100, 200, 200, 350, 500, 250),
                dfSummary.silent  = TRUE, headings = FALSE, justify = 'l')
      , method = 'render', max.tbl.height = 500)

colunas_relevantes = c("Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification", "Variant_Type", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Validation_Status", "Mutation_Status", "impact", "t_vaf", "tumor_only")
relev_data_mutations = data_mutations[, colunas_relevantes]

barplot(table(relev_data_mutations$Chromosome), main = "Cromossomos", col = "4")

pie(table(relev_data_mutation$Variant_Type), main = "Tipo de variação")

mut = table(relev_data_mutations$Tumor_Seq_Allele1)
mut = as.data.frame(mut)
indice_maximo <- which.max(mut$Freq)
nome_maximo <- mut$Var1[indice_maximo]
nome_maximo # A mutação mais comun é a alteração por C

mut2 = table(relev_data_mutations$Tumor_Seq_Allele2)
mut2 = as.data.frame(mut2)
indice_maximo <- which.max(mut2$Freq)
nome_maximo <- mut2$Var1[indice_maximo]
nome_maximo # A mutação mais comum é a alteração por T

## Filtrar para ficarmos apenas com os metadados que nos interessam para posterior análise
# aqui acho que temos que selecionar os fatores que queremos estudar 
# diria para escolher sexo, treatment_type, group, sample_site, age_at_diagnosis

variaveis_patient = c('PATIENT_ID', 'SEX', 'AGE_AT_DIAGNOSIS', 'PRIOR_CANCER', 'TREATMENT_TYPE')
variaveis_sample = c('PATIENT_ID', 'SAMPLE_ID', 'GROUP', 'SAMPLE_SITE')
subset_patient <- patient_data_rna[, variaveis_patient]
subset_sample <- sample_data_rna[, variaveis_sample]

# Agora os datasets com que vamos trabalhar são 'RNA_cpm' (dados de expressão génica), 'subset_patient' e 'subset_sample' (os dois metadados referentes às amostras cujo RNA foi analisado) e 'mutation_subset' (dados que serão usados para análises posteriores)

##  ACHO QUE PODEMOS PASSAR À ANÁLISE EXPLORATÓRIA DOS DADOS ##
# aqui acho que temos que selecionar os fatores que queremos estudar e fazer as.factor
# diria para escolher sexo, treatment_type, group, sample_site, age_at_diagnosis 





