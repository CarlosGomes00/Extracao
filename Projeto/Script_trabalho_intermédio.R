############## Trabalho de extração  #################
#Autoria de: Carlos Gomes (PG51681), Laís Carvalho (PG52536), Rita Nóbrega (PG46733)

## Packages necessários
library(readr)
library(readxl)
library(data.table)
library(tibble)
library(car)
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(summarytools)
library(limma)
library(Glimma)
library(edgeR)

###################### Carregamento dos dados###############################
RNA_cpm = read.table("data_mrna_seq_cpm.txt", header = T, sep = '\t')
data_patient = read.table("data_clinical_patient.txt", header = T, sep = '\t')
data_sample = read_tsv("data_clinical_sample.txt", na = "")
data_mutations = read.table("data_mutations.txt", header = T, sep = '\t')

##################### Veridicação dos datasets ################################
## Dados de RNA-seq 
class(RNA_cpm)
head(unlist(lapply(RNA_cpm,class)) )
dim(RNA_cpm)
head(colnames(RNA_cpm))
head(rownames(RNA_cpm))
which(duplicated(RNA_cpm[,1]))

#Marcação dos nomes de genes que estão duplicados
RNA_cpm[which(duplicated(RNA_cpm[,1])),1]
which(RNA_cpm == 'BTBD8', arr.ind = T)
RNA_cpm[which(duplicated(RNA_cpm[,1])),1] = paste(RNA_cpm[which(duplicated(RNA_cpm[,1])),1], "_duplicado", sep = "")
which(duplicated(RNA_cpm[,1]))
RNA_cpm[13532, 1]

#Associar genes às linhas
rownames(RNA_cpm) = RNA_cpm[, 1]
head(rownames(RNA_cpm))

#NA
any(is.na(RNA_cpm))
RNA_cpm$Entrez_Gene_Id <- NULL
any(is.na(RNA_cpm))

#Verificação de genes com valores de contagens CPM 0 frequentes
verificar_zero_frequente <- function(RNA_cpm) {
  genes_zero_frequente <- character(0)
  total_amostras <- ncol(RNA_cpm) - 1 
  for (i in 1:nrow(RNA_cpm)) {
    proporcao_zero <- sum(RNA_cpm[i, -1] == 0) / total_amostras
    if (proporcao_zero > 0.9) {
      genes_zero_frequente <- c(genes_zero_frequente, RNA_cpm[i, 1])
    }
  }
  return(genes_zero_frequente)
}
verificar_zero_frequente(RNA_cpm)

## Metadados referentes às amostras
#Verificação do formato do dataset
class(data_sample) 
data_sample = as.data.frame(data_sample)
class(data_sample)
str(data_sample)

#Associação das colunas às variáveis que caracterizam as amostras
colnames(data_sample)
colnames(data_sample) = data_sample[4,]
colnames(data_sample)

#Dimensões dos metadados das amostras
dim(data_sample)

#Associação do nome das amostras às linhas e modificação destes para corresponderem ao dos dados
row.names(data_sample) = data_sample$SAMPLE_ID
row.names(data_sample)  
row.names(data_sample) = gsub("-", ".", row.names(data_sample)) 
row.names(data_sample)
data_sample = data_sample[,-2]

#Verificação de duplicados
any(duplicated(data_sample$PATIENT_ID))
any(duplicated(data_sample$SAMPLE_ID))

#Filtragem para ficarmos apenas com os metadados referentes a amostras que foram sequenciadas para RNA-seq
table(data_sample$RNA_SEQ_ANALYSIS)
pie(table(data_sample$RNA_SEQ_ANALYSIS), main = "RNA-Seq analysis?")
sample_data_rna = data_sample[data_sample$RNA_SEQ_ANALYSIS == 'Yes',]
dim(sample_data_rna)
sum(row.names(sample_data_rna) %in% colnames(RNA_cpm))

## Metadados correspondentes a informações de pacientes
class(data_patient)
dim(data_patient)
unlist(lapply(data_patient,class))
str(data_patient)
row.names(data_patient)
row.names(data_patient) <- data_patient$PATIENT_ID
row.names(data_patient)
data_patient <- data_patient[,-1]
sum(is.na(data_patient))

#Filtragem para só ficarmos com metadados relativos a pacientes cujas amostras foram usadas para RNA-seq
patient_data_rna = data_patient %>% filter(row.names(data_patient) %in% sample_data_rna$PATIENT_ID)

## Dados sobre mutações
class(data_mutations) 
dim(data_mutations)
unlist(lapply(data_mutations,class))
row.names(data_mutations)
row.names(data_mutations) = make.unique(data_mutations$Hugo_Symbol) 
sum(is.na(data_mutations))

#Definir como fatores e como variáveis numéricas
str(data_mutations)
colunas_fator <- c(1, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,29,30,31,32, 33, 45, 47, 48, 51, 57, 60, 64)
colunas_numericas <- c(34,35,36,37,43, 54, 58, 63)
colunas_char <- c(2,6,7)
data_mutations[, colunas_fator]<- lapply(data_mutations[, colunas_fator], as.factor)
data_mutations[, colunas_numericas]<- lapply(data_mutations[, colunas_numericas], as.numeric)
data_mutations[, colunas_char]<- lapply(data_mutations[, colunas_char], as.character)
str(data_mutations)

#Filtragem dos genes que estão presentes nos dados de RNA-seq
mutations_subset = data_mutations %>% filter(data_mutations$Hugo_Symbol %in% RNA_cpm$Hugo_Symbol)

########################### Sumarização dos dados ##############################################
## Metadados dos pacientes
print(dfSummary(patient_data_rna, style = 'grid', graph.magnif = 1, valid.col = FALSE,
                max.distinct.values = 5, col.widths = c(100, 200, 200, 350, 500, 250),
                dfSummary.silent  = TRUE, headings = FALSE, justify = 'l')
      , method = 'render', max.tbl.height = 500)
colunas_relevantes = c("SEX", "ETHNICITY", "AGE_AT_DIAGNOSIS", "PRIOR_CANCER", "TREATMENT_TYPE", "OS_STATUS", "OS_MONTHS", "CAUSE_OF_DEATH_SOURCE", "DIAGNOSIS")

relev_data_patient <- patient_data_rna[,colunas_relevantes]
#Status
boxplot(AGE_AT_DIAGNOSIS ~ OS_STATUS, data = relev_data_patient, main = "Idade aquando do diagnostico vs Estado", xlab = "Estado", ylab = "Idade (anos)", col = c("lightblue", "lightgreen", "grey"), names = c("Unknown", "0:LIVING", "1:DECEASED"))
#Etnias
valores = table(relev_data_patient$ETHNICITY)
nome = c("AdmixedAsian", "AdmixedBlack", "AdmixedWhite", "Asian", "Black", "HispNative", "White")
pie(valores, labels = rep("", length(valores)), col = c("black",3,4,5,6,7,2), main = "Etinias") 
legend("right", legend = nome, bty = "n", cex = 0.8, fill = c("black",3,4,5,6,7,2))
#Idade a quando do diagnóstico
boxplot(AGE_AT_DIAGNOSIS ~ SEX, data = relev_data_patient, main = "Idade do diagnostico vs gênero", xlab = "Gênero", ylab = "Idade em anos", col = c(2,'lightblue'))
#Diagnóstico
d = table(relev_data_patient$DIAGNOSIS)
d = as.data.frame(d)
indice_maximo <- which.max(d$Freq)
nome_maximo <- d$Var1[indice_maximo]
nome_maximo

## Metadados das amostras
print(dfSummary(sample_data_rna, style = 'grid', graph.magnif = 1, valid.col = FALSE,
                max.distinct.values = 5, col.widths = c(100, 200, 200, 350, 500, 250),
                dfSummary.silent  = TRUE, headings = FALSE, justify = 'l')
      , method = 'render', max.tbl.height = 500)
colunas_relevantes = c("PATIENT_ID", "GROUP", "SAMPLE_SITE", "PLATFORM", "PB_LYMPHOCYTES_PERCENTAGE", "PB_MONOCYTES_PERCENTAGE")
relev_data_sample <- sample_data_rna[,colunas_relevantes]
#local de recolha das amostras
pie(table(relev_data_sample$SAMPLE_SITE), main = "Local de recolha das amostras")
#Grupo de acordo com o diagnóstico
x = table(relev_data_sample$GROUP)
x = as.data.frame(x)
indice_maximo <- which.max(x$Freq)
nome_maximo <- x$Var1[indice_maximo]
nome_maximo

##Mutações
print(dfSummary(mutations_subset, style = 'grid', graph.magnif = 1, valid.col = FALSE,
                max.distinct.values = 5, col.widths = c(100, 200, 200, 350, 500, 250),
                dfSummary.silent  = TRUE, headings = FALSE, justify = 'l')
      , method = 'render', max.tbl.height = 500)

colunas_relevantes = c("Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification", "Variant_Type", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Validation_Status", "Mutation_Status", "impact", "t_vaf", "tumor_only")
relev_data_mutations = mutations_subset[, colunas_relevantes]
#Tipo de mutação e em que cromossomas
barplot(table(relev_data_mutations$Chromosome), main = "Cromossomos", col = "4")
pie(table(relev_data_mutations$Variant_Type), main = "Tipo de variação")
#Mutação mais comum no alelo1
mut = table(relev_data_mutations$Tumor_Seq_Allele1)
mut = as.data.frame(mut)
indice_maximo <- which.max(mut$Freq)
nome_maximo <- mut$Var1[indice_maximo]
nome_maximo 
#Mutação mais coum no alelo2
mut2 = table(relev_data_mutations$Tumor_Seq_Allele2)
mut2 = as.data.frame(mut2)
indice_maximo <- which.max(mut2$Freq)
nome_maximo <- mut2$Var1[indice_maximo]
nome_maximo

############### Análise exploratória ###################
variaveis_patient = c('PATIENT_ID', 'SEX', 'AGE_AT_DIAGNOSIS', 'PRIOR_CANCER', 'TREATMENT_TYPE')
patient_data_rna$PATIENT_ID = row.names(patient_data_rna)
variaveis_sample = c('PATIENT_ID', 'SAMPLE_ID', 'GROUP', 'SAMPLE_SITE')
sample_data_rna$SAMPLE_ID = row.names(sample_data_rna)
subset_patient <- patient_data_rna[, variaveis_patient]
subset_sample <- sample_data_rna[, variaveis_sample]
subset_patient
subset_sample

#################### Análise univariada ##########################
table(patient_data_rna$SEX) #apenas dois grupos
# Testar a variância

variancia_por_sexo <- tapply(patient_data_rna$AGE_AT_DIAGNOSIS, patient_data_rna$SEX, var, na.rm = TRUE)

variancia_por_sexo
# Testar a normalidade
teste_normalidade_masculino=shapiro.test(patient_data_rna$AGE_AT_DIAGNOSIS[patient_data_rna$SEX == "Male"])
teste_normalidade_feminino=shapiro.test(patient_data_rna$AGE_AT_DIAGNOSIS[patient_data_rna$SEX == "Female"])
print(teste_normalidade_masculino)
print(teste_normalidade_feminino)
qqPlot(patient_data_rna$AGE_AT_DIAGNOSIS[patient_data_rna$SEX == "Male"], main = "Q-Q Plot - Masculino")
qqPlot(patient_data_rna$AGE_AT_DIAGNOSIS[patient_data_rna$SEX == "Female"], main = "Q-Q Plot - Feminino")

#Ver se o sexo influencia algum campo
t.test(patient_data_rna$AGE_AT_DIAGNOSIS ~ patient_data_rna$SEX, na.rm = TRUE)
qqPlot(patient_data_rna$AGE_AT_DIAGNOSIS, main = "Q-Q Plot")
idades_sem_na <- na.omit(patient_data_rna$AGE_AT_DIAGNOSIS)
# Criar gráfico de densidade sem valores ausentes
ggdensity(idades_sem_na,
          main = "Gráfico de densidade para a idade dos pacientes",
          xlab = "Idade dos pacientes")

# Avaliar a distribuição variável
sum(is.na(mutations_subset$n_vaf))
ggdensity(mutations_subset$n_vaf)
table(mutations_subset$impact)
qqPlot(mutations_subset$n_vaf) # não é normal
leveneTest(mutations_subset$n_vaf~mutations_subset$impact)
res = kruskal.test(mutations_subset$n_vaf~mutations_subset$impact, data = mutations_subset)
res 


## Análise para as variáveis "AGE_AT_PROCUREMENT" e "ELN_2017"
table(data_sample$ELN_2017)
summary(sample_data_rna$AGE_AT_PROCUREMENT)
# Criar o gráfico de dispersão 
ggplot(sample_data_rna, aes(x = ELN_2017, y = AGE_AT_PROCUREMENT, color = ELN_2017)) +
  geom_point(size = 3) +  
  labs(x = "Classificação de Risco ELN 2017", y = "Idade na Aquisição", title = "Idade na Aquisição vs. Classificação de Risco ELN 2017") +  
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Análise a normalidade e variancias do age_at_procurement
qqPlot(sample_data_rna$AGE_AT_PROCUREMENT)
shapiro.test(sample_data_rna$AGE_AT_PROCUREMENT) 
leveneTest(AGE_AT_PROCUREMENT ~ ELN_2017, data = sample_data_rna)
anova_result <- aov(AGE_AT_PROCUREMENT ~ ELN_2017, data = sample_data_rna)
summary(anova_result)
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)