############## Trabalho de extração  #################
#Autoria de: Carlos Gomes (PG51681), Laís Carvalho (PG52536), Rita Nóbrega (PG46733)

## Packages necessários
library(summarytools)
library(car)
library(readr)
library(dplyr)
library(ggplot2)
library(gplots)
library(ggpubr)
library(genefilter)
library(edgeR)
library(limma)
library(ggfortify)
library(scatterplot3d)
library(factoextra)
library(MASS)
library(caret)
library(randomForest)
library(party)
library(e1071)

###################### Carregamento dos dados###############################
RNA_cpm = read.table("Dados/data_mrna_seq_cpm.txt", header = T, sep = '\t')   
data_patient = read.table("Dados/data_clinical_patient.txt", header = T, sep = '\t')
data_sample = read_tsv("Dados/data_clinical_sample.txt", na = "", show_col_types = FALSE)
data_mutations = read.table("Dados/data_mutations.txt", header = T, sep = '\t')

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

#Tipo de dados nos metadados das amostras
fatores_sample_data = c(3:25, 27:29, 46, 49, 52, 54, 57:69, 71:74)
for (coluna in fatores_sample_data) {
  data_sample[[coluna]] = as.factor(data_sample[[coluna]])
}
colunas_numericas_sample_data = c(26, 30:45, 47, 48, 50, 51, 53, 55, 56, 70, 75)
for (coluna in colunas_numericas_sample_data) {
  data_sample[[coluna]] = as.numeric(data_sample[[coluna]])
}
unlist(lapply(data_sample,class))

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
head(row.names(data_patient))
row.names(data_patient) <- data_patient$PATIENT_ID
head(row.names(data_patient))

#Tipo de dados nos metadados dos pacientes
fatores_patient = c(3, 2, 5:19, 21:22)
for (coluna in fatores_patient) {
  data_patient[[coluna]] = as.factor(data_patient[[coluna]])
}
colunas_numericas_patient = c(4,20)
for (coluna in colunas_numericas_patient) {
  data_patient[[coluna]] = as.numeric(data_patient[[coluna]])
}
unlist(lapply(data_patient,class))

#Valores omissos nos metadados dos pacientes
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

#Merge dos metadados de pacientes e amostras num único dataframe
metadados = merge(sample_data_rna, patient_data_rna, by='PATIENT_ID', all=T)

#Tipo das vaiáveis e seleção das colunas mais interessantes
row.names(metadados) = metadados$SAMPLE_ID
variaveis_interessantes = c(1:4, 6, 8:10, 18:29, 52, 70, 73, 76:78, 90:95)
metadados = metadados[, variaveis_interessantes]

## Sumarização dos metadados
print(dfSummary(metadados, style = 'grid', graph.magnif = 1, valid.col = FALSE,
                max.distinct.values = 5, col.widths = c(100, 200, 200, 350, 500, 250),
                dfSummary.silent  = TRUE, headings = FALSE, justify = 'l')
      , method = 'render', max.tbl.height = 500)

##Criação do subset das Mutações
colunas_relevantes = c("Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification", "Variant_Type", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Validation_Status", "Mutation_Status", "impact", "n_vaf")
relev_data_mutations = mutations_subset[, colunas_relevantes]
#Sumarização do subset das mutações
print(dfSummary(mutations_subset, style = 'grid', graph.magnif = 1, valid.col = FALSE,
                max.distinct.values = 5, col.widths = c(100, 200, 200, 350, 500, 250),
                dfSummary.silent  = TRUE, headings = FALSE, justify = 'l')
      , method = 'render', max.tbl.height = 500)


## Análise exploratória
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

#local de recolha das amostras
pie(table(relev_data_sample$SAMPLE_SITE), main = "Local de recolha das amostras")
#Grupo de acordo com o diagnóstico
x = table(relev_data_sample$GROUP)
x = as.data.frame(x)
indice_maximo <- which.max(x$Freq)
nome_maximo <- x$Var1[indice_maximo]
nome_maximo

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

#################### Análise univariada ##########################
table(metadados$SEX) # Confirmação de quantos grupos de sexos há.

# Testar a homogeneidade das variâncias
leveneTest(metadados$AGE_AT_DIAGNOSIS~metadados$SEX) # Não rejeitar a hipótese de que as variâncias dos grupos são iguais. 

# Testar a normalidade dos diferentes grupos
teste_normalidade_masculino=shapiro.test(metadados$AGE_AT_DIAGNOSIS[metadados$SEX == "Male"])
teste_normalidade_feminino=shapiro.test(metadados$AGE_AT_DIAGNOSIS[metadados$SEX == "Female"])
print(teste_normalidade_masculino) # Não normal
print(teste_normalidade_feminino) # Não normal

# Visualização da normalidade com qqPlot
qqPlot(metadados$AGE_AT_DIAGNOSIS[metadados$SEX == "Male"], main = "Q-Q Plot - Masculino", ylab = "Idade")
qqPlot(metadados$AGE_AT_DIAGNOSIS[metadados$SEX == "Female"], main = "Q-Q Plot - Feminino", ylab = "Idade")-

# Normalidade geral das idades
shapiro.test(metadados$AGE_AT_DIAGNOSIS) # Não normal

# Retirada dos NAs:
idades_sem_na <- na.omit(metadados$AGE_AT_DIAGNOSIS)

# Criar gráfico de densidade sem valores ausentes
ggdensity(idades_sem_na,
          main = "Gráfico de densidade para a idade dos pacientes",
          xlab = "Idade dos pacientes")

t.test(metadados$AGE_AT_DIAGNOSIS ~ metadados$SEX, na.rm = TRUE)

# Avaliar a variável "N_vaf" quanto a presença de NAs e distribuição
sum(is.na(relev_data_mutations$n_vaf))
ggdensity(relev_data_mutations$n_vaf, main = "Distribuição da variável `n_vaf`")

# Avaliação das diferentes classes de impacto
table(relev_data_mutations$impact) # Há 3 classes: "High", "Moderate" e "Modifier"

# Avaliar a distribuição variável
sum(is.na(mutations_subset$n_vaf))
ggdensity(mutations_subset$n_vaf)
table(mutations_subset$impact)
qqPlot(mutations_subset$n_vaf) # não é normal
leveneTest(mutations_subset$n_vaf~mutations_subset$impact)
res = kruskal.test(mutations_subset$n_vaf~mutations_subset$impact, data = mutations_subset)
res 

# Avaliação das variáveis
table(metadados$ELN_2017) # Há 6 grupos de clssificações de risco
summary(metadados$AGE_AT_PROCUREMENT) # Há um NA.

# Criar o gráfico de dispersão 
ggplot(metadados, aes(x = ELN_2017, y = AGE_AT_PROCUREMENT, color = ELN_2017)) +
  geom_point(size = 3) +  
  labs(x = "Classificação de Risco ELN 2017", y = "Idade na Aquisição", title = "Idade na Aquisição vs. Classificação de Risco ELN 2017") +  
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Analisando a normalidade e homogeneidade das variâncias:

qqPlot(metadados$AGE_AT_PROCUREMENT, ylab = "idade", main = "Distribuição das idades dos pacientes na recolha das amostras")

shapiro.test(metadados$AGE_AT_PROCUREMENT) # Não é normal

# Análise da variância
leveneTest(AGE_AT_PROCUREMENT ~ ELN_2017, data = metadados) # Homogénea

# ANOVA
anova_result <- aov(AGE_AT_PROCUREMENT ~ ELN_2017, data = metadados)
summary(anova_result)

## Análise por grupo:
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)


##################### Expressão diferencial ###############################
#Preparação dos dados para DE - filtragem
selected_genes = rowMeans(RNA_cpm[,-1] > 1) > 0.5
RNA_cpm_filtered = RNA_cpm[selected_genes, ]
RNA_cpm_filtered = RNA_cpm_filtered[, -1]

logCPM = log2(RNA_cpm_filtered + 1)
head(logCPM)

y <- DGEList(counts = logCPM, genes = row.names(logCPM))

##Expressão diferencial entre sexos
table(metadados$SEX)
design_sex = model.matrix(~0+metadados$SEX, data = y$samples)
colnames(design_sex) = levels(metadados$SEX)
row.names(design_sex) = colnames(logCPM)

contrast_sex = makeContrasts(Female - Male, levels = design_sex)
fit_sex = lmFit(logCPM, design_sex)
fit_sex = contrasts.fit(fit_sex, contrast_sex)
fit_sex = eBayes(fit_sex, trend=T)

topTable(fit = fit_sex, genelist = row.names(logCPM))

summa.fit_sex = decideTests(fit_sex)
summary(summa.fit_sex)

plotMD(fit_sex)

##DE de acordo com o tipo de amostra
levels(metadados$SAMPLE_SITE) = make.names(levels(metadados$SAMPLE_SITE))

design_site = model.matrix(~0+metadados$SAMPLE_SITE, data = y$samples)
colnames(design_site) = levels(metadados$SAMPLE_SITE)
head(design_site)

contrast_site = makeContrasts( Bone_Leuk = Bone.Marrow.Aspirate - Leukapheresis, 
                           Bone_Blood = Bone.Marrow.Aspirate - Peripheral.Blood,
                           Leuk_Blood = Leukapheresis - Peripheral.Blood,
                           levels=design_site)
fit_site = lmFit(logCPM, design_site)
fit_site = contrasts.fit(fit_site, contrast_site)
fit_site = eBayes(fit_site, trend=T)

topTable(fit = fit_site, genelist = row.names(logCPM))

summa.fit_site = decideTests(fit_site)
summary(summa.fit_site)

par(mfrow= c(1, 2))
plotMD(fit_site,coef="Leuk_Blood",status=summa.fit_site[,"Leuk_Blood"], values = c(-1, 1), hl.col=c("blue","red"), main = "Leukapheresis vs Periphal Blood")
plotMD(fit_site,coef="Bone_Leuk",status=summa.fit_site[,"Bone_Leuk"], values = c(-1, 1), hl.col=c("blue","red"), main = "Leukapheresis vs Bone Marrow")
volcano plots - site}
par(mfrow= c(1, 2))
volcanoplot(fit_site,coef="Leuk_Blood",highlight=3,names=row.names(fit_site), main ="Leukapheresis vs Periphal Blood")
volcanoplot(fit_site,coef="Bone_Leuk",highlight=3,names=row.names(fit_site), main ="Leukapheresis vs Bone Marrow")

##DE de acordo com a sobrevivencia dos pacientes
table(metadados$OS_STATUS)
metadados$OS_STATUS <- gsub("0:LIVING", "Alive", metadados$OS_STATUS)
metadados$OS_STATUS <- gsub("1:DECEASED", "Deceased", metadados$OS_STATUS)
metadados$OS_STATUS <- as.character(metadados$OS_STATUS)
metadados$OS_STATUS <- ifelse(is.na(metadados$OS_STATUS) | nchar(metadados$OS_STATUS) == 0, "Unknown", metadados$OS_STATUS)
table(metadados$OS_STATUS)
metadados$OS_STATUS <- as.factor(metadados$OS_STATUS)

design_status = model.matrix(~0+metadados$OS_STATUS, data = y$samples)
colnames(design_status) = levels(metadados$OS_STATUS)


contrast_status = makeContrasts(Alive_Deceased = Alive - Deceased,
                                Alive_Unknown = Alive - Unknown,
                                Deceased_Unknown = Deceased - Unknown,
                                levels=design_status)

fit_status = lmFit(y$counts, design_status)
fit_status = contrasts.fit(fit_status, contrast_status)
fit_status = eBayes(fit_status, trend=T)
topTable(fit = fit_status, genelist = y$genes)

summa.fit_status = decideTests(fit_status)
summary(summa.fit_status)

##DE de acordo com o estado dos tratamentos cumulativos
table(metadados$CUMULATIVE_TREATMENT_STAGE_COUNT)

metadados$CUMULATIVE_TREATMENT_STAGE_COUNT <- gsub(0, "Zero", metadados$CUMULATIVE_TREATMENT_STAGE_COUNT)
metadados$CUMULATIVE_TREATMENT_STAGE_COUNT <- gsub(1, "Um", metadados$CUMULATIVE_TREATMENT_STAGE_COUNT)
metadados$CUMULATIVE_TREATMENT_STAGE_COUNT <- gsub(2, "Dois", metadados$CUMULATIVE_TREATMENT_STAGE_COUNT)
metadados$CUMULATIVE_TREATMENT_STAGE_COUNT <- gsub(3, "Três", metadados$CUMULATIVE_TREATMENT_STAGE_COUNT)
metadados$CUMULATIVE_TREATMENT_STAGE_COUNT <- gsub(4, "Quatro", metadados$CUMULATIVE_TREATMENT_STAGE_COUNT)
metadados$CUMULATIVE_TREATMENT_STAGE_COUNT <- gsub(5, "Cinco", metadados$CUMULATIVE_TREATMENT_STAGE_COUNT)
metadados$CUMULATIVE_TREATMENT_STAGE_COUNT <- gsub(6, "Seis", metadados$CUMULATIVE_TREATMENT_STAGE_COUNT)
metadados$CUMULATIVE_TREATMENT_STAGE_COUNT <- gsub(7, "Sete", metadados$CUMULATIVE_TREATMENT_STAGE_COUNT)
metadados$CUMULATIVE_TREATMENT_STAGE_COUNT <- gsub(8, "Oito", metadados$CUMULATIVE_TREATMENT_STAGE_COUNT)

metadados$CUMULATIVE_TREATMENT_STAGE_COUNT = as.factor(metadados$CUMULATIVE_TREATMENT_STAGE_COUNT)
table(metadados$CUMULATIVE_TREATMENT_STAGE_COUNT)

design_stage = model.matrix(~0+metadados$CUMULATIVE_TREATMENT_STAGE_COUNT, data = y$samples)
colnames(design_stage) = levels(metadados$CUMULATIVE_TREATMENT_STAGE_COUNT)
head(design_stage)

contrast_stage = makeContrasts('0_1'= Zero - Um, '0_2'= Zero - Dois, '0_3'= Zero - Três, '0_4'= Zero - Quatro, '0_5'= Zero - Cinco, '0_6' = Zero - Seis, '0_7' = Zero - Sete, '0_8' = Zero - Oito, '1_2'= Um - Dois, '1_3'= Um - Três, '1_4'= Um - Quatro, '1_5'= Um - Cinco,
'1_6' = Um - Seis, '1_7' = Um - Sete, '1_8' = Um - Oito, '2_3'= Dois - Três, '2_4'= Dois - Quatro, '2_5'= Dois - Cinco, '2_6' = Dois - Seis, '2_7' = Dois - Sete, '2_8' = Dois - Oito, '3_4'= Três - Quatro, '3_5'= Três - Cinco, '3_6' = Três - Seis, '3_7' = Três - Sete, '3_8' = Três - Oito, '4_5'= Quatro - Cinco, '4_6' = Quatro - Seis, '4_7' = Quatro - Sete, '4_8' = Quatro - Oito, '5_6' = Cinco - Seis, '5_7' = Cinco - Sete, '5_8' = Cinco - Oito,
'6_7' = Seis - Sete, '6_8' = Seis - Oito, '7_8' = Sete - Oito, levels=design_stage)

fit_stage = lmFit(y$counts, design_stage)
fit_stage = contrasts.fit(fit_stage, contrast_stage)
fit_stage = eBayes(fit_stage, trend=T)
topTable(fit = fit_stage, genelist = y$genes)
summa.fit_stage = decideTests(fit_stage)
summary(summa.fit_stage)

par(mfrow=c(2, 5))
volcanoplot(fit_stage,coef="2_3",highlight=5, names=row.names(fit_stage), main ="Fases de tratamento: 2 vs 3")
volcanoplot(fit_stage,coef="3_4",highlight=1, names=row.names(fit_stage), main ="Fases de tratamento: 3 vs 4")
volcanoplot(fit_stage,coef="3_7",highlight=3, names=row.names(fit_stage), main ="Fases de tratamento: 3 vs 7")
volcanoplot(fit_stage,coef="0_8",highlight=2, names=row.names(fit_stage), main ="Fases de tratamento: 0 vs 8")
volcanoplot(fit_stage,coef="1_8",highlight=2, names=row.names(fit_stage), main ="Fases de tratamento: 1 vs 8")
volcanoplot(fit_stage,coef="2_8",highlight=4, names=row.names(fit_stage), main ="Fases de tratamento: 2 vs 8")
volcanoplot(fit_stage,coef="3_8",highlight=3, names=row.names(fit_stage), main ="Fases de tratamento: 3 vs 8")
volcanoplot(fit_stage,coef="4_8",highlight=5, names=row.names(fit_stage), main ="Fases de tratamento: 4 vs 8")
volcanoplot(fit_stage,coef="5_8",highlight=2, names=row.names(fit_stage), main ="Fases de tratamento: 5 vs 8")
volcanoplot(fit_stage,coef="6_8",highlight=2, names=row.names(fit_stage), main ="Fases de tratamento: 6 vs 8")

##DE de acordo com o nº de tratamentos
table(metadados$CUMULATIVE_TREATMENT_TYPE_COUNT)

metadados$CUMULATIVE_TREATMENT_TYPE_COUNT <- gsub(0, "Zero", metadados$CUMULATIVE_TREATMENT_TYPE_COUNT)
metadados$CUMULATIVE_TREATMENT_TYPE_COUNT <- gsub(1, "Um", metadados$CUMULATIVE_TREATMENT_TYPE_COUNT)
metadados$CUMULATIVE_TREATMENT_TYPE_COUNT <- gsub(2, "Dois", metadados$CUMULATIVE_TREATMENT_TYPE_COUNT)
metadados$CUMULATIVE_TREATMENT_TYPE_COUNT <- gsub(3, "Três", metadados$CUMULATIVE_TREATMENT_TYPE_COUNT)
metadados$CUMULATIVE_TREATMENT_TYPE_COUNT <- gsub(4, "Quatro", metadados$CUMULATIVE_TREATMENT_TYPE_COUNT)
metadados$CUMULATIVE_TREATMENT_TYPE_COUNT <- gsub(5, "Cinco", metadados$CUMULATIVE_TREATMENT_TYPE_COUNT)


metadados$CUMULATIVE_TREATMENT_TYPE_COUNT = as.factor(metadados$CUMULATIVE_TREATMENT_TYPE_COUNT)
table(metadados$CUMULATIVE_TREATMENT_TYPE_COUNT)

design_ntreat = model.matrix(~0+metadados$CUMULATIVE_TREATMENT_TYPE_COUNT, data = y$samples)
colnames(design_ntreat) = levels(metadados$CUMULATIVE_TREATMENT_TYPE_COUNT)
head(design_ntreat)

contrast_ntreat = makeContrasts('0_1'= Zero - Um, '0_2'= Zero - Dois, '0_3'= Zero - Três, '0_4'= Zero - Quatro, '0_5'= Zero - Cinco,
                                '1_2'= Um - Dois, '1_3'= Um - Três, '1_4'= Um - Quatro, '1_5'= Um - Cinco,
                                '2_3'= Dois - Três, '2_4'= Dois - Quatro, '2_5'= Dois - Cinco,
                                '3_4'= Três - Quatro, '3_5'= Três - Cinco,
                                '4_5'= Quatro - Cinco,
                                levels=design_ntreat)

fit_ntreat = lmFit(y$counts, design_ntreat)
fit_ntreat = contrasts.fit(fit_ntreat, contrast_ntreat)
fit_ntreat = eBayes(fit_ntreat, trend=T)
topTable(fit = fit_ntreat, genelist = y$genes)
summa.fit_ntreat = decideTests(fit_ntreat)
summary(summa.fit_ntreat)

par(mfrow=c(2, 2))
volcanoplot(fit_ntreat,coef="0_5",highlight=2, names=row.names(fit_ntreat), main ="Tipos de tratamento: 0 vs 5")
volcanoplot(fit_ntreat,coef="1_5",highlight=4, names=row.names(fit_ntreat), main ="Tipos de tratamento: 1 vs 5")
volcanoplot(fit_ntreat,coef="2_5",highlight=4, names=row.names(fit_ntreat), main ="Tipos de tratamento: 2 vs 5")
volcanoplot(fit_ntreat,coef="3_5",highlight=4, names=row.names(fit_ntreat), main ="Tipos de tratamento: 3 vs 5")

############ Multi-Dimensional Scaling(MDS) ####################
dist_euclidiana <- dist(logCPM)
mds_result_isomds <- isoMDS(dist_euclidiana)

plot(mds_result_isomds$points, main = "Multi-Dimensional Scaling(MDS)")



################# Clustering #############
hierarquico - sexo}
# Preparação dos dados: 
DADOS = t(logCPM) #Transpor os dados
metadados$SEX = as.factor(metadados$SEX)
logCPM_metadados = cbind(DADOS, metadados)

#Clustering
logCPM_metadados.sc = scale(logCPM_metadados[,1:14405])
dist.logCPM_metadados = dist(logCPM_metadados.sc, method = "manhattan")
cl.logCPM_metadados = hclust(dist.logCPM_metadados, method = "average")
my.plot.hc = function(hclust, lab = 1:length(hclust$order), lab.col = rep(1, length(hclust$order)), hang = 0.1, ...)
  
my.plot.hc(cl.logCPM_metadados, lab.col = as.integer(logCPM_metadados$SEX)+1, cex = 0.6)

#preparação dos dados
logCPM_matrix <- as.matrix(logCPM)
class(logCPM_matrix) <- 'numeric'
max_logCPM_matrix <- apply(logCPM_matrix, 1, max)
min_logCPM_matrix <- apply(logCPM_matrix, 1, min)
vl <-max_logCPM_matrix/min_logCPM_matrix > 2
logCPM_matrix <-logCPM_matrix[vl, ]
logCPM_matrix <-na.exclude(logCPM_matrix)

##filtrar o top 50 de genes
expressao <- gsub("\\..*", "", rownames(logCPM_matrix))
rownames(logCPM_matrix) <- expressao
rownames(logCPM_matrix) <- make.names(rownames(logCPM_matrix), unique = TRUE)

test <- rowttests(logCPM_matrix)
ranking <- order(test$p.value)
p_value <- ranking[1:50]
genes_top50 = logCPM_matrix[p_value, ]

##matriz de distancias
distxy <- dist(genes_top50, method = 'euclidean')
hc = hclust(distxy, method = 'complete')
plot(hc)

hcd <- as.dendrogram(hc)
nodePar <- list(lab.cex = 0.6, pch = c(NA,20), cex = 0.8, col = "lightgreen")
plot(hcd, ylab = "Height", nodePar = nodePar, horiz = TRUE, edgePar = list(col = 2:3, lwd = 2:1))

#Heatmap
heatmap(genes_top50, labCol = sort(metadados$SEX))
heatmap(genes_top50, labCol = sort(metadados$OS_STATUS))

#Clustering
set.seed(123)
logCPM_matrix <- logCPM_matrix[1:length(metadados$OS_STATUS)]
centers <- cut(logCPM_matrix, breaks = 4)
centers_factor <- factor(centers)
live_factor <- factor(metadados$OS_STATUS)

levels(centers_factor) <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")

table <- table(centers_factor, live_factor)
print(table)


######### Redução da dimensionalidade ##################
# União dos dados com metadados relevantes para análises:
metadados_select = metadados[,c(8, 30)] # Seleciona o SAMPLE_SITE e o OS_STATUS.
merged_df <- merge(DADOS, metadados_select, by = "row.names")
row.names(merged_df) <- merged_df$Row.names
merged_df = merged_df[, -1]


# União com os outros metadados do sample_data_rna: "CUMULATIVE_TREATMENT_STAGE_COUNT" e "CUMULATIVE_TREATMENT_TYPE_COUNT": 
metadados_counts = sample_data_rna[ , c("CUMULATIVE_TREATMENT_STAGE_COUNT", "CUMULATIVE_TREATMENT_TYPE_COUNT")]
merged_df <- merge(merged_df, metadados_counts, by = "row.names")
row.names(merged_df) <- merged_df$Row.names
merged_df = merged_df[, -1]

## determinação das métricas
# Formulação do PCA dos dados:
dados_pca = merged_df[, -c(14406:14409)]
pca_result <- prcomp(dados_pca, scale = F)

summary(pca_result)
pca_result$rotation
pca_result$x

biplot(pca_result)

#Visualizar a distribuição dos pontos no PC1 e o PC2:
pca_data <- as.data.frame(pca_result$x[, 1:2])
plot(pca_data$PC1, pca_data$PC2,
     main = "Projeção nos Dois Primeiros Componentes Principais",
     xlab = "PC1",
     ylab = "PC2",
     pch = 19,  
     col = "blue") 

# determinação de qunatos PCAs explicam 90% dos dados:
prop_var <- summary(pca_result)$importance["Proportion of Variance",]
cumulative_prop_var <- cumsum(prop_var)
n_components_90 <- which(cumulative_prop_var >= 0.9)[1]
n_components_90 

# Análise da qualidade das variáveis individuais nos dois primeiros componentes principais
fviz_famd_ind(pca_result, geom = c("point"), col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), palette = "rainbow", addEllipses = FALSE, ellipse.type = "confidence",ggtheme = theme_minimal(), repel = TRUE, labels = F)  

# Análise dos dados para diferentes metadados:
# Sample_site:
autoplot(pca_result, data = as.data.frame(merged_df), colour = 'SAMPLE_SITE', size = 2, main = "SAMPLE_SITE")

# Plotegem de um gréfico com 3 dimensões
pca_data <- pca_result$x[, 1:3]
scatterplot3d(pca_data, color = as.integer(merged_df$SAMPLE_SITE), pch = 16,
              xlab = "PC1", ylab = "PC2", zlab = "PC3", main = "Representação de 3 dimensões - SAMPLE_SITE")

# "CUMULATIVE_TREATMENT_STAGE_COUNT"
autoplot(pca_result, data = as.data.frame(merged_df), colour = "CUMULATIVE_TREATMENT_STAGE_COUNT", size = 2, main = "CUMULATIVE_TREATMENT_STAGE_COUNT")

# "CUMULATIVE_TREATMENT_TYPE_COUNT"
autoplot(pca_result, data = as.data.frame(merged_df), colour = "CUMULATIVE_TREATMENT_TYPE_COUNT", size = 2, main = "CUMULATIVE_TREATMENT_TYPE_COUNT")

# "OS_STATUS"
legenda_personalizada <- c("No inf" = "gray", "0:LIVING" = "red", "1:DECEASED" = "green")
autoplot(pca_result, data = as.data.frame(merged_df), colour = "OS_STATUS", size = 2, main = "OS_STATUS") +
  scale_color_manual(values = legenda_personalizada) # Cinza são os NAs. 

scatterplot3d(pca_data, color = as.integer(merged_df$OS_STATUS), pch = 16,
              xlab = "PC1", ylab = "PC2", zlab = "PC3", main = "Representação de 3 dimensões - OS_STATUS")

#################### Machine Learning #######################
# Definir a seed para caso se queira reproduzir estes resultados
set.seed(9999)
cvcontrol = trainControl( method = 'repeatedcv', number = 10, repeats = 5)
merged_df$AGE_AT_DIAGNOSIS <- metadados[rownames(merged_df), "AGE_AT_DIAGNOSIS"]
merged_df$CANCER_TYPE <- metadados[rownames(merged_df), "CANCER_TYPE"]
merged_df$SEX <- metadados[rownames(merged_df), "SEX"]
merged_df$ETHNICITY <- metadados[rownames(merged_df), "ETHNICITY"]
ml_data = merged_df
ml_data$sobrevida_paciente = as.factor(ml_data$OS_STATUS)

set.seed(9999)
indices <- createDataPartition(ml_data$sobrevida_paciente, p = 0.7, list = FALSE)
dados_treino <- ml_data[indices, ]
dados_teste <- ml_data[-indices, ]
dados_teste$sobrevida_paciente[dados_teste$sobrevida_paciente == ""] <- NA
dados_teste_NAs <- na.omit(dados_teste)
dados_treino$sobrevida_paciente[dados_treino$sobrevida_paciente == ""] <- NA
dados_treino_NAs <- na.omit(dados_treino)
ml_data$sobrevida_paciente[dados_teste$sobrevida_paciente == ""] <- NA
ml_data_NAs <- na.omit(ml_data)

#Arvores de decisão
# Definir o modelo
modelo_dt = ctree(sobrevida_paciente ~ ., data = dados_treino_NAs)
# Visualizar o modelo
print(modelo_dt)
plot(modelo_dt)
testDT = predict(modelo_dt, dados_teste_NAs)

# Comparação dos resultados utilizando uma tabela
table(testDT, dados_teste_NAs$sobrevida_paciente)

# Calculo da precisão
precisao_DT = mean(testDT == dados_teste_NAs$sobrevida_paciente)
precisao_DT

### Naive-Bayes
modelo_nb = naiveBayes(sobrevida_paciente ~ ., data = dados_treino_NAs)
testeNB = predict(modelo_nb, dados_teste_NAs)
table(testeNB, dados_teste_NAs$sobrevida_paciente)
precisao_NB = mean(testeNB == dados_teste_NAs$sobrevida_paciente)
precisao_NB