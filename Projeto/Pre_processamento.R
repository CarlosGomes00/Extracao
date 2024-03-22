# {r Pré-processamento dos dados}
# Carregamento dos dados através do url: https://raw.githubusercontent.com/lais-carvalho/Trabalho-pratico/main/Dados/

######### RNA_CPM_ZScore #########################################################################################
# Primeiro associamos o url a uma variável
RNA_cpm_zscore_url = "https://raw.githubusercontent.com/lais-carvalho/Trabalho-pratico/main/Dados/data_mrna_seq_cpm_zscores_ref_all_samples.txt"

# Depois fazemos o download do ficheiro no url e guardamos uma cópia local 
download.file(RNA_cpm_zscore_url, destfile="RNA_seq_cpm_zscore.txt")  

# Leitura e visualização do ficheiro 
rna_cpm_zscore = read.table("RNA_seq_cpm_zscore.txt", header = T, sep = '\t')
View(rna_cpm_zscore)

# Verificar a estrutura dos dados
str(rna_cpm_zscore)
class(rna_cpm_zscore)

# Verificar se existem missing values
any(is.na(rna_cpm_zscore))

# Dimensão da tabela
dim(rna_cpm_zscore)

# Verificar se as colunas são as corretas
colnames(rna_cpm_zscore)

# Ver os primeiros valores para verificar se estão corretos
rna_cpm_zscore[1:5,1:5]

# Suspeitamos que os valores da coluna 'Entrez_gene_id' são os únicos Missing Values
# Então, a coluna será eliminada

rna_cpm_zscore$Entrez_Gene_Id <- NULL

# Verificar outra vez se ainda há missing values

any(is.na(rna_cpm_zscore))

# Já não existem NA no nosso dataset



#### Repetir o procedimento para os outros datasets
######### RNA_RPKM_ZScore ########################################################################################

RNA_rpkm_zscore_url = "https://raw.githubusercontent.com/lais-carvalho/Trabalho-pratico/main/Dados/data_mrna_seq_rpkm_zscores_ref_all_samples.txt"

download.file(RNA_rpkm_zscore_url, destfile="RNA_seq_rpkm_zscore.txt")  

rna_rpkm_zscore = read.table("RNA_seq_rpkm_zscore.txt", header = T, sep = '\t')
View(rna_rpkm_zscore)

str(rna_rpkm_zscore)
class(rna_cpm_zscore)

any(is.na(rna_rpkm_zscore))

dim(rna_rpkm_zscore)

colnames(rna_rpkm_zscore)

rna_rpkm_zscore[1:5,1:5]

rna_rpkm_zscore$Entrez_Gene_Id <- NULL

rna_rpkm_zscore[1:5,1:5]

any(is.na(rna_cpm_zscore))

# Também já não existem NA após o tratamento!!!


######### data_clinical_patient ########################################################################################

data_patient_url = "https://raw.githubusercontent.com/lais-carvalho/Trabalho-pratico/main/Dados/data_clinical_patient.txt"
download.file(data_patient_url, destfile="data_patient.txt")  

data_patient = read.table("data_patient.txt", header = T, sep = '\t')
View(data_patient)

str(data_patient)
class(data_patient)

any(is.na(data_patient))

# Existem muitos NAs no meio dos resultados, e pontos sem valores
# Ver como tratar esta tabela!!!!


######## data_clinical_sample #####################################################################################

data_sample_url = "https://raw.githubusercontent.com/lais-carvalho/Trabalho-pratico/main/Dados/data_clinical_sample.txt"
download.file(data_sample_url, destfile ="data_sample.txt")

data_sample = read.table("data_sample.txt", header = T, sep = '\t')
View(data_sample)


# Não estou a conseguir dar load desta última tabela!


#preparação
#pré-processamento



