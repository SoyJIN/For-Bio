a = readLines("C:/Users/YJ/R.practice/Data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv", n = 1)
head(a)
header = unlist(strsplit(gsub('\"', "", a), "\t"))
header = header[-1]
head(header)

length(header)
length(substr(header, 1, 12))
# substr(문자, 시작 번호, 끝 번호)
# 주어진 문자에서 시작 번호와 끝 번호에 해당하는 문자 추출

install.packages("readxl")
library(readxl)
tumortype_tcga = read_excel("C:/Users/YJ/R.practice/Data/TCGA-CDR-SupplementalTableS1.xlsx")
head(tumortype_tcga)
tumortype_tcga = as.data.frame(tumortype_tcga)
tumortype_tcga = tumortype_tcga[, c(2,3)]
tumortype_tcga_ = tumortype_tcga[,-1]
names(tumortype_tcga_) = tumortype_tcga$bcr_patient_barcode
tumortype_tcga_

brca_id = names(tumortype_tcga_)[tumortype_tcga_ == "BRCA"]
length(brca_id)

brca_exp = read.table("C:/Users/YJ/R.practice/Data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv", header = T, sep = "\t", colClasses = c("character", ifelse(substr(header, 1, 15) %in%  paste0(brca_id, "-01"), "numeric", "NULL")))
brca_exp[1:5, 1:5]
dim(brca_exp)
head(brca_exp)

rownames(brca_exp) = brca_exp[,1]
brca_exp = brca_exp[,-1]

genename = unlist(lapply(strsplit(rownames(brca_exp),"\\|"), '[[', 1))
#[[, 1은 각 나뉜 요소들 중 첫 번째 요소만을 추출하는 것
genename

brca_exp = brca_exp[!duplicated(genename), ]
dim(brca_exp)
#중복되지 않은 값은 TRUE로 표시되고, 중복된 값은 FALSE로 표시

rownames(brca_exp) = genename[!duplicated(genename)]
brca_exp[1:5, 1:5]

b = readLines("C:/Users/YJ/R.practice/Data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv", n = 1)
head(b)
header = unlist(strsplit(b, ","))
header
header = header[-(1:2)]

length(header)
length(substr(header, 1, 12))

brca_miR_exp = read.table("C:/Users/YJ/R.practice/Data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv", header = T, sep = ",", colClasses = c("character", ifelse(substr(header, 1, 15) %in% paste0(brca_id, "-01"), "numeric", "NULL")))
brca_miR_exp[1:5, 1:5]
rownames(brca_miR_exp) = brca_miR_exp[,1]
brca_miR_exp = brca_miR_exp[, -1]

save("tumortype_tcga_", "brca_id", "brca_exp", "brca_miR_exp", file = "C:/Users/YJ/R.practice/Data/Multiomics_module_1.R")
load(file = "C:/Users/YJ/R.practice/Data/Multiomics_module_1.R")

brca_exp[1:20,1:2]
brca_miR_exp[1:20, 1:2]

colnames(brca_exp) = substr(colnames(brca_exp), 1, 12)
colnames(brca_miR_exp) = substr(colnames(brca_miR_exp), 1, 12)

common_sample = intersect(colnames(brca_exp), colnames(brca_miR_exp))
length(common_sample)
brca_exp_matching = brca_exp[, common_sample]
brca_miR_exp_matching = brca_miR_exp[, common_sample]
dim(brca_exp_matching)
dim(brca_miR_exp_matching)

mRNA_miR_cor = cor(t(brca_exp_matching), t(brca_miR_exp_matching), method = "pearson")
mRNA_miR_cor[1:5, 1:5]
hist(mRNA_miR_cor, breaks = 100)

mRNA_miR_cor_random = cor(t(brca_exp_matching), t(brca_miR_exp_matching[, sample(1:length(common_sample))]), method = "pearson")

random = hist(mRNA_miR_cor_random, breaks = 100, plot = F)
#plot = F는 플롯(그래프)을 그리지 않도록 하는 옵션
random
lines(random$mids, random$counts, col = "orange", lwd = 5)
mir_target = read_excel("C:/Users/YJ/R.practice/Data/miRTarBase_MTI.xlsx")
mir_target[1:5, 1:5]
class(mir_target)
mir_target = as.data.frame(mir_target)
mir_target = mir_target[,c(2,3,4)]
head(mir_target)

save("tumortype_tcga_", "brca_id", "brca_exp", "brca_miR_exp", "common_sample", "brca_exp_matching", "brca_miR_exp_matching", "mRNA_miR_cor", "mir_target", file = "C:/Users/YJ/R.practice/Data/Multiomics_module_1.R")
load(file = "C:/Users/YJ/R.practice/Data/Multiomics_module_1.R")

cor_validatedPairs = c()
for(i in 1:nrow(mir_target)) {
  if(mir_target$miRNA[i] %in% colnames(mRNA_miR_cor) & mir_target$`Target Gene`[i] %in% rownames(mRNA_miR_cor)) {
    cor_validatedPairs = c(cor_validatedPairs, mRNA_miR_cor[match(mir_target$`Target Gene`[i], rownames(mRNA_miR_cor)), match(mir_target$miRNA[i], colnames(mRNA_miR_cor))])
  }
}
head(cor_validatedPairs)

validated = hist(cor_validatedPairs, breaks = 100, plot = F)
validated
lines(validated$mids,validated$counts*120, col = "red", lwd = 5)
