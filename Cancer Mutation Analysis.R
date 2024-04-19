tcga.mut = read.csv("C:/Users/YJ/R.practice/Data/mc3.v0.2.8.PUBLIC.maf", header = T, sep = "\t", colClasses = c("character", "NULL", "NULL", "NULL", "character", "integer", "integer", "NULL", "character", "character", "character", "NULL", "character", "NULL", "NULL", "character", "character", rep("NULL", 17), "character", "NULL", "character", "NULL", "NULL", rep("integer", 6), rep("NULL", 69)))
dim(tcga.mut)
head(tcga.mut)
# Hugo_Symbol: cancer mutation이 발생한 유전자
# Tumor_Sample_Barcode: 어느 환자에서 나타났는지
# HGVSp_Short: 해당 유전자 속 어느 아미노산에 주요 mut이 발생했는지


#### TMB: 한 환자가 갖는 mutation frequency ####

length(unique(tcga.mut$Tumor_Sample_Barcode))
# 총 10295명의 환자 데이터
head(table(tcga.mut$Tumor_Sample_Barcode))
# 각 환자가 갖는 mutation의 총 개수 확인

boxplot(log(sort(table(tcga.mut$Tumor_Sample_Barcode)), 2))

library(readxl)
a = read_excel("C:/Users/YJ/R.practice/Data/TCGA-CDR-SupplementalTableS1.xlsx")
a = as.data.frame(a)
a = a[,c(2,3)]
tumortype_tcga = a[, -1]
tumortype_tcga
names(tumortype_tcga) = a$bcr_patient_barcode
head(tumortype_tcga)

tumortype_tmb = table(tumortype_tcga)
names(tumortype_tmb)[tumortype_tmb > 500]
#환자가 500명 이상인 tumor type 10개 확보

par(mfrow = c(1, 10))

for( tumortype in names(tumortype_tmb)[tumortype_tmb > 500]) {
  samples = names(tumortype_tcga)[tumortype_tcga == tumortype]
  samples_im_mut = tcga.mut$Tumor_Sample_Barcode[substr(tcga.mut$Tumor_Sample_Barcode,1,15) %in% paste0(samples, "-01")]
  plot(log(sort(table(samples_im_mut)),2), type = "p", ylab = "", xlab = "", main = tumortype, xaxt="n")
}

#### 유전자 기준 mutation abundance ####

table(tcga.mut$Hugo_Symbol)
tail(sort(table(tcga.mut$Hugo_Symbol)))
# mutation이 많이 발생한 유전자 6개
#TTN은 가장 긴 유전자, 따라서 mutation frequency가 높음
# cancer driver gene 찾는 것이 포인트 (Mutsig 알고리즘 사용)
# GenePattern -> Use GenePattern -> AWS -> module -> mutsigCV 검색

if(!requireNamespace("BiocManager", quietly = T)) {
  install.packages("BiocManager")
}
library(BiocManager)
BiocManager::install("GenVisR")
library(GenVisR)

head(tumortype_tcga)
LUSC_sample = names(tumortype_tcga)[tumortype_tcga == "LUSC"]
length(LUSC_sample)
LUSC_sample
head(tcga.mut)
LUSC_mut = tcga.mut[substr(tcga.mut$Tumor_Sample_Barcode,1,15) %in% paste0(LUSC_sample, "-01"),]
head(LUSC_mut)
colnames(LUSC_mut)[12]
colnames(LUSC_mut)[12] = "amino.acid.change"
LUSC_mut_filtered = LUSC_mut[LUSC_mut$Hugo_Symbol %in% c("KEAP1", "PTEN", "TP53", "CDKN2A", "MLL2", "NFE2L2", "RB1", "FBXW7", "NOTCH1", "PIK3CA"),]
LUSC_mut_filtered
waterfall(LUSC_mut_filtered, fileType = "MAF")

# amino acid 수준의 변화로 기능적 변화 유추

table(tcga.mut$Variant_Classification)
par(mar = c(12, 6, 2, 2))
plot(table(tcga.mut$Variant_Classification), type = 'h', las = 2)
# Missense_MUT이 가장 많고 다음은 Silent_MUT
 
#### LUAD(Lung Adenocarcionmas) and EGFR mutations ####

LUAD_sample = names(tumortype_tcga)[tumortype_tcga == "LUAD"]
LUAD_mut = tcga.mut[substr(tcga.mut$Tumor_Sample_Barcode, 1, 15) %in% paste0(LUAD_sample, "-01"),]
LUAD_mut_EGFR = LUAD_mut[LUAD_mut$Hugo_Symbol == "EGFR", ]
table(LUAD_mut_EGFR$HGVSp_Short)
# ex) p.G901V -> protein 수준에서 901번째 a.a G가 V로 바뀐 환자 수

write.table(LUAD_mut_EGFR[, c(2,3,4,7,8)], "C:/Users/YJ/R.practice/Data/LUAD_EGFR_mutations_forANNOVAR.txt", sep = "\t", quote = F, row.names = F, col.names = F)
LUAD_mut_EGFR[, c(2,3,4,7,8)]
# http://wannovar.wglab.org/ 에서 삭제된 Annot 확인 가능
# 즉, 아미노산 변화 확인 가능
# Sample Identifier에 "TCGA-LUAD-EGFR" 입력

#### mutation 빈도, mutation hotspot 그리기 ####

pik3ca_missense = tcga.mut[tcga.mut$Hugo_Symbol == "PIK3CA" & tcga.mut$Variant_Classification == "Missense_Mutation", ]
pik3ca_missense$HGVSp_Short
a = regexpr("[0-9]+", pik3ca_missense$HGVSp_Short, perl = T)
# 정규 표현식은 [0-9]+로, 이는 문자열에서 연속된 숫자를 찾는 것을 의미
# perl = TRUE 옵션은 Perl 스타일의 정규 표현식을 사용하도록 지정
# 결과는 패턴이 발견된 위치의 첫 번째 인덱스를 반환, 만약 해당 패턴이 문자열에 없다면 -1을 반환
b = regmatches(pik3ca_missense$HGVSp_Short, a)
# regmatches 함수는 regexpr 또는 gregexpr 함수를 사용하여 찾은 패턴에 대해 일치하는 부분 문자열을 추출하는 데 사용
hist(as.numeric(b), breaks = 100)

bp = boxplot(as.numeric(b))
bp$stats
as.numeric(b)[1068]
# PIK3ca gene에서 542번 아미노산에 Missense mutation이 가장 많이 발생

tiff(file="C:/Users/YJ/R.practice/MUT_Hotspot_in_Gene.tiff", width=50, height=40, units="cm", res=300)
par(mfrow = c(6,1), mar = c(2,2,2,1))
genes = c("PIK3CA", "CTNNB1", "BRAF", "KRAS", "EGFR", "ERBB2")
for(gene in genes) {
  gene_missense = tcga.mut[tcga.mut$Hugo_Symbol == gene & tcga.mut$Variant_Classification == "Missense_Mutation", ]
  a = regexpr("[0-9]+", gene_missense$HGVSp_Short, perl = T)
  b = regmatches(gene_missense$HGVSp_Short, a)
  barplot(table(as.numeric(b)), col = "red", main = gene)
}
dev.off()