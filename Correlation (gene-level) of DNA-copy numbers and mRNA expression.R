library(readxl)
xl = read_excel("C:/Users/YJ/R.practice/Data/TCGA-CDR-SupplementalTableS1.xlsx")
xl = as.data.frame(xl)
xl = xl[, c(2, 3)]
tumortype_tcga = xl[, -1]
names(tumortype_tcga) = xl[, 1]
head(tumortype_tcga)

#GBM type tumor 선정

a = readLines("C:/Users/YJ/R.practice/Data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv", n = 1)
head(a)
header = unlist(strsplit(gsub('\"', "", a), "\t"))
head(header)
header = header[-1]

gbm_id = names(tumortype_tcga)[tumortype_tcga == "GBM"]
gbm_id

gbm_exp = read.table("C:/Users/YJ/R.practice/Data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv", sep = "\t", header = T, colClasses = c("character", ifelse(substr(1, 15, header) %in% paste0(gbm_id, "-01"), "numeric", "NULL")))
gbm_exp[1:5, 1:5]
rownames(gbm_exp) = gbm_exp[,1]
gbm_exp = gbm_exp[, -1]

genename = gsub("\\|.*","", rownames(gbm_exp))
genename

gbm_exp = gbm_exp[!duplicated(genename),]
rownames(gbm_exp) = genename[!duplicated(genename)]

#6분

tcga.seg = read.table("C:/Users/YJ/R.practice/Data/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg", header = T, sep = "\t")
head(tcga.seg)

tcga.seg.sample = unique(tcga.seg$Sample)
tcga.seg.sample
#unique함수는 중복값 제거

tumortype_seg = tumortype_tcga[substr(tcga.seg.sample,1,12)]
tumornormal_seg = ifelse(substr(tcga.seg.sample, 14, 15) == "01", "tumor", "others")
tumornormal_seg
sampleAnnotation = cbind(tcga.seg.sample, tumortype_seg, tumornormal_seg)
head(sampleAnnotation)
colnames(sampleAnnotation) = c("sample", "tumor type", "tumor normal")
write.table(sampleAnnotation, "C:/Users/YJ/R.practice/Data/tcga_seg_sampleAnnotation.txt", sep = "\t", quote = F, row.names = F)

# IGV browser 확인 결과:
# 7번 chromosome(EGFR gene)에서 copy number 증폭
# 9번 chromosome(CDKN2A gene)에서 copy number 감소
# 10번 chromosome에서 copy number 감소

gbm_tumor_seg = tcga.seg[substr(tcga.seg.sample, 1, 12) %in% gbm_id & substr(tcga.seg.sample, 14, 15) == "01",]
head(gbm_tumor_seg)
write.table(gbm_tumor_seg, "C:/Users/YJ/R.practice/Data/GBMseg_forGISTIC.seg", sep = "\t", quote = F, row.names = F)
# GenePattern에서 GIstic 알고리즘으로 driver gene 찾기

#### purity, ploid 분석 ####

tcga_absolute = read.table("C:/Users/YJ/R.practice/Data/TCGA_mastercalls.abs_tables_JSedit.fixed.txt", header = T, sep = "\t")
tcga_absolute = tcga_absolute[substr(tcga_absolute$array,14,15) == "01", ]
dim(tcga_absolute)
tumor_type = c("LUAD", "LUSC", "HNSC", "KIRC", "BRCA", "BLCA", "UCEC", "GBM", "OV")
samples = names(tumortype_tcga)[tumortype_tcga == tumor_type]
length(samples)
length(tumortype_tcga[match(samples, names(tumortype_tcga))])
names(samples) = tumortype_tcga[match(samples, names(tumortype_tcga))]
pp_data = tcga_absolute[, c(1,4,5,6)]
pp_data
table(pp_data$Genome.doublings)
pp_data = pp_data[pp_data$array %in% paste0(samples, "-01"), ]
rownames(pp_data) = pp_data[,1]
pp_data = pp_data[, -1]
hist(pp_data$purity, breaks = 30, ylim = c(0, 60), xlab = "PURITY", ylab = "Sample", main = "All_Lineages")

par(mfrow = c(1, 9))
for(tumor in tumor_type) {
  sam = names(tumortype_tcga)[tumortype_tcga == tumor]
  data = pp_data[rownames(pp_data) %in% paste0(sam, "-01"),]
  boxplot(data$purity, main = tumor, xlab = "PURITY", ylab = "Sample")
}

par(mfrow = c(1,9))
for(tumor in tumor_type) {
  sam = names(tumortype_tcga)[tumortype_tcga == tumor]
  data = pp_data[rownames(pp_data) %in% paste0(sam, "-01"),]
  plot(data$ploidy, pch = 16, main = tumor, ylim = c(1,5), ylab = "PLOIDY", xlab = "", col =
         ifelse(data$Genome.doublings < 1, "purple",
                ifelse(data$Genome.doublings < 2, "darkgreen", "red")))
}