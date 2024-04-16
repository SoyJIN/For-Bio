if(!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("pheatmap", "RColorBrewer", "edgeR", "DESeq2", "msigdbr", "fgsea", "xlsx"))
library(pheatmap)
library(RColorBrewer)
library(edgeR)
library(DESeq2)
library(msigdbr)
library(fgsea)
library(xlsx)

adr = "C:/Users/YJ/R.practice/Data"
list.files(adr)

s_info = read.csv(file = sprintf("%s/sample_info.csv", adr), header = T, check.names = F, stringsAsFactors = F)
s_info[1:5, 1:5]

file_name = sprintf("%s/CCLE_expression_proteincoding_genes_expected_count.csv", adr)
rcm = read.csv(file = file_name, header = T, check.names = F, row.names = 1)
#read count matrix
rcm[1:5, 1:5]

tpm = read.csv(file = sprintf("%s/CCLE_expression.csv", adr), header = T, row.names = 1, check.names = F)
#유전자 발현 수준 matrix
tpm[1:5, 1:5]

dim(rcm)
dim(tpm)
rcm = t(rcm)
rcm[1:5, 1:5]
tpm = t(tpm)
tpm[1:5, 1:5]

rownames(rcm) = gsub(" \\(.*", "", rownames(rcm))
rownames(tpm) = gsub(" \\(.*", "", rownames(tpm))

rcm = round(rcm, digits = 0)

file_name2 = sprintf("%s/Drug_sensitivity_AUC_(CTD^2).csv", adr)
ds = read.csv(file = file_name2, header = T, check.names = F, row.names = 1, stringsAsFactors = F)
ds[1:5, 1:5]
dim(ds)
colnames(ds) = gsub(" \\(CTRP.*", "", colnames(ds))

identical(colnames(rcm), colnames(tpm))
icell = intersect(colnames(rcm), rownames(ds))
length(icell)

rcm = rcm[ , match(icell, colnames(rcm))]
tpm = tpm[ , match(icell, colnames(tpm))]
ds = ds[match(icell, rownames(ds)), ]
s_info = s_info[match(icell, s_info$DepMap_ID), ]

save(s_info, rcm, tpm, ds, file = sprintf("%s/Proceed_Data", adr))
load(file = sprintf("%s/Proceed_Data", adr))

## Read Count
rcm[1:5, 1:5]
ri = sample(1:ncol(rcm), 9)
ri #암세포 9개 임의 추출

par(mfrow = c(3,3))
for (i in ri){
  hist(log2(rcm[,i]+1), xlab = "Read Count", main = colnames(rcm)[i], breaks = 30)
}

rcm_sd = apply(log2(rcm+1), 1, sd)
hist(rcm_sd, col = "skyblue", xlab = "SD of Read Count")

## 유전자 발현 수준
tpm[1:5, 1:5]
par(mfrow = c(3,3))
for (i in ri){
  hist(tpm[,i], xlab = "log2(TPM+1)", main = colnames(tpm)[i], breaks = 30,  xlim = c(0, 10))
}

tpm_sd = apply(tpm, 1, sd)
hist(tpm_sd, col = "pink", xlab = "SD of TPM")

## 약물 반응성
ds[1:5, 1:5]
dim(ds)

rd = sample(1:ncol(ds), 9)
rd #약물 9개 임의 추출
par(mfrow = c(3,3))
for (i in rd){
  hist(ds[,i], breaks = 30,  main = colnames(ds)[i], xlab = "AUC")
}

auc_sd = apply(ds, 2, sd, na.rm = T)
hist(auc_sd, col = "purple", xlab = "SD of AUC", breaks = 20)

which(auc_sd == max(auc_sd))
dsv = ds[,"OLIGOMYCIN A"]
dsv
names(dsv) = rownames(ds)
dsv = dsv[!is.na(dsv)]

hist(dsv, breaks = 30, col = "orange", xlab = "AUC", main = "Oligomycin A", prob = T)
lines(density(dsv))

q2 = quantile(dsv, prob = c(0.01, 0.99))
abline(v = q2, col = "black", lty = 5, lwd = 2)
text(x=q2[2], y=0.14, "<R", xpd = T, col = "red")
text(x=q2[1], y=0.14, "S>", xpd = T, col = "blue")

data.frame(response = "S", CCL = names(dsv)[q2[1]>=dsv])
data.frame(response = "R", CCL = names(dsv)[dsv>=q2[2]])

cinfo = rbind(data.frame(response = "S", CCL = names(dsv)[q2[1]>=dsv]), data.frame(response = "R", CCL = names(dsv)[dsv>=q2[2]]))
cinfo
colnames(s_info)
cinfo = cbind(cinfo,s_info[match(cinfo$CCL,s_info$DepMap_ID), c( "primary_disease", "lineage", "cell_line_name")])
cinfo

tpm[1:5, 1:5]
stpm = tpm[,match(cinfo$CCL, colnames(tpm))]
dim(stpm)

rcm[1:5, 1:5]
srcm = rcm[,match(cinfo$CCL,colnames(rcm))]
dim(srcm)

save(cinfo, stpm, srcm, file = sprintf("%s/Oligomycin_SR_transcriptome.R",adr))
load(file = sprintf("%s/Oligomycin_SR_transcriptome.R", adr))

stpm[1:5, 1:5]
dim(stpm)
sds = apply(stpm, 1, sd)
hist(sds, breaks = 30, col = "skyblue")

xx = stpm[which(sds>=0.5), ]
xx[1:5, 1:5]
dim(xx)

color = c("darkblue", "darkred")
names(color) = c("S", "R")
cls = color[cinfo$response]
pca = prcomp(t(xx), scale = T)
summary(pca)
pca
plot(pca$x[,1], pca$x[,2], col = cls, pch = 19, cex = 1, xlab = "PCA1", ylab = "PCA2", main = paste("PCA - ", length(rownames(xx)), "Genes"), cex.main = 1.2)
text(x = pca$x[,1], y = pca$x[,2]+3, cinfo$primary_disease, col = cls, cex = 1)
legend("topright", names(color), col = color, pch = 19)

cm = cor(xx, method = "pearson")
cm
cinfo
annot = cinfo[,c(1,3,4,5)]
annot
rownames(annot) = cinfo$CCL

hcls = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100)
pheatmap(cm, col = hcls, annotation=annot[,c(1,3)],labels_col=annot$cell_line_name, labels_row=annot$cell_line_name,  clustering_distance_rows = "correlation", clustering_distance_cols = "correlation")

grp = factor(cinfo$response, levels = c("S", "R"))
keep = filterByExpr(srcm, group = grp)
keep
srcm
srcm = srcm[which(keep),]

coldata = data.frame(Group = grp)
rownames(coldata) = colnames(srcm)
coldata

dds = DESeqDataSetFromMatrix(srcm, coldata, design = ~Group)
dds = DESeq(dds)
degm = results(dds, pAdjustMethod = "fdr", independentFiltering = F)
degm = as.data.frame(degm)
degm

identical(rownames(degm), rownames(srcm))
srcm

degm = data.frame(srcm, degm, check.names = F)
degm
degm = degm[order(degm$padj),]

par(mfrow = c(2,3))
for(g in rownames(degm)[1:6]){
  plotCounts(dds, gene = g, intgroup = "Group")
}

save(degm, file=sprintf("%s/DESeq2_DEG.R", adr))
write.table(degm, file = sprintf("%s/DESeq2_DEG.txt", adr), col.names=T, sep="\t", quote=F)

degm$logFDR = -log10(degm$padj)
plot(degm$log2FoldChange,degm$logFDR, xlab = "log2FC", ylab = "log10FDR", pch = 21, col = "black", bg = "darkgray", xlim = c(-10, 10))
abline(v = c(-log2(1.5), log2(1.5)), lty = 3, col ="red")
abline(h = 2, lty = 3, col = "red")

uidx = degm$log2FoldChange >= log2(1.5) & degm$logFDR >= 2
didx = degm$log2FoldChange <= -log2(1.5) & degm$logFDR >= 2

points(x=degm$log2FoldChange[uidx], y=degm$logFDR[uidx], col = "black", bg = "red", pch =21)
points(x=degm$log2FoldChange[didx], y=degm$logFDR[didx], col = "black", bg = "blue", pch = 21)
legend("topleft", c("UP-DEGs", "Down-DEGs"), pt.bg = c("red", "blue"), pch = 21)

uidx = degm$log2FoldChange >= log2(1.5) & degm$logFDR >= 3
didx = degm$log2FoldChange <= -log2(1.5) & degm$logFDR >= 3
text(x = degm$log2FoldChange[uidx], y = degm$logFDR[uidx] + 0.3, rownames(degm)[uidx], col = "red", cex = 0.7)
text(x = degm$log2FoldChange[didx], y = degm$logFDR[didx] + 0.3, rownames(degm)[didx], col = "blue", cex = 0.7)

degs = rownames(degm)[uidx|didx]
xx = stpm[rownames(stpm) %in% degs, ]
annot = cinfo[, c(1,3,4,5)]
rownames(annot) = cinfo$CCL
hcls = colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(100)
pheatmap(xx, col=hcls, scale="row", annotation = annot[,c(1,3)], clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", labels_col=annot$cell_line_name, fontsize=8)
