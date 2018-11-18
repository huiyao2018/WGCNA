####WGCNA test
#第一步：对表达矩阵进行处理，由于这个是测试数据，我就不进行筛选差异基因了。
#虽然我觉得不应该筛选差异基因再去做WGCNA，但是有些论文以及一些期刊是这么做的。。。
#我觉得只要进行变异系数过滤，低表达量过滤等等就够了。。。
#下面我就同博主视频里面的代码进行过滤，实在是原始数据基因数有点多了。。然后取3000个基因来分析
## 首先在GSE48213_RAW目录里面生成tmp.txt文件
#awk '{print FILENAME"\t"$0}' * | grep -v EnsEMBL_Gene_ID >tmp.txt
## 然后把tmp.txt导入R语言里面用reshape2处理即可！
setwd('F:/WGCNA')
# 第一步：对表达矩阵进行处理，由于这个是测试数据，我就不进行筛选差异基因了。
# 虽然我觉得不应该筛选差异基因再去做WGCNA，但是有些论文以及一些期刊是这么做的。。。
# 我觉得只要进行变异系数过滤，低表达量过滤等等就够了。。。下面我就同博主视频里面的代码
# 进行过滤，实在是原始数据基因数有点多了。。然后取3000个基因来分析

#安装管理生信分析包的BiocManager包
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("ape")


#使用*****NCBI*****上面下载的GSE_RAW数据（备注：sangerbox中的TGCA简易下载工具v14没搜出来数据所以不好下载)
### 首先在GSE48213_RAW(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48213)
### 目录里面生成tmp.txt文件,使用下面的linux命令
### awk '{print FILENAME"\t"$0}' * | grep -v EnsEMBL_Gene_ID > tmp.txt
data <- read.table(file = "tmp.txt", sep = "\t", stringsAsFactors = F)
library(reshape2)
RNAseq <- dcast(data, formula = V2~V1)
rownames(RNAseq) <- RNAseq[,1]
RNAseq <- RNAseq[,-1]
names(RNAseq) <- unlist(lapply(names(RNAseq), function(x){
  tmp <- strsplit(x, '_')[[1]][2]
  tmp <- strsplit(tmp, '\\.')[[1]][1]
})
)

####先使用2018-10-30日https://bioconductor.org/packages/devel/bioc/vignettes/CVE/inst/doc/WGCNA_from_TCGA_RNAseq.html
####网站更新的分析方法进行测试
load(file = "SKCM_expr_df.Rda") #导入来自tcga数据库的数据
expr_df[-which(substr(expr_df$gene_id,1,2) != "EN"),] #删除后面几行不正确的数据
RNAseq <- expr_df
#Further preprocessing included the removal of control samples (for more information
#see the TCGA Wiki) and expression estimates with counts in less than 20% of cases.
RNAseq = RNAseq[apply(RNAseq,1,function(x) sum(x==0)) < ncol(RNAseq)*0.8,]
# RNAseq data was normalised with the  voom methodology33 Charity W Law et al. 
# voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. 
# In: Genome biology 15.2 (Jan. 2014), R29–R29.. The voom method estimates the 
# mean-variance of the log-counts and generates a precision weight for each observation. 
# This way, a comparative analysis can be performed with all bioinformatic workflows 
# originally developed for microarray analyses.
library(limma)
RNAseq_voom = voom(RNAseq)$E
##transpose matrix to correlate genes in the following
#取前面5000个基因
WGCNA_matrix = t(RNAseq[order(apply(RNAseq,1,mad), decreasing = T)[1:5000],])


#2. 选择合适的阀值
#similarity measure between gene profiles: biweight midcorrelation
library(WGCNA)
s = abs(bicor(WGCNA_matrix))
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9 #建议值为0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
#calculation of adjacency matrix
beta = 3
a = s^beta
#dissimilarity measure
w = 1-a
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# 左边那张图表示，红线附近的值都是可以作为阈值的，R^2已经在0.9左右了。
# 如果想知道WGCNA建议的阈值，可以输入sft，列表的第一项值就是估计的最优值。。。
# 也就是后续代码里面的power值。

#分步构建网络
#create gene tree by average linkage hierarchical clustering 
geneTree = hclust(as.dist(w), method = 'average')
#module identification using dynamic tree cut algorithm
modules = cutreeDynamic(dendro = geneTree, distM = w, deepSplit = 4, 
                        pamRespectsDendro = FALSE,minClusterSize = 30)
#assign module colours
module.colours = labels2colors(modules)
#plot the dendrogram and corresponding colour bars underneath
plotDendroAndColors(geneTree, module.colours, 'Module colours', dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main='')

# The relation between the identified co-expression modules can be visualized by 
# a dendrogram of their eigengenes. The module eigengene is defined as the 
# first principal component of its expression matrix. It could be shown that the 
# module= eigengene is highly correlated with the gene that has the highest 
# intramodular connectivity.
library(ape)
#calculate eigengenes
MEs = moduleEigengenes(WGCNA_matrix, colors = module.colours, excludeGrey = FALSE)$eigengenes
#calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
#cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = 'average')
#plot the result with phytools package
par(mar=c(2,2,2,2))
plot.phylo(as.phylo(METree),type = 'fan',show.tip.label = FALSE, main='')
tiplabels(frame = 'circle',col='black', text=rep('',length(unique(modules))), 
          bg = levels(as.factor(module.colours)))

####？？？？到下面的需要临床clinical数据那一步就蒙逼了！！！
#Relation of co-expression modules to sample traits
#load clinical metadata. Make sure that patient barcodes are in the same format 
#create second expression matrix for which the detailed clinical data is available 
WGCNA_matrix2 = WGCNA_matrix[match(clinical$Name, rownames(WGCNA_matrix)),]
#CAVE: 1 sample of detailed clinical metadata is not in downloaded data (TCGA-GN-A269-01')
not.available = which(is.na(rownames(WGCNA_matrix2))==TRUE)
WGCNA_matrix2 = WGCNA_matrix2[-not.available,]
str(WGCNA_matrix2)
#hence it needs to be removed from clinical table for further analysis
clinical = clinical[-not.available,]






# 第三步：
# 筛选到合适的阈值后，我们就可以构建网络图了（也就是 模块图）。
# 官网给出了两种方法，一个是一步法：简单明了，一步出结果。
# 另一个是分步法：可以细调一系列的参数，以达到满意的结果。
# 官网也没说哪种适合什么条件，都可以用吧。具体参数我就不列了，
# 可以看这个博客：http://blog.sina.com.cn/s/blog_61f013b80101lcpr.html 讲比我详细多了。

# 也可以使用一步法网络构建：One-step network construction and module detection
net = blockwiseModules(WGCNA_matrix, power = 6, maxBlockSize = 5000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.45,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "AS-green-FPKM-TOM",
                       verbose = 3)
table(net$colors)
# table(net$colors)可以查看不同模块的颜色以及每个模块里的基因数目。
# 我为了减少模块数，所以mergeCutHeight用了0.45，这个根据自身需求来吧，
# 一般都是0.25就行了。然后下面的代码就可以展示出网络图了，如图2所示。

#4. 绘画结果展示
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# 然后我们需要把一些参数保留，以用于后续的分析，如下:
#5.结果保存
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]

# 第四步：
# 接下来的就是不同的策略的选择了。
# 方法1：最简单的就是对于每个模块进行GO以及KEGG富集，因为从理论上来说，
# 每个模块都应该都有生物学意义的，因此可以看每个模块的GO以及KEGG注释的结果，
# 看看都富集到哪个通路，再结合每个模块的特征值（查看MEs即可）来确定哪个模块的
# 基因作为后续研究的切入点。
# 方法2：如果再有前人研究或者自己的调查的结果，知道某几基因是自己研究方向的重要基因，
# 可以通过这几个基因在哪些或者某个模块，来确定接下来的选择哪个模块作为研究重点。
# 方法3：就像官网说的，结合模块与样本性状的相关性来分析，但一般都没有样本性状数据的，
# 我猜。。。
# 方法4：虽然我觉得有点扯，但是我还是会选择用，至少在论文里看到过。。
# 就是结合模块与样本间的相关系数。。
MEs0 <- moduleEigengenes(WGCNA_matrix, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
datTraits <- diag(56)
datTraits <- as.data.frame(datTraits)
rownames(datTraits) <- rownames(WGCNA_matrix)
colnames(datTraits) <- rownames(WGCNA_matrix)
nSamples = nrow(WGCNA_matrix)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
jpeg(file="labeledHeatmap.jpeg")
par(mar = c(3.5, 1.5, 1.5, 1))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
# y轴是模块名，x轴是样本名。。每个模块都跟每个样本有一个空格，
# 空格的上部分是模块与样本的相关系数（代码里面我构建了一个标量矩阵代表样本，
# 然后求了模块和样本的相关系数），空格的下部分是一个p值，我猜是代表相关的显著性吧
# （不知道理解对不对，但是反正都是挑p值越小的越好）。。由于我不知道这个测试数据的样
# 本信息各个代表咋样的实验，但是如果是有生物学重复实验的话，有人是这样做的：
# 在3个生物学重复中有2个样本在某个模块是显著相关的话，就挑选作为后续研究。。。。
# 有的人是这样的：关注与某一个或者几个样本，看看有哪个模块与这些样本是显著相关的。。。。
# 还有的人。。我也不知道了。。文献还没怎么看。。

# 第五步：调选完模块后，我们可以做个该模块的特征值所对应的柱状图和聚类图，
# 我这里就不做了（其实是我没查到WGCNA里有无该内置函数，反正别人是做了。。
# 实在找不到函数的话，只好自己画了）然后就可以做一张TOM图
#1. 可视化全部基因网络
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(WGCNA_matrix, power = 12)
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

# 第六步：
# 导出网络图到cytoscape，记得在代码里面修改自己所需要模块（颜色）。。
# 然后再筛选下txt表格里面的weight值（应该就是TOM值？），可以选择硬阈值筛，
# 还有种说法是用软阈值筛。。等我搞懂了再来修改。。。最后就导入cytoscape作图即可。。
#6. 导出网络到Cytoscape
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(WGCNA_matrix, power = 6)
# Read in the annotation file
# annot = read.csv(file = "GeneAnnotation.csv");
# Select modules需要修改，选择需要导出的模块颜色
modules = c("turquoise")
# Select module probes选择模块探测
probes = colnames(WGCNA_matrix)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("AS-green-FPKM-One-step-CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("AS-green-FPKM-One-step-CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
