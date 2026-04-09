library(Seurat)
library(TITAN)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(bluster)
library(rstatix)
library(effectsize)
library(reshape2)
library(GGally)
library(viridis)
#library(factoextra)
#library(FactoMineR)
library(umap)
library(ggrepel)
library(igraph)

SO <- readRDS("data/cosmx_filteredCells_noLiver_200counts_newAnno.rds") #Liver cells were removed with help of pathologist annotations

cancer_meta <- read.table("data/cosmx_Patients_meta_clones.tsv.gz", sep = "\t", header = T, row.names = 1) #Clones were identified by CNA analysis

SO$Clone <- cancer_meta[rownames(SO@meta.data), "Clone"]

# Keeping only cancer cells

SO <- subset(SO, Clone %in% c("PatA_Clone1", "PatA_Clone2", "PatA_Clone3", "PatB_Clone1", "PatB_Clone2",
                              "PatB_Clone3", "PatC_Clone1", "PatC_Clone2", "PatC_Clone3", "PatD_Clone1",
                              "PatD_Clone2", "PatD_Clone3"))

#runLDA(SO, ntopics = seq(10,100, by=10), parallel = T, outDir = "allCosmx_cancer_Models", cores = 10)

#pdf("allCosmx_cancer_Models/TITAN_elbowPlot.pdf")
#LDAelbowPlot("allCosmx_cancer_Models/", SO) + theme_bw() + theme(panel.border = element_blank())
#dev.off()

model <- readRDS("data/allCosmx_cancer_Models/Model_30topics.rds")

top_genes <- TopTopicGenes(model, ngenes = 50)

gene_scores <- GeneScores(model)

SO <- addTopicsToSeuratObject(model, SO)

PA_SO <- readRDS("data/Patient_A_SOwTimepoints.rds")
PB_SO <- readRDS("data/Patient_B_SOwTimepoints.rds")
PC_SO <- readRDS("data/Patient_C_SOwTimepoints.rds")
PD_SO <- readRDS("data/Patient_D_SOwTimepoints.rds")

allSO_wTimepoints <- merge(PA_SO, c(PB_SO, PC_SO, PD_SO))

SO$Timepoint <- allSO_wTimepoints@meta.data[colnames(SO), "Timepoint"]  # Adding cell timepoint information

colors23 <- c("#58041e", "#880d22", "#db0844", "#b7210e", "#ce5118",
              "#ff943d", "#f9ba0e", "#a19026", "#6c6747", "#425229",
              "#62da26", "#659b5c", "#41665e", "#284650", "#6fa1ca",
              "#7774d2", "#5d3ac0", "#49296e", "#6b4d73", "#d19adc",
              "#913291") 

HBCA_colors <- colors23[-c(1,3,5,7,9,11,13,14,17,19,20)]

HeatmapTopic <- function(Object,
                    topics,
                    AnnoVector,
                    AnnoName,
                    clusterTopics = F, cols) {
  #Create a dataframe with the annotation information and corresponding colors
  anno_col           <- data.frame(row.names = colnames(Object),
                         Column1=AnnoVector)
  colnames(anno_col) <- AnnoName
  num_colors         <- length(unique(anno_col[,1]))
  anno_colors        <- cols
  names(anno_colors) <- sort(unique(anno_col[,1]))
  anno_colors        <- list(Cluster = anno_colors)
  names(anno_colors) <- AnnoName
  #Add annotation color information to topics
  topics <- data.matrix(topics[order(anno_col[,1]),])
  #plot Heatmap
  if (clusterTopics == F) {
    p1 <- pheatmap(topics,
                   hclustfun = function(x) hclust(x, method="ward.D2"),
                   scale = "row",
                   cluster_cols = F,
                   cluster_rows = F,show_rownames = F,
                   col=colorRampPalette(rev(brewer.pal(11, "RdBu"))[c(1:4,8:11)])(256),
                   annotation_row = anno_col,
                   annotation_names_row = T,
                   annotation_colors = anno_colors,
                   cex=1)
  } else {
    p1 <- pheatmap(topics,
                   hclustfun = function(x) hclust(x, method="ward.D2"),
                   scale = "row",
                   cluster_cols = T,
                   cluster_rows = F,show_rownames = F,
                   col=colorRampPalette(rev(brewer.pal(11, "RdBu"))[c(1:4,8:11)])(256),
                   annotation_row = anno_col,
                   annotation_names_row = T,
                   annotation_colors = anno_colors,
                   cex=1)
  }
  return(p1)
}

SO <- RunUMAP(SO, reduction = "lda", dims = 1:30, reduction.name = "TITANumap")

# Figure 1B (All Malignant Epithelial Cells Section)

png("plots/allCosmx_cancer_30topics_UMAP_byTimepoint.png")
DimPlot(SO, reduction = "TITANumap", group.by = "Timepoint")
dev.off()

png("plots/allCosmx_cancer_30topics_UMAP_byPatient.png")
DimPlot(SO, reduction = "TITANumap", group.by = "Patient")
dev.off()

topics <- SO@meta.data[,24:53]
anno_col           <- data.frame(row.names = colnames(SO), Column1=SO$Patient, Column2=SO$Timepoint)
colnames(anno_col) <- c("Patient", "Timepoint")
num_colors         <- 4
anno_colors        <- gg_color_hue(num_colors)[c(1,2,3,4)]
names(anno_colors) <- sort(unique(anno_col[,2]))
anno_colors2       <- HBCA_colors[c(1,3,5,8)]
names(anno_colors2) <- sort(unique(anno_col[,1]))
anno_colors        <- list(Timpoint = anno_colors, Cell_type=anno_colors2)
names(anno_colors) <- c("Timepoint", "Patient")
topics <- data.matrix(topics[order(anno_col[,1], anno_col[,2]),])

p1 <- pheatmap(topics,
                   hclustfun = function(x) hclust(x, method="ward.D2"),
                   scale = "row",
                   cluster_cols = T,
                   cluster_rows = F,show_rownames = F,
                   col=colorRampPalette(rev(brewer.pal(11, "RdBu"))[c(1:4,8:11)])(256),
                   annotation_row = anno_col,
                   annotation_names_row = T,
                   annotation_colors = anno_colors,
                   cex=1)

saveRDS(SO, "allCosmx_cancer_wTopics.rds")

meta_df <- SO@meta.data[,c(24:53,6,54)]

#write.table(meta_df, "allCosmx_TITAN_cancerOnly_topics_andAnno_Fig4.tsv", quote = F, sep = "\t")

plot_df <- meta_df %>% group_by(Patient, Timepoint) %>% summarise_all("mean")

anno_col <- data.frame(row.names = rownames(plot_df), Patient = plot_df$Patient, Timepoint = plot_df$Timepoint)
anno_colors <- gg_color_hue(4)[c(1,2,3,4)]
names(anno_colors) <- sort(unique(anno_col[,2]))
anno_colors2 <- HBCA_colors[c(1,3,5,8)]
names(anno_colors2) <- sort(unique(anno_col[,1]))
anno_colors <- list(Timepoint = anno_colors, HBCA_anno=anno_colors2)
#plot_df <- data.matrix(plot_df[order(anno_col[,1], anno_col[2,]),])

hm_mat <- as.matrix(as.data.frame(plot_df)[,3:32])
rownames(hm_mat) <- rownames(plot_df)

p1 <- pheatmap(hm_mat,
                   hclustfun = function(x) hclust(x, method="ward.D2"),
                   scale = "row",
                   cluster_cols = T,
                   cluster_rows = F,show_rownames = F,
                   col=colorRampPalette(rev(brewer.pal(11, "RdBu"))[c(1:4,8:11)])(256),
                   annotation_row = anno_col,
                   annotation_names_row = T,
                   annotation_colors = anno_colors,
                   cex=1)

#saveRDS(hm_mat, "allCosmx_cancer_condensedTopicMat.rds")

topics <- SO@meta.data[,24:53]

cor_mat <- cor(topics)

model <- readRDS("data/allCosmx_cancer_Models/Model_30topics.rds")

topic_gene_mat <- GeneScores(model)

topic.umap <- umap(t(topic_gene_mat), n_neighbors=10, min.dist = 0.5, spread = 0.5)

colnames(topic.umap$layout) <- c("Dim1", "Dim2")

umap_plot_df <- as.data.frame(topic.umap$layout)
umap_plot_df$Topic <- rownames(umap_plot_df)

ggplot(umap_plot_df, mapping = aes(x=Dim1, y=Dim2)) + geom_point() + theme_classic()

# Clustering topics into Topic clusters

make.knn.graph <- function(D,k) {
    dist <- as.matrix(dist(D))
    edges <- mat.or.vec(0,2)
    for (i in 1:nrow(dist)) {
        matches <- setdiff(order(dist[i,], decreasing = F)[1:(k+1)],i)
        edges <- rbind(edges, cbind(rep(i,k), matches))
        edges <- rbind(edges, cbind(matches, rep(i,k)))
    }
    graph <- graph_from_edgelist(edges, directed = F)
    V(graph)$frame.color <- NA
    set.seed(1)
    g.layout <- layout_with_fr(graph)
    return(list(graph=graph, layout=g.layout))
}

set.seed(10)

topic_graph <- make.knn.graph(t(topic_gene_mat), k=5)
topic_clusters <- cluster_louvain(topic_graph$graph, resolution = 0.9)

table(topic_clusters$memberships[1,])

cluster_assignment <- as.factor(topic_clusters$memberships[1,])
names(cluster_assignment) <- colnames(topic_gene_mat)

raw_topics <- model$topics

threshold <- .01
nms       <- as.character()
value     <- as.numeric()
n_clusters = length(unique(topic_clusters$memberships[1,]))

test <- topic_gene_mat

for(n in 1:n_clusters) {
    for(i in 1:ncol(t(test)[row.names(t(test)) %in% names(cluster_assignment[cluster_assignment == n]),])){
        print(i)
        print(n)
     if(mean(t(test)[row.names(t(test)) %in% names(cluster_assignment[cluster_assignment == n]),i]) > threshold){
         nms <- c(nms, colnames(t(test)[row.names(t(test)) %in% names(cluster_assignment[cluster_assignment == n]),])[i])
         value <- c(value,mean(t(test)[row.names(t(test)) %in% names(cluster_assignment[cluster_assignment == n]),i]))
     }
        print(length(value))
        print(length(nms))
    }
    if(n == 1) {
        result <- data.frame(GENE = nms, VALUE = value, CLUSTER = rep(n,length(nms)))
        rm(nms)
        nms <- as.character()
        rm(value)
        value <- as.numeric()
        }
    if(n >  1) {
        result <- rbind(result, data.frame(GENE = nms, VALUE = value, CLUSTER = rep(n,length(nms))))
        rm(nms)
        nms <- as.character()
        rm(value)
        value <- as.numeric()
        }
}

table(result$CLUSTER)

top.markers <- as.data.frame(result %>% group_by(CLUSTER) %>% top_n(n=10, wt = VALUE))

#write.table(top.markers, "allCosmx_cancer_topTopicCluster_genes.tsv", sep = "\t", row.names = F, quote = F)

plot_df <- data.frame(Topics = sub("Topic_", "", rownames(topic.umap$layout)), Dim1 = topic.umap$layout[,1], Dim2 = topic.umap$layout[,2], Cluster = as.character(topic_clusters$memberships[1,]))

plot_df$oldCluster <- plot_df$Cluster
plot_df["Topic_17", "Cluster"] <- "6" # Creating new cluster based off of biological pathway investigation

plot_df <- plot_df[-(which(rownames(plot_df) %in% c("Topic_11", "Topic_13", "Topic_18"))),]

#saveRDS(plot_df, "allCosmx_cancerTITAN_UMAPofTopic_Coordinates.rds")
#plot_df <- readRDS("allCosmx_cancerTITAN_UMAPofTopic_Coordinates.rds")

#saveRDS(hm_mat, "allCosmx_cancerTITAN_HMmat.rds")
#hm_mat <- readRDS("allCosmx_cancerTITAN_HMmat.rds")

ggplot(plot_df, mapping = aes(x=Dim1, y=Dim2, color=Cluster)) + geom_point() + theme_classic() + scale_color_manual(values=colors23[c(1,3,5,7,9,11,13,15,17)])

anno_clusters = data.frame(row.names=rownames(plot_df), Cluster = as.character(plot_df$Cluster))

hm_mat <- hm_mat[,rownames(plot_df)[order(plot_df$Cluster)]]
anno_clusters <- data.frame(row.names = rownames(plot_df)[order(plot_df$Cluster)], Cluster = as.character(plot_df$Cluster)[order(plot_df$Cluster)])

plot_df <- meta_df %>% group_by(Patient, Timepoint) %>% summarise_all("mean")

anno_col <- data.frame(row.names = rownames(plot_df), Patient = plot_df$Patient, Timepoint = plot_df$Timepoint)
anno_colors <- c(rgb(1,0,0,0.5),
                 rgb(0,1,0,0.5),
                 rgb(0,0,1,0.5),
                 rgb(1,1,0,0.5))
names(anno_colors) <- sort(unique(anno_col[,1]))
anno_colors2 <- HBCA_colors[c(1,3,5,8)]
names(anno_colors2) <- sort(unique(anno_col[,2]))
anno_colors3 <- colors23[c(1,3,5,7,9,11,13,15,17,19)]
names(anno_colors3) <- sort(unique(anno_clusters[,1]))
anno_colors <- list(Timepoint = anno_colors2, Patient=anno_colors, Cluster=anno_colors3)

p1 <- pheatmap(hm_mat,
               hclustfun = function(x) hclust(x, method="ward.D2"),
               scale = "row",
               cluster_cols = F,
               cluster_rows = F,
               show_rownames = F,
               col=colorRampPalette(rev(brewer.pal(11, "RdBu"))[c(1:4,8:11)])(256),
               annotation_row = anno_col,
               annotation_col = anno_clusters,
               annotation_names_row = T,
               annotation_names_col = F,
               annotation_colors = anno_colors,
               fontsize = 8,
               gaps_col = c(6, 13, 19, 22, 26),
               cex=1)

# Extended Figure 5A

pdf("plots/allCosmx_cancerTITAN_condensedHeatmap_wTopicClusterLabels_wGaps.pdf")
p1
dev.off()

#saveRDS(SO, "allCosmx_cancer_wTopics.rds")
#SO <- readRDS("allCosmx_cancer_wTopics.rds")

cor_mat <- cor(SO@meta.data[,24:53])

anno_colors3 <- colors23[c(1,3,5,7,9,11,13,15,17,19)]
names(anno_colors3) <- sort(unique(anno_clusters[,1]))
anno_colors <- list(Cluster=anno_colors3)

cor_mat <- cor_mat[rownames(anno_clusters)[order(anno_clusters$Cluster)], rownames(anno_clusters)[order(anno_clusters$Cluster)]]

paletteLength <- 256

myBreaks <- c(seq(-1, 0, length.out = ceiling(paletteLength/2) + 1),
              seq(max(cor_mat)/paletteLength,max(cor_mat), length.out = floor(paletteLength/2)))

p1 <- pheatmap(cor_mat,
               hclustfun = function(x) hclust(x, method="ward.D2"),
#               scale = "row",
               cluster_cols = F,
               cluster_rows = F,
               show_rownames = T,
               col=colorRampPalette(rev(brewer.pal(11, "RdBu"))[c(1:4,8:11)])(256),
               annotation_row = anno_clusters,
               annotation_col = anno_clusters,
               annotation_names_row = T,
               annotation_names_col = F,
               annotation_colors = anno_colors,
               fontsize = 8,
               breaks = myBreaks,
#               gaps_col = c(15, 27, 32, 38, 44, 51, 58),
               cex=1)


clus1_topics <- rownames(anno_clusters)[which(anno_clusters$Cluster == "1")]
clus2_topics <- rownames(anno_clusters)[which(anno_clusters$Cluster == "2")]
clus3_topics <- rownames(anno_clusters)[which(anno_clusters$Cluster == "3")]
clus4_topics <- rownames(anno_clusters)[which(anno_clusters$Cluster == "4")]
clus5_topics <- rownames(anno_clusters)[which(anno_clusters$Cluster == "5")]
clus6_topics <- rownames(anno_clusters)[which(anno_clusters$Cluster == "6")]


SO@meta.data$Topic_Cluster1 <- rowMeans(SO@meta.data[,clus1_topics])
SO@meta.data$Topic_Cluster2 <- rowMeans(SO@meta.data[,clus2_topics])
SO@meta.data$Topic_Cluster3 <- rowMeans(SO@meta.data[,clus3_topics])
SO@meta.data$Topic_Cluster4 <- rowMeans(SO@meta.data[,clus4_topics])
SO@meta.data$Topic_Cluster5 <- rowMeans(SO@meta.data[,clus5_topics])
SO@meta.data$Topic_Cluster6 <- SO@meta.data[,clus6_topics]

cluster_radar <- SO@meta.data[,c(54:60,6)]
cluster_radar[which(cluster_radar < 0, arr.ind = T)] <- 0

head(cluster_radar)

saveRDS(cluster_radar, "data/allCosmx_cancer_radarDF.rds")

