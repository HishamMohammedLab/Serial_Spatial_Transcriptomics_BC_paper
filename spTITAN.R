library(Seurat)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(stats)
library(pdist)
library(parallel)
library(lda)
library(text2vec)
library(LICORS)
library(TITAN)

SO <- readRDS("data/Patient_A_SOwTimepoints.rds")
PA_raw_spat <- readRDS("data/PA_raw_spatLocs.rds")

rownames(PA_raw_spat) <- PA_raw_spat$cell_ID

slide2_cells <- PA_raw_spat$cell_ID[PA_raw_spat[,1] < 17] #Look only at cells on one physical slide

PA_slide2 <- subset(SO, cells = slide2_cells)

SO <- PA_slide2

SO[["image"]] <- new(Class = "SlideSeq", assay = "Spatial", coordinates = PA_raw_spat[colnames(SO),2:1])

pos_df <- PA_raw_spat[colnames(SO),1:2]

set.seed(8)

km <- kmeans(pos_df, 200, nstart = 25)   #Run kmeans clustering to identify centroids

pos_df$points <- "data"

center_df <- as.data.frame(km$centers)
center_df$points <- "centroid"
rownames(center_df) <- paste0("centroid_", 1:nrow(center_df))

total_df <- rbind(pos_df, center_df)

# Extended Data Figure I

pdf("plots/Centroids_on_PAslide2.pdf")
ggplot(total_df, mapping=aes(y=-(sdimx), x=sdimy, color=points)) + geom_point() + theme_classic()
dev.off()

# Calculate distance for each cell to centroids

dist_df <- as.matrix(pdist(pos_df[,1:2], center_df[,1:2]))
rownames(dist_df) <- rownames(pos_df)
colnames(dist_df) <- rownames(center_df)

inv_dists <- max(dist_df) - dist_df

test <- t(inv_dists)

for (i in 1:ncol(test)) {
  temp <- test[,i]
  cutoff <- quantile(temp, probs = 0.95) 
  temp[which(temp < cutoff)] <- 0
  test[,i] <- temp
}

test <- test/10 #scale distance values

Object.sparse <- GetAssayData(SO, layer = "counts",assay = "RNA")
final_mat <- rbind(Object.sparse, test)

temp_SO <- CreateSeuratObject(counts=final_mat)

# Run TITAN's underlying topic model on concatenated gene/distance matrix
# Very similar to runLDA function, with slight tweaks to accounts for the distance matrix

assayName <- "RNA"
normalizationMethod <- "CLR"
ntopics <- seq(10,200, by = 10)
alpha <- 50
beta <- 0.1
varFeatures <- 5000
iterations <- 500
burnin <- 250
seed.number <- 8
parallel <- T
outDir <- "plots/test_centroidModels_5quantile"
cores <- 10

Object        <- NormalizeData(temp_SO, assay = "RNA", normalization.method = "CLR")
#Object        <- FindVariableFeatures(Object, assay = assayName, nfeatures = 5000)
Object.sparse <- GetAssayData(Object, layer = "data",assay = "RNA")
#Object.sparse <- Object.sparse[VariableFeatures(Object, assay = assayName),]

#convert data into the proper input format for lda.collapsed.gibbs.sampler
data.use      <- Matrix::Matrix(Object.sparse, sparse = T)

data.use      <- data.use * 10
data.use      <- round(data.use)
data.use      <- Matrix::Matrix(data.use, sparse = TRUE)
sumMat        <- Matrix::summary(data.use)
cellList      <- split(as.integer(data.use@i),
                       sumMat$j)
ValueList     <- split(as.integer(sumMat$x),
                       sumMat$j
)
cellList      <- mapply(rbind, cellList, ValueList, SIMPLIFY=F)
Genes         <- rownames(data.use)
cellList      <- lapply(cellList, function(x) {colnames(x) <- Genes[x[1,]+1];x})

#Run model
model_maker <- function(topics) {
  selected.Model <- lda.collapsed.gibbs.sampler(
    cellList,
    topics,
    Genes,
    num.iterations=iterations,
    alpha=alpha,
    eta=beta,
    compute.log.likelihood=TRUE,
    burnin=burnin)[-1]
  if (parallel) {
    if (!dir.exists(outDir)) {
      dir.create(outDir)
    }
    saveRDS(selected.Model, paste0(outDir, "/Model_", as.character(topics), "topics.rds"))
  } else {
    return(selected.Model)
  }
}

if (parallel) {
  mclapply(ntopics, model_maker, mc.cores = cores)
} else {
  Model <- model_maker(ntopics)
  return(Model)
}

# Grab the data to calculate elbow plot and determine best model

files <- list.files(path = outDir, pattern = "Model_")

# Get model input data
Object <- temp_SO

#Normalize and extract the gene expression data from the Seurat Object
Object        <- NormalizeData(Object, assay = assayName, normalization.method = "CLR")
#Object        <- FindVariableFeatures(Object, assay = assayName, nfeatures = varFeatures)
Object.sparse <- GetAssayData(Object, layer = "data",assay = assayName)
#Object.sparse <- Object.sparse[VariableFeatures(Object, assay = assayName),]

#convert data into the proper input format for lda.collapsed.gibbs.sampler
data.use      <- Matrix::Matrix(Object.sparse, sparse = T)

data.use <- data.use * 10
data.use <- round(data.use)

#initialize necessary variables
perp_list     <- NULL
topic_numbers <- NULL
RPC           <- NULL
files         <- files[order(nchar(files), files)]

for (model_file in files) {
  topic_num     <- as.numeric(gsub("[^0-9]+([0-9]+).*", "\\1", model_file))
  topic_numbers <- c(topic_numbers, topic_num)
  model         <- readRDS(paste0(outDir, "/", model_file))
  
  #extract document-term matrix
  docterMat     <- t(as.matrix(data.use))
  docterMat     <- as(docterMat, "sparseMatrix")
  
  #calculate topic word distribution
  topworddist   <- normalize(model$topics, byrow = T)
  
  #calculate document topic distribution
  doctopdist    <- normalize(t(model$document_sums), byrow = T)
  
  #calculate perpelexity
  perp          <- perplexity(docterMat, topworddist, doctopdist)
  perp_list     <- c(perp_list, perp)
  
  #calculate RPC (rate of perplexity change)
  if (length(perp_list) > 1) {
    RPC_temp <- abs((perp_list[length(perp_list)] - perp_list[length(perp_list) - 1]) / (topic_numbers[length(topic_numbers)]        - topic_numbers[length(topic_numbers) - 1]))
    RPC      <- c(RPC, RPC_temp)
  }
}

#build plot dataframe and create ggplot object
plot_df           <- as.data.frame(cbind(topic_numbers[-1], RPC))
colnames(plot_df) <- c("Topics", "RPC")
p                 <- ggplot(data = plot_df, aes(x = Topics, y = RPC, group = 1)) + geom_line() + geom_point()

pdf(paste0(outDir,"/test_elbowPlot.pdf"))
p
dev.off()

#Use 80 topics

model <- readRDS("data/test_centroidModels_5quantile/Model_80topics.rds")

model_topGenes <- TopTopicGenes(model, n = 100)

SO@meta.data <- SO@meta.data[,-c(9:48)]

SO <- addTopicsToSeuratObject(model, SO)

# Extended data Figure J

pdf("plots/PA_slide2_newSpatTITAN_5quantile_Topic74_spatialFeaturePlot.pdf")
SpatialFeaturePlot(SO, features = "Topic_74", stroke = NA, min.cutoff = 0) + scale_fill_gradient(low = "gray90", high = "darkblue")
dev.off()

pdf("PA_slide2_Topic46_IFITM1_spatialExpression.pdf")
SpatialFeaturePlot(SO, features = "IFITM1", stroke = NA, min.cutoff = 0, max.cutoff = 2) + scale_fill_gradient(low = "gray90", high = "forestgreen")
dev.off()
