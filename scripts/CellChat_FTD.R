############################################################
##### Cell-cell communication using CellChat - Disease #####
############################################################

setwd("")

library(anndata)
library(Seurat)
library(CellChat)
library(dplyr)

###### Load data ######
data <- read_h5ad("FTD_adata.h5ad")

data$X
counts = data$layers["counts"]
count_matrix = t(as.matrix(counts))

### Normalize count data
norm_count_matrix = NormalizeData(count_matrix, scale.factor = 10000, do.log = TRUE)

###### Load metadata ######
metadata = data.frame(data$obs)

# See if rownames of metdata match column names of count matrix
identical(colnames(norm_count_matrix), rownames(metadata)) 

unique(metadata$annotation_04142026) # This is the cell annotation column, which will be used in 'group.by' below. What groups you want to find ligands and receptors in
colnames(metadata)[colnames(metadata) == "SampleID"] <- "samples"

### Create CellChat object
cellchat <- createCellChat(object = norm_count_matrix, meta = metadata, group.by = "annotation_04142026")

cellchat <- setIdent(cellchat, ident.use = "annotation_04142026") # set "labels" as default cell identity

levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

cellchat <- setIdent(cellchat, ident.use = "annotation_04142026") # set "labels" as default cell identity

levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

### Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

### Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

### Inference of cell-cell communication network
# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)

### Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

### Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file="cellChat_FTD.rds")
