############################################################
##### Cell-cell communication using CellChat - Merged ######
############################################################

setwd("")

library(CellChat)
library(patchwork)
library(ggplot2)
library(svglite)

##### Merge cellchat objects for comparison - https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html#load-cellchat-object-of-each-dataset-and-merge-them-together #####
cellchat.Control <- readRDS("cellChat_Control.rds")
cellchat.FTD <- readRDS("cellChat_FTD.rds")
object.list <- list(Control = cellchat.Control, FTD = cellchat.FTD)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Export the merged CellChat object and the list of the two separate objects for later use
save(object.list, file = "cellchat_object.list_FTD_vs_Control.RData")
save(cellchat, file = "cellchat_merged_FTD_vs_Control.RData")

# Look at number of interactions and strength of interactions between the two datasets (FTD and Control)
gg1 <- compareInteractions(cellchat, show.legend = F) 
gg1$mapping$fill = gg1$mapping$x
gg1 = gg1 + scale_fill_manual(values = c("Control" = "#9999FF", 
                                         "FTD" = "#669900"))

gg2 <- compareInteractions(cellchat, show.legend = F, measure = "weight")
gg2$mapping$fill = gg2$mapping$x
gg2 = gg2 + scale_fill_manual(values = c("Control" = "#9999FF", 
                                         "FTD" = "#669900"))

gg1 + gg2

##### Differential interaction strength among different cell populations across two datasets ###
#Where red or blue colored edges represent increased or decreased signaling in the second dataset compared to the first one.
#change to measure = "count" for differential number of interactions

# To make comparable across all plots below
weight.max <- getMaxWeight(object.list, 
                           slot.name = "net",
                           attribute = "weight")

netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c("Control", "FTD"), sources.use = "InNeu SOX14+", targets.use = NULL, measure = "weight", edge.weight.max = weight.max, label.edge = TRUE) 
netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c("Control", "FTD"), sources.use = "InNeu SOX14-", targets.use = NULL, measure = "weight", edge.weight.max = weight.max, label.edge = TRUE) 
netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c("Control", "FTD"), sources.use = "ExNeu RNF220+", targets.use = NULL, measure = "weight", edge.weight.max = weight.max, label.edge = TRUE) 
netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c("Control", "FTD"), sources.use = "ExNeu RNF220-", targets.use = NULL, measure = "weight", edge.weight.max = weight.max, label.edge = TRUE) 

### Heatmap ###
cell_order <- c("ExNeu RNF220+", "ExNeu RNF220-", 
                "InNeu SOX14-", "InNeu SOX14+", 
                "Lymphocyte", "BVC", "CP",
                "OPC", "OL", "Astrocyte", "Microglia", "PVM")

svglite("heatmap_cellchat_FTDvControl_strength.svg",
        width = 5, height = 5,
        system_fonts = list(sans = "Arial"))
netVisual_heatmap(cellchat, 
                  comparison = c("Control", "FTD"), 
                  title.name = "FTD vs. FTD-Control - Interaction strength", 
                  measure = "weight", # Change to "count" for differential number of interactions
                  col.show = cell_order,
                  row.show = cell_order) 
dev.off()

### Specific cell types in pathway ###
pdf("GABA_InNeu.pdf", width = 12, height = 4)

cells_of_interest <- c("InNeu SOX14-", "InNeu SOX14+")
pathways.show <- c("GABA-A") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:2) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]), sources.use = cells_of_interest)
} 

dev.off()





