rm(list=ls())
graphics.off()

library(ggplot2)
library(networkD3) #data visualization
#devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey) #data visualzation
library(dplyr) #data manipulation
library(htmlwidgets) #data visualization

#load the CMS results
CRIS <- read.delim2(file = "data/CRIS_COADREAD_prediction_result.xls", header = TRUE)
DeepCC <- read.csv(file = "data/DeepCC_result.csv", header = TRUE)

#clean the dataset
# Function to remove all characters after the third dot
remove_after_third_dot <- function(x) {
  parts <- unlist(strsplit(x, "\\.")) # Split the string into parts
  result <- paste(parts[1:3], collapse = ".") # Concatenate the first three parts
  return(result)
}
# apply the function to the column
CRIS$sample.names <- sapply(CRIS$sample.names, remove_after_third_dot)

#print(CRIS$sample.names)
#print(DeepCC$Row.names)

#match
merged_df <- merge(CRIS, DeepCC, by.x = "sample.names", by.y = "Row.names", all.x = TRUE)

#piechart table
CRIS <- as.data.frame(table(merged_df$predict.label2))
CMS <- as.data.frame(table(merged_df$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples))
DeepCC <- as.data.frame(table(merged_df$y))
#piechart
#CRIS
examplePlot <- ggplot(CRIS, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  geom_text(aes(label = Freq),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(title = "CRIS Data Distribution")

#CMS
examplePlot2 <- ggplot(CMS, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  geom_text(aes(label = Freq),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(title = "CMS Data Distribution")

#DeepCC
examplePlot3 <- ggplot(DeepCC, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  geom_text(aes(label = Freq),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(title = "DeepCC Data Distribution")

merged_df <- merged_df[, -c(3:7)]

#create a new column for sankey df
selected_columns <- c(4,2)
sankey_df <- as.data.frame(merged_df[, selected_columns])

#rename
colnames(sankey_df)[1] <- "DeepCC"
colnames(sankey_df)[2] <- "CRIS"

#NA
sankey_df[which(is.na(sankey_df$DeepCC)), "DeepCC"] <-  "NA"
table(sankey_df$DeepCC)
table(sankey_df$CRIS)

# assign numerical values to 'DeepCC_Column' based on the mapping
sankey_df <- mutate(sankey_df, DeepCC_val = case_when(
  DeepCC == "CMS1" ~ 0,
  DeepCC == "CMS2" ~ 1,
  DeepCC == "CMS3" ~ 2,
  DeepCC == "CMS4" ~ 3,
  DeepCC == "NA" ~ 4
))

deepcc_mapping <- c("CMS1" = 0, "CMS2" = 1, "CMS3" = 2, "CMS4" = 3, "NA" = 4)

sankey_df <- mutate(sankey_df, CRIS_val = case_when(
  CRIS == "CRIS-A" ~ 5,
  CRIS == "CRIS-B" ~ 6,
  CRIS == "CRIS-C" ~ 7,
  CRIS == "CRIS-D" ~ 8,
  CRIS == "CRIS-E" ~ 9,
  TRUE ~ NA_real_
))

sankey_df$assigned_val <- rep(1, nrow(sankey_df))

#rename col names
colnames(sankey_df) <- c("source", "target", "source.pos","target.pos","value")

###My dataset###
# CMS AND CRIS CLASSIFICATION
CMSpred <- data.frame("Sample" = c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9"), 
                      "CMS classification" = c(
                        "CMS1", "CMS3", "CMS2", "CMS1", "CMS4", "CMS4", "CMS2", NA, "CMS2"))
CRISpred <- data.frame("Sample" = c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9"), 
                       "CRIS classification" = c(
                         "CRIS-B", "CRIS-C", "CRIS-A", "CRIS-D", "CRIS-E", "CRIS-E", "CRIS-B", "CRIS-B", "CRIS-B"))
CMSpred[which(is.na(CMSpred$CMS.classification)), "CMS.classification"] <-  "NA"
CRISpred[which(is.na(CRISpred$CRIS.classification)), "CRIS.classification"] <-  "NA"
CMSpred.node <- data.frame("name" = levels(factor(CMSpred$CMS.classification)), 
                           "node" = c(0:4))
CRISpred.node <- data.frame("name" = levels(factor(CRISpred$CRIS.classification)), 
                            "node" = c(5:9))
Node <- data.frame(rbind(CMSpred.node, CRISpred.node))

# DEFINE LINK
Link <- data.frame("source" = c(unique(CMSpred$CMS.classification)), 
                   "source.pos" = c(0:4),
                   "target" = c(unique(CRISpred$CRIS.classification)), 
                   "target.pos" = c(5:9),
                   "value" = 1)
CRCsubtypes <- list("links" = sankey_df, "nodes" = Node)

Group <- as.factor(Node$name)
CRCsubtypes$nodes$group <-  Group

#plot
plot <- sankeyNetwork(Links = CRCsubtypes$links, 
                      Nodes = CRCsubtypes$nodes, 
                      Source = 'source.pos', Target = 'target.pos', Value = 'value', 
                      NodeID = 'name', fontSize = 14, 
                      nodeWidth = 20, iterations = 0,
                      units = 'TWh',
                      NodeGroup = "group", 
                      fontFamily = "arial")
                      #width = 900, height = 600,
                      #margin = list("left" = 200, "right" = 10)
plot

#heatmap
COADREAD <- as.matrix(sapply(COADREAD, as.numeric))
heatmap(COADREAD)