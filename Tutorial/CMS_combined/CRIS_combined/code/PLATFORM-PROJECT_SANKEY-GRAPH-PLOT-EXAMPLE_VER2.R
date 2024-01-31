library(magrittr)

# CMS AND CRIS CLASSIFICATION
CMSpred <- data.frame("Sample" = c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10"), 
                      "CMS classification" = c(
  "CMS1", "CMS3", "CMS2", "CMS1", "CMS4", "CMS4", "CMS2", NA, NA, "CMS2"))
CRISpred <- data.frame("Sample" = c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10"), 
                       "CRIS classification" = c(
  "CRIS-B", "CRIS-C", "CRIS-A", "CRIS-D", "CRIS-E", "CRIS-E", "CRIS-B", "CRIS-B", "CRIS-B", NA))
CMSpred[which(is.na(CMSpred$CMS.classification)), "CMS.classification"] <-  "NA"
CRISpred[which(is.na(CRISpred$CRIS.classification)), "CRIS.classification"] <-  "NA"

# LOAD TREATMENT BASED ON CMS AND CRIS CLASSIFICATION
CMSdrug <- read.csv("data/CMS-drug.csv")
CMSdrug[which(is.na(CMSdrug$Subtype)), ] <- "NA"
CRISdrug <- read.csv("data/CRIS-drug.csv")
CRISdrug[which(is.na(CRISdrug$Treatment_list)), "Treatment_list"] <- "NA"
CRISdrug[which(is.na(CRISdrug$Subtype)), ] <- "NA"

# DEFINE NODE
CMSdrug.node <- data.frame("name" = unique(CMSdrug$Treatment_list), 
                           "node" = seq(1:length(unique(CMSdrug$Treatment_list))) - 1)
CMSpred.node <- data.frame("name" = levels(factor(CMSpred$CMS.classification)), 
                           "node" = seq(1:length(levels(factor(CMSpred$CMS.classification)))) + 
                             max(CMSdrug.node$node))
CRISpred.node <- data.frame("name" = levels(factor(CRISpred$CRIS.classification)), 
                            "node" = seq(1:length(levels(factor(CRISpred$CRIS.classification)))) + 
                              max(CMSpred.node$node))
CRISdrug.node <- data.frame("name" = unique(CRISdrug$Treatment_list), 
                            "node" = seq(1:length(unique(CRISdrug$Treatment_list))) + 
                              max(CRISpred.node$node))
Node <- data.frame(rbind(CMSdrug.node, CMSpred.node, CRISpred.node, CRISdrug.node))

CMSPred.list <- vector()
CMSDrug.list <- vector()
CMSSamp.list <- vector()
for (i in 1:length(CMSpred$CMS.classification)) {
  subtype.list <- CMSdrug[which(CMSdrug$Subtype == CMSpred$CMS.classification[i]), "Subtype"]
  treat.list <- CMSdrug[which(CMSdrug$Subtype == CMSpred$CMS.classification[i]), "Treatment_list"]
  sample.list <- CMSpred$Sample[i]
  
  CMSPred.list <- append(CMSPred.list, subtype.list)
  CMSDrug.list <- append(CMSDrug.list, treat.list)
  CMSSamp.list <- append(CMSSamp.list, rep(sample.list, length(subtype.list)))
}
CMSPred.Drug <- data.frame("Sample" = CMSSamp.list, 
                           "Subtype" = CMSPred.list, 
                           "Treatment_list" = CMSDrug.list) 

CMSPred.pos <- vector()
CMSDrug.pos <- vector()
for (i in 1:length(CMSPred.Drug$Sample)) {
  treat.pos <- CMSdrug.node$node[which(CMSdrug.node$name %in% CMSPred.Drug$Treatment_list[i])]
  subtype.pos <- CMSpred.node$node[which(CMSpred.node$name %in% CMSPred.Drug$Subtype[i])]
  CMSDrug.pos <- append(CMSDrug.pos, treat.pos)
  CMSPred.pos <- append(CMSPred.pos, subtype.pos)
}
CMSPred.Drug <- CMSPred.Drug %>% 
  dplyr::mutate("Subtype.pos" = CMSPred.pos) %>% 
  dplyr::mutate("Treatment_list.pos" = CMSDrug.pos)

CRISPred.list <- vector()
CRISDrug.list <- vector()
CRISSamp.list <- vector()
for (i in 1:length(CRISpred$CRIS.classification)) {
  sutype.list <- CRISdrug[which(CRISdrug$Subtype == CRISpred$CRIS.classification[i]), "Subtype"]
  treat.list <- CRISdrug[which(CRISdrug$Subtype == CRISpred$CRIS.classification[i]), "Treatment_list"]
  sample.list <- CRISpred$Sample[i]
  CRISPred.list <- append(CRISPred.list, sutype.list)
  CRISDrug.list <- append(CRISDrug.list, treat.list)
  CRISSamp.list <- append(CRISSamp.list, sample.list)
}
CRISPred.Drug <- data.frame("Sample" = CRISSamp.list, "Subtybe" = CRISPred.list, 
                            "Treatment_list" = CRISDrug.list)
#CRISPred.Drug_1 <- data.frame()
#for (i in 1:length(CMSPred.Drug$Sample)) {
 # CRISPred.Drug_1[i, c("Sample", "Subtybe", "Treatment_list")] <- 
  #  CRISPred.Drug[which(CRISPred.Drug$Sample == CMSPred.Drug$Sample[i]), ]
#}

CRISPred.pos <- vector()
CRISDrug.pos <- vector()
for (i in 1:length(CRISPred.Drug$Sample)) {
  treat.pos <- CRISdrug.node$node[which(CRISdrug.node$name %in% CRISPred.Drug$Treatment_list[i])]
  subtype.pos <- CRISpred.node$node[which(CRISpred.node$name %in% CRISPred.Drug$Subtybe[i])]
  CRISDrug.pos <- append(CRISDrug.pos, treat.pos)
  CRISPred.pos <- append(CRISPred.pos, subtype.pos)
}
CRISPred.Drug <- CRISPred.Drug %>% 
  dplyr::mutate("Subtype.pos" = CRISPred.pos) %>% 
  dplyr::mutate("Treatment_list.pos" = CRISDrug.pos)

# DEFINE LINK
Link <- data.frame("source" = c(CMSPred.Drug$Treatment_list, CMSPred.Drug$Subtype, 
                                CRISPred.Drug$Subtybe), 
                   "source.pos" = c(CMSPred.Drug$Treatment_list.pos, CMSPred.Drug$Subtype.pos, 
                                    CRISPred.Drug$Subtype.pos),
                   "target" = c(CMSPred.Drug$Subtype, CRISPred.Drug$Subtybe, 
                                CRISPred.Drug$Treatment_list), 
                   "target.pos" = c(CMSPred.Drug$Subtype.pos, CRISPred.Drug$Subtype.pos, 
                                    CRISPred.Drug$Treatment_list.pos),
                   "value" = 1)

# Combine data in a list
CRCsubtypes <- list("links" = Link, "nodes" = Node)

# Add a factor column to list in the nodes 
Group <- as.factor(Node$name)
CRCsubtypes$nodes$group <-  Group
CRCsubtypes$nodes$name_1 <- Node$name
CRCsubtypes$nodes$name[1:length(CRCsubtypes$nodes$name)] = " "

# Call library for plotting the network
library(networkD3)
library(htmlwidgets)
#node_color <- 'd3.scaleOrdinal() .domain(["CMS1", "CMS2", "CMS3", "CMS4", "CRIS-A", "CRIS-B", "CRIS-C", #"CRIS-D", "CRIS-E", "Not Defined"]) .range(["#FFA9A9", "#D7BEFF", "#9FE2BF", "#FFE493", "#D7263D", "#F46036", "#2E294E", "#1B998B", "#00B7E0", "#bfc0c0"])'

node_color <- 'd3.scaleOrdinal() .domain(["CMS1", "CMS2", "CMS3", "CMS4", "CRIS-A", "CRIS-B", "CRIS-C", "CRIS-D", "CRIS-E", "Response:Cetuximab", "Response:FOLFOX", "Response:FOLFOX,Cetuximab", "Unknown", "Unresponsiveness:Cetuximab", "Unresponsiveness:5 FU/Response:Immunotherapy", "Unresponsiveness:FOLFOX,Cetuximab/Response:FOLFIRI", "NA"]) .range(["#FFA9A9", "#D7BEFF", "#9FE2BF", "#FFE493", "#D7263D", "#F46036", "#2E294E", "#1B998B", "#00B7E0", "#2E294E", "#9FE2BF", "#D7BEFF", "#F46036", "#D7263D", "#FFA9A9", "#FFE493",  "#bfc0c0"])'

# Plot
plot <- sankeyNetwork(Links = CRCsubtypes$links, 
                      Nodes = CRCsubtypes$nodes, 
                      Source = 'source.pos', Target = 'target.pos', Value = 'value', 
                      NodeID = 'name', fontSize = 14, 
                      nodeWidth = 20, iterations = 0,
                      colourScale = node_color, units = 'TWh',
                      NodeGroup = "group", 
                      fontFamily = "arial", 
                      #width = 900, height = 600,
                      #margin = list("left" = 200, "right" = 10)
                      ) 
plot[["x"]][["nodes"]][["name_1"]] <- CRCsubtypes[["nodes"]][["name_1"]]
plot
htmlwidgets::onRender(plot, '
                      function(el, x) {
                        d3.select(el).selectAll(".node text").each(function(d){
                          var arr, val, anc
                          arr = d.name_1 ;
                          arr = arr.split("/");
                          val = d3.select(this).attr("x");
                          anc = d3.select(this).attr("text-anchor"); 
                          for(i = 0; i < arr.length; i++) {
                            d3.select(this).append("tspan")
                                .text(arr[i])
                                .attr("dy", i ? "1.2em" : 0)
                                .attr("x", val)
                                .attr("text-anchor", anc)
                                .attr("class", "tspan" + i);
                          }
                        })
                      }')

# https://stackoverflow.com/questions/44700596/put-line-break-in-node-labels-in-networkd3-sankey-diagram
# https://stackoverflow.com/questions/72404933/how-to-plot-sankey-graph-with-r-networkd3-values-and-percentage-below-each-node
