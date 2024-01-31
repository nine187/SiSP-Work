library(DeepCC)
library(keras)

# get functional spectra from gene expression profiles and the geneSets will use MSigDBv7 by default. Add parameter `geneSets = "MSigDBvx"` (x = 5/6/7) or `geneSets = MSigDBr` (MSigDBr attained from R package `msigdbr` with function get_msigdbr() in DeepCC) for different version of MSigDB.
fs <- getFunctionalSpectra(eps, geneSets = "MSigDBv7")

# train DeepCC model
deepcc_model <- train_DeepCC_model(fs, labels)

# obtain deep features 
df <- get_DeepCC_features(deepcc_model, fs)