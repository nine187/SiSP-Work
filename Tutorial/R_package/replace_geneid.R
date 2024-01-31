library(Biobase)
x.entrez <- Biobase::exprs(crcTCGAsubset)
# rownames are entrez
head(rownames(x.entrez))
# translate to gene symbols
x.symbol <- replaceGeneId(x.entrez, id.in="entrez", id.out="symbol")
head(rownames(x.symbol))
x.entrez2 <- replaceGeneId(x.symbol, id.in="symbol", id.out="entrez")
# translations are not cycle consistent
table(rownames(x.entrez2) == rownames(x.entrez))
# matrix values are not changed
all(x.entrez == x.entrez2)