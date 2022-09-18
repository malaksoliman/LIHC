if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')

devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)

ens <- rownames(data)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = "ENSMBL")
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(data), names(symbols))]
rownames(data) <- symbols
keep <- !is.na(rownames(data))
data <- data[keep,]


res.df= as.data.frame(res)
res.symbols <- mapIds(org.Hs.eg.db, keys = rownames(res.df),
                  column = c('SYMBOL'), keytype = 'ENSEMBL')

library(org.Hs.eg.db)
res$symbols <- mapIds(org.Hs.eg.db, keys = rownames(res.df),
                  column = c('SYMBOL'), keytype = " ENSMBL")
EnhancedVolcano(res.df,
                lab = res.df$symbols,
                x = 'log2FoldChange',
                y = 'pvalue')

EnhancedVolcano(res.df,
                lab = res.df$symbol,
               x= 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.01,
                FCcutoff= 2)
EnhancedVolcano(res.df,
                lab = res.df$symbols,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = '',
                pCutoff = 0.01,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 6.0)
selected= c("PKM","IL6","HIF1A","HK2","SLC2A1","GPD1","ENO2","GCK","PFKFB3","PFKP","SLC2A2","GLS2","CAMK2B","ENO3")
EnhancedVolcano(res.df,
                lab = res.df$symbols,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.01,
                FCcutoff= 2,
                selectLab = selected)

EnhancedVolcano(res.df,
                lab = res.df$symbols,
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = selected,
                xlim = c(-20, 20),
                ylim = c(0,30),
                title = '',
                pCutoff = 0.01,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 2.0,
                boxedLabels = TRUE,
                colAlpha = 1,
                cutoffLineType = 'blank',
                cutoffLineCol = 'black',
                cutoffLineWidth = 0.8,
                hline = c(10e-20,
                          10e-20 * 10e-30,
                          10e-20 * 10e-60,
                          10e-20 * 10e-90),
                hlineCol = c('pink', 'hotpink', 'purple', 'black'),
                hlineType = c('solid', 'longdash', 'dotdash', 'dotted'),
                hlineWidth = c(1.0, 1.5, 2.0, 2.5),
                gridlines.major = FALSE,
                gridlines.minor = FALSE)
p1

p1 +
  ggplot2::coord_cartesian(xlim=c(-6, 6)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 1))
 
p1

EnhancedVolcano(res.df,
                lab = res.df$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.01,
                FCcutoff = 2,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-value & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0)

EnhancedVolcano(res.df,
                lab = res.df$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = "HIF1A",
                ylim = c(0, 50),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.01,
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')

lab_italics <- paste0("italic('", res.df$symbols, "')")
selectLab_italics = paste0(
  "italic('",
  c("PKM","IL6","HIF1A","HK2","SLC2A1","GPD1","ENO2","GCK","PFKFB3","PFKP","SLC2A2","GLS2","CAMK2B","ENO3"),
  "')")

EnhancedVolcano(res,
                lab = lab_italics,
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = selectLab_italics,
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.01,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                col = c('black', 'pink', 'purple', 'red3'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + coord_flip()
ens <- rownames(data)

library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(data), names(symbols))]
rownames(data) <- symbols
keep <- !is.na(rownames(data))
airway <- data[keep,]

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  res$log2FoldChange < -2, 'royalblue',
  ifelse(res$log2FoldChange > 2, 'gold',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'gold'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'low'

EnhancedVolcano(res.df,
                lab = res.df$symbols,
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = selected,
                
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Custom colour over-ride',
                pCutoff = 0.01,
                FCcutoff = 2.0,
                pointSize = 3.5,
                labSize = 4.5,
                shape = c(6, 4, 2, 11),
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'left',
                legendLabSize = 15,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                arrowheads = FALSE,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black')




p1 <- EnhancedVolcano(res.df,
                      lab = res.df$symbols,
                      x = "log2FoldChange",
                      y = "pvalue",
                      pCutoff = 0.01,
                      FCcutoff = 2,
                      ylim = c(0, -log10(10e-12)),
                      pointSize = c(ifelse(res$log2FoldChange>2, 8, 1)),
                      labSize = 6.0,
                      selectLab = selected,
                      shape = c(6, 6, 19, 16),
                      title = "DESeq2 results",
                      subtitle = "interesting genes",
                      caption = bquote(~Log[2]~ "fold change cutoff, 2; p-value cutoff, 10e-4"),
                      legendPosition = "right",
                      legendLabSize = 14,
                      colAlpha = 0.9,
                      colGradient = c('red3', 'royalblue'),
                      drawConnectors = TRUE,
                      hline = c(10e-8),
                      widthConnectors = 0.5)

p1

p1 +
  ggplot2::coord_cartesian(xlim=c(-6, 6)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 1))



# Define the genes of interest.
goi <- selected2
stopifnot(all(goi %in% names(dds)))
goi<- c("ENSG00000058404", "ENSG00000067057", "ENSG00000067225")

library(dplyr)
library(tidyr)
library(ggplot2)
library(DESeq2)
library(airway)
tcounts <- t(log2((counts(dds.run[goi, ], normalized=FALSE, replaced=FALSE)+.5))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression,(ncol(.)-length(goi)+1):ncol(.))
box<- tcounts %>% 
  select(Row.names,  gene, expression) %>% 
 
  knitr::kable()
selected2<- c("ENSG00000067225", "ENSG00000136244", "ENSG00000100644", "ENSG00000159399")

tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

ggplot(tcounts, aes(expression, fill=fake)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y") + 
  labs(x="Dexamethasone treatment", 
       y="Expression (log normalized counts)", 
       fill="(Some made \nup variable)", 
       title="Top Results")

p <-plotCounts(dds, gene="ENSG00000067225", intgroup="Pheno")

library(tidyverse)
library(hrbrthemes)
library(viridis)
data %>%
  ggplot(p, aes(x=Pheno, y=count, fill=Pheno)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Basic boxplot") +
  xlab("")

