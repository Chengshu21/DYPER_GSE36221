##################################################################
# Gene Set Analysis of transcriptome data for human lung samples #
#################  Xie chengshu write in Oct,2019 ################
##################################################################
########################### GSVA analysis ########################
##################################################################

## load packages dependencies 
suppressPackageStartupMessages(library(GSVAdata));
suppressPackageStartupMessages(library(GSEABase));
suppressPackageStartupMessages(library(GSVA));

### this file could be downloaded form [MsigDB](http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C2).
KEGG_genesetcollection = readRDS("data/KEGG_genesetcollection.rds")

								
run.gsva = function(matrix, KEGG_genesetcollection){
  
    result.gsva = gsva(as.matrix(matrix), KEGG_genesetcollection, method = "gsva", 
                       mx.diff = FALSE, parallel.sz=1, abs.ranking = FALSE, verbose=TRUE)
    result.plage = gsva(as.matrix(matrix), KEGG_genesetcollection, method = "plage", 
                        mx.diff = FALSE, parallel.sz=1, abs.ranking = FALSE, verbose=TRUE)
    result.ssgsea = gsva(as.matrix(matrix), KEGG_genesetcollection, method = "ssgsea", 
                         mx.diff = FALSE, parallel.sz=1, abs.ranking = FALSE, verbose=TRUE)
    result.zscore = gsva(as.matrix(matrix), KEGG_genesetcollection, method = "zscore", 
                         mx.diff = FALSE, parallel.sz=1, abs.ranking = FALSE, verbose=TRUE)
    result = list("result.ssgsea" = result.ssgsea, "result.plage" = result.plage, 
                  "result.gsva" = result.gsva, "result.zscore" = result.zscore)
}

#### all gene set results
## result_gsva = lapply(data_new, function(x) result.gsva(x, KEGG_genesetcollection))
## names(result_gsva) = c("result.gsva.mon0", "result.gsva.mon6", "result.gsva.mon30")

### results under different time_point
gsva.mon0 = run.gsva(as.matrix(data_new$mon0[,which(colnames(data_new$mon0) %in% rownames(pdata_mon0))]), KEGG_genesetcollection)
gsva.mon6 = run.gsva(as.matrix(data_new$mon6[,which(colnames(data_new$mon6) %in% rownames(pdata_mon6))]), KEGG_genesetcollection)
gsva.mon30 = run.gsva(as.matrix(data_new$mon30[,which(colnames(data_new$mon30) %in% rownames(pdata_mon30))]), KEGG_genesetcollection)


##################################################################
######################### DE gene sets ###########################
##################################################################
run_limma_gsva = function(exprSet, pdata, contrast){
  
  if (contrast == "tre1-tre0"){
    exprSet = exprSet[,which(pdata[,40] != 2)]
    pdata = pdata[which(pdata[,40] != 2),]
    f = factor(pdata[,40], levels=c("0", "1"))
    design = model.matrix(~0+f)
    colnames(design) = c("tre0", "tre1")
    rownames(design) = rownames(pdata[which(pdata[,40] !=2 ), ]) 
    contrast.matrix = makeContrasts(tre1-tre0, levels=design)
  } 
  if (contrast == "tre2-tre0"){
    exprSet = exprSet[,which(pdata[,40] != 1)]
    pdata = pdata[which(pdata[,40] != 1),]
    f = factor(pdata[,40], levels=c("0", "2"))
    design = model.matrix(~0+f)
    colnames(design) = c("tre0", "tre2")
    rownames(design) = rownames(pdata[which(pdata[,40] !=1 ), ]) 
    contrast.matrix = makeContrasts(tre2-tre0, levels=design)
  }
  fit = lmFit(exprSet, design)
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2)
  results = topTable(fit2, adjust="BH", n = Inf)
  results$Term = rownames(results)  
  colnames(results) = c("logFC","AveExpr","t","P.value","adj.P.Val","B","Term" )
  results
  #write.table(results,paste0("limma_DEG.",contrast,".txt"),quote = F,sep = '\t',row.names = T)
}

### HERE ONLY USE "gsva" METHOD
### To check other method, change the result data:
### example in "ssgsea": run_limma(gsva.mon0$result.gsva, pdata_mon0, contrast = "tre1-tre0")
### mon 0 
DEgenesets_mon0_tre10 = run_limma_gsva(gsva.mon0$result.gsva, pdata_mon0, contrast = "tre1-tre0")
DEgenesets_mon0_tre20 = run_limma_gsva(gsva.mon0$result.gsva, pdata_mon0, contrast = "tre2-tre0")
### mon 6 
DEgenesets_mon6_tre10 = run_limma_gsva(gsva.mon6$result.gsva, pdata_mon6, contrast = "tre1-tre0")
DEgenesets_mon6_tre20 = run_limma_gsva(gsva.mon6$result.gsva, pdata_mon6, contrast = "tre2-tre0")
### mon 30 
DEgenesets_mon30_tre10 = run_limma_gsva(gsva.mon30$result.gsva, pdata_mon30, contrast = "tre1-tre0")
DEgenesets_mon30_tre20 = run_limma_gsva(gsva.mon30$result.gsva, pdata_mon30, contrast = "tre2-tre0")



##################################################################
############################# plot ###############################
##################################################################
suppressPackageStartupMessages(library(gplots));
suppressPackageStartupMessages(library(ggplot2));

## manage data and check the head DE gene sets
gsva.tre10_epmon0 = cbind(DEgenesets_mon0_tre10, time_point = rep(0, 186))[,c(4,7:8)]
gsva.highlight10_mon0 = head(gsva.tre10_epmon0[order(gsva.tre10_epmon0$P.value, decreasing=FALSE),])
gsva.tre10_epmon6 = cbind(DEgenesets_mon6_tre10, time_point = rep(6, 186))[,c(4,7:8)]
gsva.highlight10_mon6 = head(gsva.tre10_epmon6[order(gsva.tre10_epmon6$P.value, decreasing=FALSE),])
gsva.tre10_epmon30 = cbind(DEgenesets_mon30_tre10, time_point = rep(30, 186))[,c(4,7:8)]
gsva.highlight10_mon30 = head(gsva.tre10_epmon30[order(gsva.tre10_epmon30$P.value, decreasing=FALSE),])

gsva.tre20_epmon0 = cbind(DEgenesets_mon0_tre20, time_point = rep(0, 186))[,c(4,7:8)]
gsva.highlight20_mon0 = head(gsva.tre20_epmon0[order(gsva.tre20_epmon0$P.value,decreasing=FALSE),] )
gsva.tre20_epmon6 = cbind(DEgenesets_mon6_tre20, time_point = rep(6, 186))[,c(4,7:8)]
gsva.highlight20_mon6 = head(gsva.tre20_epmon6[order(gsva.tre20_epmon6$P.value,decreasing=FALSE),] )
gsva.tre20_epmon30 = cbind(DEgenesets_mon30_tre20, time_point = rep(30, 186))[,c(4,7:8)]
gsva.highlight20_mon30 = head(gsva.tre20_epmon30[order(gsva.tre20_epmon30$P.value,decreasing=FALSE),] )


gsva.tre10 = rbind(gsva.tre10_epmon0, gsva.tre10_epmon6,gsva.tre10_epmon30)
gsva.tre20 = rbind(gsva.tre20_epmon0, gsva.tre20_epmon6, gsva.tre20_epmon30)
# write.csv(gsva.tre10, "E:/gsva.tre10.csv")
# write.csv(gsva.tre20, "E:/gsva.tre20.csv")


### manage and check 
### which pathway consider to be significant
gsva.high10.pathway0 = head(gsva.highlight10_mon0$Term)
gsva.high10.pathway6 = head(gsva.highlight10_mon6$Term)
gsva.high10.pathway30 = head(gsva.highlight10_mon30$Term)

gsva.high20.pathway0 = head(gsva.highlight20_mon0$Term)
gsva.high20.pathway6 = head(gsva.highlight20_mon6$Term)
gsva.high20.pathway30 = head(gsva.highlight20_mon30$Term)
## gsva.high10.pathway6 have more pathways related to COPD
## "KEGG_NITROGEN_METABOLISM", "KEGG_ALDOSTERONE_REGULATED_SODIUM_REABSORPTION", "KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS"              "KEGG_HOMOLOGOUS_RECOMBINATION"                 
## "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY", "KEGG_PPAR_SIGNALING_PATHWAY" 
## color6 
## "Blue", "Green",  "Black", "Purple", "DeepPink","DarkRed"
## "lightgrey" means gene sets are out of high.pathway, that is "others"

### source("data_mutate.R")
gsva.tre10_temp = gsva.tre10
gsva.tre20_temp = gsva.tre20
gsva.t1 = data_mutate(gsva.tre10_temp, gsva.high10.pathway6)
gsva.t2 = data_mutate(gsva.tre20_temp, gsva.high10.pathway6)

color = c("Blue", "Green",  "Black", "Purple", "DeepPink","LightGray","DarkRed")
gsva.plot.spaghetti = function(data, color){
  plot <- ggplot(data,aes(x = time_point, y = P.value, group = Term, 
                          colour = color, size = factor(highlight))) +
    geom_line() +   
    ### Caution:
    ### the color of each pathway will be  changed randomly
    scale_color_manual(values = color) +
    scale_size_manual(values = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 1.5)) +
    ### ggtitle("Trend of pathways during 0~30 months under treatment1_treatment0") +
    ### labels of each line:
    ### In fact, it's better not to make them into one function
    ### because we should adjust the parameters according to the plot
    geom_label( x=20, y=0.6, label="KEGG_NITROGEN_METABOLISM", size=2, color="Blue") +
    geom_label( x=3, y=0.2, label="KEGG_PPAR_SIGNALING_PATHWAY", size=2, color="DarkRed") +
    geom_label( x=25, y=0.5, label="KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY", size=2, color="DeepPink") +
    geom_label( x=6, y=0.6, label="KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS", size=2, color="Black") +
    geom_label( x=20, y=0.2, label="KEGG_ALDOSTERONE_REGULATED_SODIUM_REABSORPTION", size=2, color="Green") +
    geom_label( x=25, y=0.4, label="KEGG_HOMOLOGOUS_RECOMBINATION", size=2, color="Purple") +
    theme(legend.position="none", 
          axis.title.x=element_text(size=15,face="bold",hjust=0.5),
          axis.title.y=element_text(size=15,hjust=0.5),
          axis.text.x=element_text(size=14),
          axis.text.y=element_text(size=14))
  plot 
}

jpeg(file = "plots/gsva.naive1.jpg")
#pdf(file = "plots/gsva.naive1.pdf")
gsva.p1 = gsva.plot.spaghetti(gsva.t1,color)
gsva.p1
dev.off()
jpeg(file = "plots/gsva.naive2.jpg")
#pdf(file = "plots/gsva.naive2.pdf")
gsva.p2 = gsva.plot.spaghetti(gsva.t2,color)
gsva.p2
dev.off()


### save.image("RData/SSonly.RData")
