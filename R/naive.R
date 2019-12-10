##################################################################
# Gene Set Analysis of transcriptome data for human lung samples #
#################  Xie chengshu write in Oct,2019 ################
##################################################################
########################## naive analysis ########################
##################################################################
### load packages
suppressPackageStartupMessages(library(limma));
suppressPackageStartupMessages(library(enrichR));
suppressPackageStartupMessages(library(tidyverse));
suppressPackageStartupMessages(library(ggplot2));
suppressPackageStartupMessages(library(ggvis));
suppressPackageStartupMessages(library(gridExtra));  ## required to arrange ggplot2 plots in a grid

### functions to get DE genes
run_limma = function(exprSet, pdata, contrast){
  
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
  DEgenes = results[results$P.Value <=0.05, ]
  DEgenes
  #write.table(results,paste0("limma_DEG.",contrast,".txt"),quote = F,sep = '\t',row.names = T)
}


##################################################################
########################### DE genes #############################
##################################################################
### mon 0 
DEgenes_mon0_tre10 = run_limma(data_new$mon0, pdata_mon0, contrast = "tre1-tre0")
DEgenes_mon0_tre20 = run_limma(data_new$mon0, pdata_mon0, contrast = "tre2-tre0")
### mon 6 
DEgenes_mon6_tre10 = run_limma(data_new$mon6, pdata_mon6, contrast = "tre1-tre0")
DEgenes_mon6_tre20 = run_limma(data_new$mon6, pdata_mon6, contrast = "tre2-tre0")
### mon 30 
DEgenes_mon30_tre10 = run_limma(data_new$mon30, pdata_mon30, contrast = "tre1-tre0")
DEgenes_mon30_tre20 = run_limma(data_new$mon30, pdata_mon30, contrast = "tre2-tre0")

#####################################
############ DE pathways ############
#####################################
dbs <- listEnrichrDbs()
dbs = "KEGG_2019_Human"
### pathway results
enriched_pathway_mon0_tre10 = enrichr(rownames(DEgenes_mon0_tre10), dbs)
enriched_pathway_mon0_tre20 = enrichr(rownames(DEgenes_mon0_tre20), dbs)
enriched_pathway_mon6_tre10 = enrichr(rownames(DEgenes_mon6_tre10), dbs)
enriched_pathway_mon6_tre20 = enrichr(rownames(DEgenes_mon6_tre20), dbs)
enriched_pathway_mon30_tre10 = enrichr(rownames(DEgenes_mon30_tre10), dbs)
enriched_pathway_mon30_tre20 = enrichr(rownames(DEgenes_mon30_tre20), dbs)

## manage data and check the head DE gene sets
tre10_epmon0 = cbind(as.data.frame(enriched_pathway_mon0_tre10$KEGG_2019_Human), time_point = rep(0, 219))[,c(1,3,10)]
highlight10_mon0 = head(tre10_epmon0[order(tre10_epmon0$P.value,decreasing=FALSE),] )
tre10_epmon6 = cbind(as.data.frame(enriched_pathway_mon6_tre10$KEGG_2019_Human), time_point = rep(6, 294))[,c(1,3,10)]
highlight10_mon6 = head(tre10_epmon6[order(tre10_epmon6$P.value,decreasing=FALSE),])
tre10_epmon30 = cbind(as.data.frame(enriched_pathway_mon30_tre10$KEGG_2019_Human), time_point = rep(30, 263))[,c(1,3,10)]
highlight10_mon30 = head(tre10_epmon30[order(tre10_epmon30$P.value,decreasing=FALSE),])

tre20_epmon0 = cbind(as.data.frame(enriched_pathway_mon0_tre20$KEGG_2019_Human), time_point = rep(0, 308))[,c(1,3,10)]
highlight20_mon0 = head(tre20_epmon0[order(tre20_epmon0$P.value,decreasing=FALSE),] )
tre20_epmon6 = cbind(as.data.frame(enriched_pathway_mon6_tre20$KEGG_2019_Human), time_point = rep(6, 303))[,c(1,3,10)]
highlight20_mon6 = head(tre20_epmon6[order(tre20_epmon6$P.value,decreasing=FALSE),] )
tre20_epmon30 = cbind(as.data.frame(enriched_pathway_mon30_tre20$KEGG_2019_Human), time_point = rep(30, 305))[,c(1,3,10)]
highlight20_mon30 = head(tre20_epmon30[order(tre20_epmon30$P.value,decreasing=FALSE),] )


tre10 = rbind(tre10_epmon0, tre10_epmon6,tre10_epmon30)
tre20 = rbind(tre20_epmon0, tre20_epmon6, tre20_epmon30)
## write.csv(tre10, "~/tre10.csv")
## write.csv(tre20, "~/tre20.csv")
###heatmap data
###tre10_t = gather(tre10_tt, "time_point", "P.value", -Term)
###tre10_t$time_point = as.numeric(tre10_t$time_point)


### manage and check 
### which pathway consider to be significant
high10.pathway0 = head(highlight10_mon0$Term)
high10.pathway6 = head(highlight10_mon6$Term)
high10.pathway30 = head(highlight10_mon30$Term)

high20.pathway0 = head(highlight20_mon0$Term)
high20.pathway6 = head(highlight20_mon6$Term)
high20.pathway30 = head(highlight20_mon30$Term)
## high10.pathway6 have more pathways related to COPD
## "p53 signaling pathway", "Cell cycle", "Systemic lupus erythematosus",
## "Pancreatic cancer", "Viral carcinogenesis", "Melanoma"  
## color6 
## "Blue", "Green",  "Black", "Purple", "DeepPink","DarkRed"
## "lightgrey" means gene sets are out of high.pathway, that is "others"
data_mutate = function(data, high.pathway){
  
  tr1 = data[-which(data$Term %in% high.pathway),];
  tr1$highlight = rep("others",nrow(tr1));
  tr1$color = rep("LightGray",nrow(tr1));
  re = c();
  res = c();
  
  
  for (i in 1:6){
    color <- switch(i,
                    "Blue", 
                    "Green",
                    "Black",
                    "Purple",
                    "DeepPink",
                    "DarkRed")
    re[[i]] = data[which(data$Term == high.pathway[[i]]),] %>%
      mutate( highlight=ifelse(Term==high.pathway[[i]], high.pathway[[i]])) %>%
      mutate( color=ifelse(Term==high.pathway[[i]], color))
      
  }
  res = rbind(re[[1]], re[[2]],re[[3]],re[[4]],re[[5]],re[[6]], tr1)
  
}
tre10_temp = tre10
tre20_temp = tre20
t1 = data_mutate(tre10_temp, high10.pathway6)
t2 = data_mutate(tre20_temp, high10.pathway6)

color = c("Blue", "Green",  "Black", "Purple", "DeepPink","LightGray","DarkRed")
plot.spaghetti = function(data, color){
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
    geom_label( x=20, y=0.6, label="P53 signaling pathway", size=4, color="Blue") +
    geom_label( x=3, y=0.2, label="Melanoma", size=4, color="DarkRed") +
    geom_label( x=25, y=0.5, label="Viral carcinogenesis", size=4, color="DeepPink") +
    geom_label( x=6, y=0.6, label="Systemic lupus erythematosus", size=4, color="Black") +
    geom_label( x=20, y=0.2, label="Cell cycle", size=4, color="Green") +
    geom_label( x=25, y=0.4, label="Pancreatic cancer", size=4, color="Purple") +
    theme(legend.position="none", 
          axis.title.x=element_text(size=15,face="bold",hjust=0.5),
          axis.title.y=element_text(size=15,hjust=0.5),
          axis.text.x=element_text(size=14),
          axis.text.y=element_text(size=14))
  plot 
}

jpeg(file = "plots/naive1.jpg")
#pdf(file = "plots/naive1.pdf")
p1 = plot.spaghetti(t1,color)
p1
dev.off()
jpeg(file = "plots/naive2.jpg")
# pdf(file = "plots/naive2.pdf")
p2 = plot.spaghetti(t2,color)
p2
dev.off()


save.image("RData/naive.RData")
