##################################################################
# Gene Set Analysis of transcriptome data for human lung samples #
#################  Xie chengshu write in Oct,2019 ################
##################################################################
######################## clustering test #########################
##################################################################
#### prepare data
patient_ID = pdata_de[,38];
newID = names(which(table(patient_ID) == 3));
data_conv_gsva = function(data, pdata, cond, KEGG_genesetcollection){
  if (cond == "tre1-tre0"){
    npdata = pdata[which(pdata[,40] !=2 & pdata[,38] %in% newID ),];
    gsva_input_data = data[,which(colnames(data) %in% rownames(npdata))];
    colnames(gsva_input_data) = npdata[order(rownames(npdata)), ][,38]  ;
  } else if (cond == "tre2-tre0"){
    npdata = pdata[which(pdata[,40] !=1 & pdata[,38] %in% newID ),];
    gsva_input_data = data[,which(colnames(data) %in% rownames(npdata))];
    colnames(gsva_input_data) = npdata[order(rownames(npdata)), ][,38];
  } else if (is.NULL(cond)){
    print("Please add augument value: cond");
  }
  
  suppressPackageStartupMessages(library(GSVAdata));
  suppressPackageStartupMessages(library(GSEABase));
  suppressPackageStartupMessages(library(GSVA));
  
  result.gsva = gsva(as.matrix(gsva_input_data), KEGG_genesetcollection, method = "gsva", 
                     mx.diff = FALSE, parallel.sz=1, abs.ranking = FALSE, verbose=TRUE);
  result.gsva
}


#### get gsva results
tre10_mon0.gsvares = data_conv_gsva(data_new$mon0, pdata_mon0, cond = "tre1-tre0", KEGG_genesetcollection);
tre10_mon6.gsvares = data_conv_gsva(data_new$mon6, pdata_mon6, cond = "tre1-tre0", KEGG_genesetcollection);
tre10_mon30.gsvares = data_conv_gsva(data_new$mon30, pdata_mon30, cond = "tre1-tre0", KEGG_genesetcollection);

tre20_mon0.gsvares = data_conv_gsva(data_new$mon0, pdata_mon0, cond = "tre2-tre0", KEGG_genesetcollection);
tre20_mon6.gsvares = data_conv_gsva(data_new$mon6, pdata_mon6, cond = "tre2-tre0", KEGG_genesetcollection);
tre20_mon30.gsvares = data_conv_gsva(data_new$mon30, pdata_mon30, cond = "tre2-tre0", KEGG_genesetcollection);

clu.gsva.data = list("c1" = tre10_mon0.gsvares, "c2" = tre10_mon6.gsvares, "c3" = tre10_mon30.gsvares,
                     "c4" = tre20_mon0.gsvares, "c5" = tre20_mon6.gsvares, "c6" = tre20_mon30.gsvares)

##### function to get cluster results
clu.res = function(data){
  t = dist(scale(t(data)), method = "euclidean");
  t.clust_ward = hclust(t, method = "ward.D2");
  ## plot(t.clust_ward, hang =-1, main = "ward.D2", xlab = "tre10_mon0")
  ## rect.hclust(t.clust_ward, k=4)
  label <- cutree(t.clust_ward,k=4);
  
  clust.name = list();
  for (i in 1:4){
    temp = c();
    for(j in 1:length(label)){
      if (label[j] == i){
        temp[[j]] = names(label[j]);
        temp = na.omit(temp);
        clust.name[[i]] = as.list(temp);
        names(clust.name[[i]]) = as.character(temp);
      }
    }
  }
  names(clust.name) = c("cluster1", "cluster2", "cluster3", "cluster4");
  return(clust.name)
}

#### cluster results of patients
clu.gsva.res = lapply(clu.gsva.data, clu.res)
names(clu.gsva.res) = c("C1", "C2", "C3","C4","C5","C6");

##### function to compute jaccard index
jaccard_index =  function(data1, data2, cond){
  
  simiID.num = c();
  diffID.num = c();
  jaccard.num = data.frame();
  
  for(i in 1:4){
    for(j in 1:4){
    simiID.num[i] = length(intersect(names(data1[[i]]), names(data2[[j]])));
    diffID.num[i] = length(union(names(data1[[i]]), names(data2[[j]])));
    jaccard.num[i,j] = simiID.num[i]/diffID.num[i];
    }
  }
  rownames(jaccard.num) = names(data1);
  colnames(jaccard.num) = names(data2);
  jaccard.num$sum1 = apply(jaccard.num, 1, sum);
  jaccard.num = rbind(jaccard.num, "sum2" = t(apply(jaccard.num, 2, sum)));
  jaccard.re = list("jaccard.num" = jaccard.num,
                    "Row VS Col" = paste("data1 VS data2"));
  return(jaccard.re)
}

###### 
jaccard.res_C1_C2 = jaccard_index(clu.gsva.res$C1, clu.gsva.res$C2);
jaccard.res_C3_C2 = jaccard_index(clu.gsva.res$C3, clu.gsva.res$C2);
jaccard.res_C4_C5 = jaccard_index(clu.gsva.res$C4, clu.gsva.res$C5);
jaccard.res_C6_C5 = jaccard_index(clu.gsva.res$C6, clu.gsva.res$C5)


save.image("E:/COPD_GSE36221/results/clustering test.RData")
