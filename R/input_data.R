################
# Get GEO data(GSE36221)
# AUTHOR: Chengshu Xie 
# Created: Sep.2019
# Last update: Dec, 2019
# GEO: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36221
# PLATFORM: GPL6244
################

################################################
###### Convert probe IDs into gene symbols #####
################################################
suppressPackageStartupMessages(library(hugene10sttranscriptcluster.db));
suppressPackageStartupMessages(library(R.utils));
suppressPackageStartupMessages(library(GEOquery));
suppressPackageStartupMessages(library(tidyverse));

get_GEO_data = function(geo_accession, bioc_annotation_package){  

	suppressPackageStartupMessages(library(GEOquery));
	suppressPackageStartupMessages(library(tidyverse));
	
	geo_data = getGEO(geo_accession, destdir=".",getGPL = F);
	exprSet = as.data.frame(exprs(geo_data[[1]]));
	pdata = as.data.frame(pData(geo_data[[1]]));
	exprSet = exprSet[which(colnames(exprSet) %in% rownames(pdata))];
	exprSet = exprSet[, order(colnames(exprSet))];
	pdata = pdata[order(rownames(pdata)), ];
	
	#### convert probe ids into 
	exprSet$probe_id =  as.character(rownames(exprSet));
	
	if (bioc_annotation_package == "hugene10sttranscriptcluster.db") {
		probe2symbol_df = toTable(get("hugene10sttranscriptclusterSYMBOL"))};
  
    newmatrixdata = exprSet %>% 
                     inner_join(probe2symbol_df,by="probe_id") %>% #merge informatio of probe_id
					 select(-probe_id) %>% #remove extra column         
					 select(symbol, everything()) %>% #re-arrange
					 mutate(rowMean = rowMeans(.[grep("GSM", names(.))])) %>% #calculate means
					 arrange(desc(rowMean))  %>% #order the means
					 distinct(symbol,.keep_all = T) %>% # keep the first symbol information
					 select(-rowMean) %>% #remove the new colmun
					 tibble::column_to_rownames(colnames(.)[1]); #convert the first column into rownames and remove it  
    newmatrixdata = as.data.frame(newmatrixdata);
	result = list("pdata" = pdata, "exprsdata" = newmatrixdata)
	result
}


origin_data = get_GEO_data("GSE36221", "hugene10sttranscriptcluster.db")

################################################
################# tidy data ####################
################################################
pdata = origin_data$pdata
pdata_de = subset(pdata, pdata[,40] != 3)
GSE36221 = origin_data$exprsdata[, which(colnames(origin_data$exprsdata) %in% rownames(pdata_de))]
GSE36221 = GSE36221[, order(colnames(GSE36221))]

pdata_mon0 = pdata_de[which(pdata_de[,39] == 0),]
pdata_mon6 = pdata_de[which(pdata_de[,39] == 6),]
pdata_mon30 = pdata_de[which(pdata_de[,39] == 30),]
data_new = list()
data_new$mon0 = GSE36221[, colnames(GSE36221) %in% rownames(pdata_mon0)]
data_new$mon6 = GSE36221[, colnames(GSE36221) %in% rownames(pdata_mon6)]
data_new$mon30 = GSE36221[, colnames(GSE36221) %in% rownames(pdata_mon30)]

## save.image("RData/input_data.RData")