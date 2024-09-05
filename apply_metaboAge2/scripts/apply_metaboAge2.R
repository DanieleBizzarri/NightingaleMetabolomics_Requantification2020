#This script is done to apply the metabolic age model trained in the Nightingale Metabolomics data quantified afte 2020
#We built two models for it: a linear model and an ElasticNET model. They both give similar results.

#Libraries
require(purrr)
library(readr)
################################
## Load the functions and data##
################################
# ! To modify by the user: Load metabolic features
metabo_mat <- as.data.frame(read_csv("~/Desktop/LUMC/Requantification/Data/LLS_requantified.csv"))
#Please be careful that you have rownames of the metabolomics dataset
rownames(metabo_mat)<-metabo_mat$sampleid
#Set the right path to where the scripts are:
setwd("~/Desktop/LUMC/metaboAge/MetaboAge2/apply_metaboAge2")
#Load the functions
source('functions_metaboAge2.R')
#Load the PARAMETERS to compute the metaboAge
PARAM_MetaboAge2<-readRDS("PARAM_MetaboAge2_2022_04_26.RData")
#Metabolic names translator
metabo_names_translator<-readRDS("metabolomic_feature_translator_BBMRI.rds")


######################################################################
## Correcting common differences in Nightingale Health raw datasets ##
######################################################################
# Correct the unit of crea
if((mean(metabo_mat$crea)/BBMRI_summaries["crea","Means"])>500){
  metabo_mat$crea<-metabo_mat$crea*0.001
}
#Calculate faw6_faw3 if missing
if(length(which(colnames(metabo_mat)=="faw6_faw3"))==0){
  metabo_mat$faw6_faw3<-metabo_mat$faw6/metabo_mat$faw3
}

######################################
## Find the metabolic features names##
######################################
#avoid case-sensitive alternative names
colnames(metabo_mat)<-tolower(colnames(metabo_mat))
#Looking for alternative names
nam<-find_BBMRI_names(colnames(metabo_mat))
i<-which(nam$BBMRI_names %in% metabo_names_translator$BBMRI_names)
metabo_mat<-metabo_mat[,i]
colnames(metabo_mat)<-nam$BBMRI_names[i]

####################################
## Apply metaboAge to the new data##
####################################
# Apply the quality control:
new_mat <- QCprep(as.matrix(metabo_mat[,PARAM_MetaboAge2$MET]),PARAM=PARAM_MetaboAge2)
#Apply the ElasticNET model (RECOMMENDED)
metaboAge_EN <- apply.fit.metaboAge(new_mat,PARAM=PARAM_MetaboAge2, model_type="EN")

#Apply the Linear model
metaboAge_LM <- apply.fit.metaboAge(new_mat,PARAM=PARAM_MetaboAge2, model_type="LM")

