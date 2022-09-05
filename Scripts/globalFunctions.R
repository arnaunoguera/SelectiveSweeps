#########################################
######## LIBRARIES ######################
#########################################

library(dplyr)
library(data.table)
library(tidyverse)
library(DBI)
library(dbplyr)
library(RMySQL)
library(stringr)  
library(ggplot2)
library(lawstat)



#########################################
######## VARIALBES ######################
#########################################

# For each population, which metapop they belong to, and associated color
popsData <- data.table(
	"pops" =c("YRI", "LWK", "GWD", "MSL", "ESN", "ACB", "ASW", 
		"CEU", "TSI", "FIN", "GBR", "IBS",
		"CHB", "JPT", "CHS", "CDX", "KHV", 
		"GIH", "PJL", "BEB", "STU", "ITU"),
	"meta"=c(rep("AFR",7),rep("EUR",5),rep("EAS",5),rep("SAS",5)),
	"colors" =c("#F5F39F","#F0F56B","#F7F14A","#F5F287","#F5EC05",
		"#FAFAC8","#FAFABE","#2D74B2","#ACD1F2","#6FA6D6","#5691C4",
		"#8AB9E3","#33B033","#6DD66D","#37B337","#008F00","#90DE90",
		"#A965BA","#D5A3E3","#8E3EA3","#B79BBD","#C587D6"))


# For each metapopulation, associated color
metapopsData <- data.table(
	"meta"=c('AFR', 'EAS' , 'EUR', 'SAS'),
	"colors"=c("#F7F14A", "#33B033", "#5691C4", "#A965BA"))

# PopsbyMetapopts
AFRpops<-c('ACB','ASW','ESN','GWD','LWK','MSL','YRI')
EASpops<-c('CDX','CHB','CHS','JPT','KHV')
EURpops<-c('CEU','FIN','GBR','IBS','TSI')
SASpops<-c('BEB','GIH','ITU','PJL','STU')

#Vector de colors i poblacions
popPal<- c("YRI"="#F5F39F", "LWK"="#F0F56B", "GWD"="#F7F14A", 
           "MSL"="#F5F287", "ESN"="#F5EC05", "ACB"="#FAFAC8", 
           "ASW"="#FAFABE", "CLM"="#AB3B35", "MXL"="#BA5852", 
           "PEL"="#CC7570", "PUR"="#DE9D99", "CEU"="#2D74B2", 
           "TSI"="#ACD1F2", "FIN"="#6FA6D6", "GBR"="#5691C4", 
           "IBS"="#8AB9E3", "CHB"="#33B033", "JPT"="#6DD66D", 
           "CHS"="#37B337", "CDX"="#008F00", "KHV"="#90DE90", 
           "GIH"="#A965BA", "PJL"="#D5A3E3", "BEB"="#8E3EA3", 
           "STU"="#B79BBD", "ITU"="#C587D6")


#########################################
####### FUNCTIONS #######################
#########################################

# ie <- 'chr2:109450405-109606617'
# For a chr, start & end, the info on one category statistic info
# dbNames[ihs, ihs_pval, nsl, nsl_pval, isafe, isafe_pval, functional, geva]

sqlToDf <- function(chr, coordStart, coordEnd, dbName){
        con <- dbConnect(RMySQL::MySQL(),
                user='shiny',
                password='shinyServer?21',
                dbname=dbName,
                host='localhost')
        myQuery <- paste0("SELECT * FROM ",
                "chr", chr,
                " WHERE physicalPos BETWEEN ", coordStart,
                " AND ", coordEnd)
        df <- dbGetQuery(con, myQuery)
        dbDisconnect(con)
        return(df %>% as.data.table)
}


# Get isafe + geva at once
getInfoRegion <- function(region){
  regionSep <- strsplit(region, ":")[[1]]
  CC <- as.numeric(regionSep[1])
  SS <- as.numeric(strsplit(regionSep[2], "-")[[1]][1])
  EE <- as.numeric(strsplit(regionSep[2], "-")[[1]][2])
  isafeData <- sqlToDf(CC, SS, EE, 'isafe')
  gevaData <- sqlToDf (CC, SS, EE, 'geva')
  #Faig inner_join perquè les variants han d'estar necessàriament a les dues taules
  return(dplyr::inner_join(isafeData, gevaData, by = 'rsid')) 
}
