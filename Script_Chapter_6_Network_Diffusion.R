# Chapter 6: Zebrafish Network Diffusion ####
# Author: R. Liscovsky (rliscovs@ed.ac.uk)

# 1. Install & Load Packages ####

if (!require("netdiffuseR")) install.packages("netdiffuseR")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("readr")) install.packages("readr")
if (!require("dummies")) install.packages("dummies")
if (!require("devtools")) install.packages("devtools")
if (!require("broom")) install.packages("broom")
if (!require("MASS")) install.packages("MASS")
if (!require("rpart")) install.packages("rpart")
if (!require("rpart.plot")) install.packages("rpart.plot")

library(netdiffuseR)
library(ggplot2)
library(dplyr)
library(readr)
library(dummies)
library(devtools)
library(broom)
library(MASS)
library(rpart)
library(rpart.plot)


# 2. Read Data In ####

edgelist <- read.csv("Data/Edgelist_INST_non_weighted(1996-2017).csv", 
                     header=TRUE, 
                     fileEncoding="UTF-8-BOM", 
                     stringsAsFactors=TRUE, 
                     sep = ",")

attributes <- read.csv("Data/Attribute_INST(1996-2017).csv", 
                       header=TRUE, 
                       fileEncoding="UTF-8-BOM", 
                       stringsAsFactors=TRUE, 
                       sep = ",")

## Inspect ####
glimpse(edgelist)
glimpse(attributes)

## Transformations ####
attributes$pub_no <- as.numeric(as.character(attributes$pub_no))
attributes$ego_net_size_no <- as.numeric(as.character(attributes$ego_net_size_no))
attributes$collab_country_no <- as.numeric(as.character(attributes$collab_country_no))
attributes$n_refs <- as.numeric(as.character(attributes$n_refs))
attributes$zebra_cited_no <- as.numeric(as.character(attributes$zebra_cited_no))
attributes$zebra_cited_per <- as.numeric(as.character(attributes$zebra_cited_per))


# 3. Build Diffnet object ####

diffnet <- edgelist_to_diffnet(
  edgelist = edgelist[,1:2],
  t0 = edgelist$time,
  t1 = edgelist$time,
  dat = attributes,
  idvar = "id",
  toavar = "toa",
  timevar = "time",
  undirected = getOption("diffnet.undirected", FALSE), # When TRUE only the lower triangle of the adjacency matrix will considered (faster
  keep.isolates = FALSE, # When FALSE, rows with NA/NULL values (isolated vertices unless have autolink) will be droped (see details).
  warn.coercion = FALSE # by default unserveyed individuals are excluded from the adjacency matrix by setting variables in netvars equal to NA when the nominated id can't be found in idvar.
  )

## Get summary statistics ####
summary(diffnet)

# Static and Dynamic Attributes ####

## Extract the graph data and create a diffnet object from scratch ####

graph <- diffnet$graph
ids <- diffnet$meta$ids
graph <- Map(function(g) {
  dimnames(g) <- list(ids,ids)
  g
}, g=graph)
attrs <- diffnet.attrs(diffnet, "vertex", "static", as.df=TRUE) # they have the format of dynamic but they are static (losts of NA)
toa <- diffnet.toa(diffnet)

## Extract dynamic attributes ####

attrs_dyn <- diffnet.attrs(diffnet, "vertex", "dyn", as.df=TRUE) # The real dynamic attributes

## And export this as a table into SQL to get the countries of each afid ####
write.csv(attrs, file = "attrs_new.csv") 

## Create the static attributes ####
# Export the new table from SQL open it in Excel, delete the first column 
# and "per" (keep only: toa, id, country) and remove duplicates

attrs.static <- read.csv("Data/attrs._static_with_mobility&ref.csv", 
                         header=TRUE, 
                         fileEncoding="UTF-8-BOM", 
                         stringsAsFactors=TRUE, 
                         sep = ",")
head(attrs.static)

## Some additional transformations... ####
attrs.static$inst_size <- as.numeric(as.character(attrs.static$inst_size))
attrs.static$inst_mob_size <- as.numeric(as.character(attrs.static$inst_mob_size))
attrs.static$mob_share <- as.numeric(as.character(attrs.static$mob_share))

# 4. Create a New Diffnet Object with the attributes ####

v <- new_diffnet(graph, 
                           toa=toa, 
                           vertex.static.attrs=attrs.static, 
                           vertex.dyn.attrs = attrs_dyn, 
                           id.and.per.vars = c("id", "per"))
# Summary statistics:
summary(diffnet_new)

# 5. Calculate Network Exposure ####

## Simple exposure: raw count of neighbors; weights are not considered ####

e1.count <- exposure(diffnet_new, valued=FALSE, normalize=FALSE)

# Write as .csv
write.csv(e1.count, file = "Data/whole_exp_count_GLOBAL.csv")

## Normalized exposure (%) ####
e1.perc <- exposure(diffnet_new, valued=FALSE, normalize=TRUE)

# Write as .csv 
write.csv(e1.perc, file = "Data/whole_exp_per_GLOBAL.csv")

## Lagged exposure (lag = 1) ####

lagged_expo <- exposure(diffnet_new, valued=FALSE, normalize=FALSE, lags = 1L)

# Write as .csv 
write.csv(lagged_expo, file = "Data/whole(lagged)_exp_count_GLOBAL.csv")

## Lagged exposure normalized ####
lagged_expo_per <- exposure(diffnet_new, valued=FALSE, normalize=TRUE, lags = 1L)

write.csv(lagged_expo_per, file = "Data/whole(lagged)_exp_per_GLOBAL.csv")

# Netdiffuser automatically identifies whether the input is dynamic or not.
diffnet_new[["lagged_expo"]] <- lagged_expo
diffnet_new[["adopted"]]  <- toa_mat(diffnet_new)$cumadopt

##  Getting the network diffusion data as a frame ####
mydata <- as.data.frame(diffnet_new)
head(mydata)

# 6. Adopters classification  ####
ALL_adopters <- cbind(as.data.frame(classify(diffnet_new)), diffnet_new$toa)
write.csv(ALL_adopters, file = "Data/Adopters/Adopters_ALL.csv")

# 7. Sub-graph Partitions #### 
countries <- unique(diffnet_new[["country"]])
View(countries)

## North America #### 
# USA & Canada (2 countries)

diffnet_NorthAmerica <- diffnet_new[diffnet_new[["country"]] %in% countries[c(6,5)]]
dim(diffnet_NorthAmerica)

## LAZEN #### 
# (Arg, Bra, Chl, Col, Ecu, Mex, Per, Ury) (8 countries) 

diffnet_LAZEN <- diffnet_new[diffnet_new[["country"]] %in% countries[c(46,25,27,43,105,83,70,47)]]
dim(diffnet_LAZEN)

## EU-COST (core) ####   
# (Deu, Gbr, Fra, Nld) (4 countries)

diffnet_EU_core_COST <- diffnet_new[diffnet_new[["country"]] %in% countries[c(4,15,8,22)]]
dim(diffnet_EU_core_COST)

## EU-COST (periphery) ####
# (alb, aut, bel, bgr, hrv, cze, dnk, est, fin, grc, hun, isl, irl, ita, lva, lux, mkd, nor, pol, prt,
#rou, srb, svk, svn, esp, swe, che, tur, isr) With >2 pubs
# (29 countries)

diffnet_EU_peri_COST <- diffnet_new[diffnet_new[["country"]] %in% countries[c(108, 34, 14, 85, 79, 61, 38, 48, 40, 32, 
                                                                              26, 77, 7, 3, 97, 82, 109, 17, 39, 30, 75, 
                                                                              50, 73, 41, 16, 21, 18, 66, 33)]]
dim(diffnet_EU_peri_COST)

## Asia-Pacific (core) ####
# (chn, jpn, aus, twn, kor) > 60% total output (5 countries)

diffnet_ASIA_OC_core <- diffnet_new[diffnet_new[["country"]] %in% countries[c(9, 2, 10, 28, 20)]]

dim(diffnet_ASIA_OC_core)

## Asia-Pacific (periphery) #### 
# (sgp, ind, khg, nzl, sau, mys, mac, irn, tha, pak, are, phl, bgd, qat, vnm, idn, jor, 
# lbn, kwt, lka, omn, npl, fji, mng, brn) (25 countries)

diffnet_ASIA_OC_peri <- diffnet_new[diffnet_new[["country"]] %in% countries[c(12, 11, 31, 36, 24, 42, 65, 64, 55, 1,
                                                                              54, 74, 81, 114, 19, 62, 93, 95, 89, 44,
                                                                              101, 71, 104, 110, 107)]]

dim(diffnet_ASIA_OC_peri)


## All six communities combined #### 

diffnet_ALL <- diffnet_new[diffnet_new[["country"]] %in% countries[c(6,5,
                                                                     46,25,27,43,105,83,70,47,
                                                                     4,15,8,22, 
                                                                     108, 34, 14, 85, 79, 61, 38, 48, 40, 32, 
                                                                     26, 77, 7, 3, 97, 82, 109, 17, 39, 30, 75, 
                                                                     50, 73, 41, 16, 21, 18, 66, 33,
                                                                     9, 2, 10, 28, 20,
                                                                     12, 11, 31, 36, 24, 42, 65, 64, 55, 1,
                                                                     54, 74, 81, 114, 19, 62, 93, 95, 89, 44,
                                                                     101, 71, 104, 110, 107)]]
dim(diffnet_ALL)

# 8. Adoptions per Community #### 

cols <- c("green","red", "blue4", "blue", "yellow", "yellow3")
adoption_comp <- {
  plot_adopters(diffnet_NorthAmerica, bg = cols[1], include.legend = FALSE, what="cumadopt", main = "Cumulative Adopters Compared")
  plot_adopters(diffnet_LAZEN, bg = cols[2], add=TRUE, what="cumadopt")
  plot_adopters(diffnet_EU_core_COST, bg = cols[3], add=TRUE, what="cumadopt")
  plot_adopters(diffnet_EU_peri_COST, bg = cols[4], add=TRUE, what="cumadopt")
  plot_adopters(diffnet_ASIA_OC_core, bg = cols[5], add=TRUE, what="cumadopt")
  plot_adopters(diffnet_ASIA_OC_peri, bg = cols[6], add=TRUE, what="cumadopt")
}

# Add a legend
legend("topleft", bty="n", legend = c("North America","LAZEN", "EU-COST (core)", "EU-COST (periphery)", 
                                      "Asia-Pacific (core)", "Asia-Pacific (periphery)"), fill=cols)

setwd("C:/Users/s1529697/Documents/R/Network exposure (NEW)/Adopters")

# Write as .csv
adopters_6_commu

# 9. Exposures per Community (raw) ####

# North America:
exp_reg_count_NorthAmer <- exposure(diffnet_NorthAmerica, valued=FALSE, normalize=FALSE)
write.csv(exp_reg_count_NorthAmer, file = "Data/Regional exposure (unclean)/Normal/reg_exp_count_North America.csv")

lagged_exp_reg_count_NorthAmer <- exposure(diffnet_NorthAmerica, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_reg_count_NorthAmer, file = "Data/Regional exposure (unclean)/Lagged/lagged_reg_exp_count_North America.csv")

## LAZEN ####
exp_reg_count_LAZEN <- exposure(diffnet_LAZEN, valued=FALSE, normalize=FALSE)
write.csv(exp_reg_count_LAZEN, file = "Data/Regional exposure (unclean)/Normal/reg_exp_count_LAZEN.csv")

lagged_exp_reg_count_LAZEN <- exposure(diffnet_LAZEN, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_reg_count_LAZEN, file = "Data/Regional exposure (unclean)/Lagged/lagged_reg_exp_count_LAZEN.csv")

## EU-COST (core) ####
exp_reg_count_EU_core <- exposure(diffnet_EU_core_COST, valued=FALSE, normalize=FALSE)
write.csv(exp_reg_count_EU_core, file = "Data/Regional exposure (unclean)/Normal/reg_exp_count_EU (core).csv")

lagged_exp_reg_count_EU_core <- exposure(diffnet_EU_core_COST, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_reg_count_EU_core, file = "Data/Regional exposure (unclean)/Lagged/lagged_reg_exp_count_EU (core).csv")

## EU-COST (periphery) ####
exp_reg_count_EU_peri <- exposure(diffnet_EU_peri_COST, valued=FALSE, normalize=FALSE)
write.csv(exp_reg_count_EU_peri, file = "Data/Regional exposure (unclean)/Normal/reg_exp_count_EU (periphery).csv")

lagged_exp_reg_count_EU_peri <- exposure(diffnet_EU_peri_COST, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_reg_count_EU_peri, file = "Data/Regional exposure (unclean)/Lagged/lagged_reg_exp_count_EU (periphery).csv")

## Asia-Pacific (core) ####
exp_reg_count_ASOC_core <- exposure(diffnet_ASIA_OC_core, valued=FALSE, normalize=FALSE)
write.csv(exp_reg_count_ASOC_core, file = "Data/Regional exposure (unclean)/Normal/reg_exp_count_AsiaPacific (core).csv")

lagged_exp_reg_count_ASOC_core <- exposure(diffnet_ASIA_OC_core, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_reg_count_ASOC_core, file = "Data/Regional exposure (unclean)/Lagged/lagged_reg_exp_count_AsiaPacific (core).csv")

## Asia-Pacific (periphery) ####
exp_reg_count_ASOC_peri <- exposure(diffnet_ASIA_OC_peri, valued=FALSE, normalize=FALSE)
write.csv(exp_reg_count_ASOC_peri, file = "Data/Regional exposure (unclean)/Normal/reg_exp_count_AsiaPacific (periphery).csv")

lagged_exp_reg_count_ASOC_peri <- exposure(diffnet_ASIA_OC_peri, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_reg_count_ASOC_peri, file = "Data/Regional exposure (unclean)/Lagged/lagged_reg_exp_count_AsiaPacific (periphery).csv")

## EU-COST (all) ####
diffnet_EU_COST_all <- diffnet_new[diffnet_new[["country"]] %in% countries[c(4,15,8,22, 108, 34, 14, 85, 79, 61, 38, 48, 40, 32, 
                                                                             26, 77, 7, 3, 97, 82, 109, 17, 39, 30, 75, 
                                                                             50, 73, 41, 16, 21, 18, 66, 33)]]

exp_reg_count_EU_COST_all <- exposure(diffnet_EU_COST_all, valued=FALSE, normalize=FALSE)
write.csv(exp_reg_count_EU_COST_all, file = "Data/Regional exposure (unclean)/Normal/reg_exp_count_EU-COST(all).csv")

lagged_exp_reg_count_EU_COST_all <- exposure(diffnet_EU_COST_all, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_reg_count_EU_COST_all, file = "Data/Regional exposure (unclean)/Lagged/lagged_reg_exp_count_EU-COST(all).csv")


## Asia-Pacific (all) ####

diffnet_ASIA_OC_all <- diffnet_new[diffnet_new[["country"]] %in% countries[c(9, 2, 10, 28, 20, 12, 11, 31, 36, 24, 42, 65, 64, 55, 1,
                                                                             54, 74, 81, 114, 19, 62, 93, 95, 89, 44,
                                                                             101, 71, 104, 110, 107)]]

exp_reg_count_ASOC_all <- exposure(diffnet_ASIA_OC_all, valued=FALSE, normalize=FALSE)
write.csv(exp_reg_count_ASOC_all, file = "Data/Regional exposure (unclean)/Normal/reg_exp_count_Asia-Pacific(all).csv")

lagged_exp_reg_count_ASOC_all <- exposure(diffnet_ASIA_OC_all, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_reg_count_ASOC_all, file = "Data/Regional exposure (unclean)/Lagged/lagged_reg_exp_count_Asia-Pacific(all).csv")


# 10. Exposures per Country (raw) ####

## North America ####
diffnet_USA <- diffnet_new[diffnet_new[["country"]] %in% countries[6]]
exp_nat_count_USA <- exposure(diffnet_USA, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_USA, file = "Data/National exposure/North America/Normal/nat_exp_count_USA.csv")

lagged_exp_nat_count_USA <- exposure(diffnet_USA, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_USA, file = "Data/National exposure/North America/Lagged/lagged_nat_exp_count_USA.csv")

diffnet_CAN <- diffnet_new[diffnet_new[["country"]] %in% countries[5]]
exp_nat_count_CAN <- exposure(diffnet_CAN, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_CAN, file = "Data/National exposure/North America/Normal/nat_exp_count_CAN.csv")

lagged_exp_nat_count_CAN <- exposure(diffnet_CAN, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_CAN, file = "Data/National exposure/North America/Lagged/lagged_nat_exp_count_CAN.csv")

## LAZEN  ####

diffnet_ARG <- diffnet_new[diffnet_new[["country"]] %in% countries[46]]
exp_nat_count_ARG <- exposure(diffnet_ARG, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_ARG, file = "Data/National exposure/LAZEN/Normal/nat_exp_count_ARG.csv")

lagged_exp_nat_count_ARG <- exposure(diffnet_ARG, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_ARG, file = "Data/National exposure/LAZEN/Lagged/nlagged_nat_exp_count_ARG.csv")

diffnet_BRA <- diffnet_new[diffnet_new[["country"]] %in% countries[25]]
exp_nat_count_BRA <- exposure(diffnet_BRA, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_BRA, file = "Data/National exposure/LAZEN/Normal/nnat_exp_count_BRA.csv")

lagged_exp_nat_count_BRA <- exposure(diffnet_BRA, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_BRA, file = "Data/National exposure/LAZEN/Lagged/lagged_nat_exp_count_BRA.csv")

diffnet_CHL <- diffnet_new[diffnet_new[["country"]] %in% countries[27]]
exp_nat_count_CHL <- exposure(diffnet_CHL, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_CHL, file = "Data/National exposure/LAZEN/Normal/nnat_exp_count_CHL.csv")

lagged_exp_nat_count_CHL <- exposure(diffnet_CHL, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_CHL, file = "Data/National exposure/LAZEN/Lagged/lagged_nat_exp_count_CHL.csv")

diffnet_COL <- diffnet_new[diffnet_new[["country"]] %in% countries[43]]
exp_nat_count_COL <- exposure(diffnet_COL, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_COL, file = "Data/National exposure/LAZEN/Normal/nnat_exp_count_COL.csv")

lagged_exp_nat_count_COL <- exposure(diffnet_COL, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_COL, file = "Data/National exposure/LAZEN/Lagged/lagged_nat_exp_count_COL.csv")

diffnet_ECU <- diffnet_new[diffnet_new[["country"]] %in% countries[105]]
exp_nat_count_ECU <- exposure(diffnet_ECU, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_ECU, file = "Data/National exposure/LAZEN/Normal/nnat_exp_count_ECU.csv")

lagged_exp_nat_count_ECU <- exposure(diffnet_ECU, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_ECU, file = "Data/National exposure/LAZEN/Lagged/lagged_nat_exp_count_ECU.csv")

diffnet_MEX <- diffnet_new[diffnet_new[["country"]] %in% countries[83]]
exp_nat_count_MEX <- exposure(diffnet_MEX, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_MEX, file = "Data/National exposure/LAZEN/Normal/nnat_exp_count_MEX.csv")

lagged_exp_nat_count_MEX <- exposure(diffnet_MEX, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_MEX, file = "Data/National exposure/LAZEN/Lagged/lagged_nat_exp_count_MEX.csv")

diffnet_PER <- diffnet_new[diffnet_new[["country"]] %in% countries[70]]
exp_nat_count_PER <- exposure(diffnet_PER, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_PER, file = "Data/National exposure/LAZEN/Normal/nnat_exp_count_PER.csv")

lagged_exp_nat_count_PER <- exposure(diffnet_PER, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_PER, file = "Data/National exposure/LAZEN/Lagged/lagged_nat_exp_count_PER.csv")

diffnet_URY <- diffnet_new[diffnet_new[["country"]] %in% countries[47]]
exp_nat_count_URY <- exposure(diffnet_URY, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_URY, file = "Data/National exposure/LAZEN/Normal/nnat_exp_count_URY.csv")

lagged_exp_nat_count_URY <- exposure(diffnet_URY, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_URY, file = "Data/National exposure/LAZEN/Lagged/lagged_nat_exp_count_URY.csv")

## EU-COST (core)   ####

diffnet_DEU <- diffnet_new[diffnet_new[["country"]] %in% countries[15]]
exp_nat_count_DEU <- exposure(diffnet_DEU, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_DEU, file = "Data/National exposure/EU-COST (core)/Normal/nat_exp_count_DEU.csv")

lagged_exp_nat_count_DEU <- exposure(diffnet_DEU, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_DEU, file = "Data/National exposure/EU-COST (core)/Lagged/lagged_nat_exp_count_DEU.csv")

diffnet_GBR <- diffnet_new[diffnet_new[["country"]] %in% countries[4]]
exp_nat_count_GBR <- exposure(diffnet_GBR, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_GBR, file = "Data/National exposure/EU-COST (core)/Normal/nat_exp_count_GBR.csv")

lagged_exp_nat_count_GBR <- exposure(diffnet_GBR, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_GBR, file = "Data/National exposure/EU-COST (core)/Lagged/lagged_nat_exp_count_GBR.csv")

diffnet_FRA <- diffnet_new[diffnet_new[["country"]] %in% countries[8]]
exp_nat_count_FRA <- exposure(diffnet_FRA, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_FRA, file = "Data/National exposure/EU-COST (core)/Normal/nat_exp_count_FRA.csv")

lagged_exp_nat_count_FRA <- exposure(diffnet_FRA, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_FRA, file = "Data/National exposure/EU-COST (core)/Lagged/lagged_nat_exp_count_FRA.csv")

diffnet_NLD <- diffnet_new[diffnet_new[["country"]] %in% countries[22]]
exp_nat_count_NLD <- exposure(diffnet_NLD, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_NLD, file = "Data/National exposure/EU-COST (core)/Normal/nat_exp_count_NLD.csv")

lagged_exp_nat_count_NLD <- exposure(diffnet_NLD, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_NLD, file = "Data/National exposure/EU-COST (core)/Lagged/lagged_nat_exp_count_NLD.csv")

# 9.4 EU-COST (periphery):

diffnet_ALB <- diffnet_new[diffnet_new[["country"]] %in% countries[108]]
exp_nat_count_ALB <- exposure(diffnet_ALB, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_ALB, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_ALB.csv")

lagged_exp_nat_count_ALB <- exposure(diffnet_ALB, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_ALB, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_ALB.csv")

diffnet_AUT <- diffnet_new[diffnet_new[["country"]] %in% countries[34]]
exp_nat_count_AUT <- exposure(diffnet_AUT, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_AUT, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_AUT.csv")

lagged_exp_nat_count_AUT <- exposure(diffnet_AUT, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_AUT, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_AUT.csv")

diffnet_BEL <- diffnet_new[diffnet_new[["country"]] %in% countries[14]]
exp_nat_count_BEL <- exposure(diffnet_BEL, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_BEL, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_BEL.csv")

lagged_exp_nat_count_BEL <- exposure(diffnet_BEL, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_BEL, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_BEL.csv")

diffnet_BGR <- diffnet_new[diffnet_new[["country"]] %in% countries[85]]
exp_nat_count_BGR <- exposure(diffnet_BGR, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_BGR, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_BGR.csv")

lagged_exp_nat_count_BGR <- exposure(diffnet_BGR, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_BGR, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_BGR.csv")

diffnet_HRV <- diffnet_new[diffnet_new[["country"]] %in% countries[79]]
exp_nat_count_HRV <- exposure(diffnet_HRV, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_HRV, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_HRV.csv")

lagged_exp_nat_count_HRV <- exposure(diffnet_HRV, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_HRV, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_HRV.csv")

diffnet_CZE <- diffnet_new[diffnet_new[["country"]] %in% countries[61]]
exp_nat_count_CZE <- exposure(diffnet_CZE, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_CZE, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_CZE.csv")

lagged_exp_nat_count_CZE <- exposure(diffnet_CZE, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_CZE, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_CZE.csv")

diffnet_DNK <- diffnet_new[diffnet_new[["country"]] %in% countries[38]]
exp_nat_count_DNK <- exposure(diffnet_DNK, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_DNK, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_DNK.csv")

lagged_exp_nat_count_DNK <- exposure(diffnet_DNK, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_DNK, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_DNK.csv")

diffnet_EST <- diffnet_new[diffnet_new[["country"]] %in% countries[48]]
exp_nat_count_EST <- exposure(diffnet_EST, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_EST, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_EST.csv")

lagged_exp_nat_count_EST <- exposure(diffnet_EST, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_EST, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_EST.csv")

diffnet_FIN <- diffnet_new[diffnet_new[["country"]] %in% countries[40]]
exp_nat_count_FIN <- exposure(diffnet_FIN, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_FIN, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_FIN.csv")

lagged_exp_nat_count_FIN <- exposure(diffnet_FIN, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_FIN, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_FIN.csv")

diffnet_GRC <- diffnet_new[diffnet_new[["country"]] %in% countries[32]]
exp_nat_count_GRC <- exposure(diffnet_GRC, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_GRC, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_GRC.csv")

lagged_exp_nat_count_GRC <- exposure(diffnet_GRC, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_GRC, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_GRC.csv")

diffnet_HUN <- diffnet_new[diffnet_new[["country"]] %in% countries[26]]
exp_nat_count_HUN <- exposure(diffnet_HUN, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_HUN, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_HUN.csv")

lagged_exp_nat_count_HUN <- exposure(diffnet_HUN, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_HUN, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_HUN.csv")

diffnet_ISL <- diffnet_new[diffnet_new[["country"]] %in% countries[77]]
exp_nat_count_ISL <- exposure(diffnet_ISL, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_ISL, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_ISL.csv")

lagged_exp_nat_count_ISL <- exposure(diffnet_ISL, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_ISL, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_ISL.csv")

diffnet_IRL <- diffnet_new[diffnet_new[["country"]] %in% countries[7]]
exp_nat_count_IRL <- exposure(diffnet_IRL, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_IRL, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_IRL.csv")

lagged_exp_nat_count_IRL <- exposure(diffnet_IRL, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_IRL, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_IRL.csv")

diffnet_ITA <- diffnet_new[diffnet_new[["country"]] %in% countries[3]]
exp_nat_count_ITA <- exposure(diffnet_ITA, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_ITA, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_ITA.csv")

lagged_exp_nat_count_ITA <- exposure(diffnet_ITA, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_ITA, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_ITA.csv")

diffnet_LVA <- diffnet_new[diffnet_new[["country"]] %in% countries[97]]
exp_nat_count_LVA <- exposure(diffnet_LVA, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_LVA, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_LVA.csv")

lagged_exp_nat_count_LVA <- exposure(diffnet_LVA, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_LVA, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_LVA.csv")

diffnet_LUX <- diffnet_new[diffnet_new[["country"]] %in% countries[82]]
exp_nat_count_LUX <- exposure(diffnet_LUX, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_LUX, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_LUX.csv")

lagged_exp_nat_count_LUX <- exposure(diffnet_LUX, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_LUX, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_LUX.csv")

diffnet_MKD <- diffnet_new[diffnet_new[["country"]] %in% countries[109]]
exp_nat_count_MKD <- exposure(diffnet_MKD, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_MKD, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_MKD.csv")

lagged_exp_nat_count_MKD <- exposure(diffnet_MKD, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_MKD, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_MKD.csv")

diffnet_NOR <- diffnet_new[diffnet_new[["country"]] %in% countries[17]]
exp_nat_count_NOR <- exposure(diffnet_NOR, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_NOR, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_NOR.csv")

lagged_exp_nat_count_NOR <- exposure(diffnet_NOR, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_NOR, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_NOR.csv")

diffnet_POL <- diffnet_new[diffnet_new[["country"]] %in% countries[39]]
exp_nat_count_POL <- exposure(diffnet_POL, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_POL, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_POL.csv")

lagged_exp_nat_count_POL <- exposure(diffnet_POL, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_POL, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_POL.csv")

diffnet_PRT <- diffnet_new[diffnet_new[["country"]] %in% countries[30]]
exp_nat_count_PRT <- exposure(diffnet_PRT, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_PRT, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_PRT.csv")

lagged_exp_nat_count_PRT <- exposure(diffnet_PRT, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_PRT, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_PRT.csv")

diffnet_ROU <- diffnet_new[diffnet_new[["country"]] %in% countries[75]]
exp_nat_count_ROU <- exposure(diffnet_ROU, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_ROU, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_ROU.csv")

lagged_exp_nat_count_ROU <- exposure(diffnet_ROU, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_ROU, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_ROU.csv")

diffnet_SRB <- diffnet_new[diffnet_new[["country"]] %in% countries[50]]
exp_nat_count_SRB <- exposure(diffnet_SRB, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_SRB, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_SRB.csv")

lagged_exp_nat_count_SRB <- exposure(diffnet_SRB, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_SRB, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_SRB.csv")

diffnet_SVK <- diffnet_new[diffnet_new[["country"]] %in% countries[73]]
exp_nat_count_SVK <- exposure(diffnet_SVK, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_SVK, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_SVK.csv")

lagged_exp_nat_count_SVK <- exposure(diffnet_SVK, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_SVK, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_SVK.csv")

diffnet_SVN <- diffnet_new[diffnet_new[["country"]] %in% countries[41]]
exp_nat_count_SVN <- exposure(diffnet_SVN, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_SVN, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_SVN.csv")

lagged_exp_nat_count_SVN <- exposure(diffnet_SVN, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_SVN, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_SVN.csv")

diffnet_ESP <- diffnet_new[diffnet_new[["country"]] %in% countries[16]]
exp_nat_count_ESP <- exposure(diffnet_ESP, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_ESP, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_ESP.csv")

lagged_exp_nat_count_ESP <- exposure(diffnet_ESP, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_ESP, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_ESP.csv")

diffnet_SWE <- diffnet_new[diffnet_new[["country"]] %in% countries[21]]
exp_nat_count_SWE <- exposure(diffnet_SWE, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_SWE, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_SWE.csv")

lagged_exp_nat_count_SWE <- exposure(diffnet_SWE, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_SWE, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_SWE.csv")

diffnet_CHE <- diffnet_new[diffnet_new[["country"]] %in% countries[18]]
exp_nat_count_CHE <- exposure(diffnet_CHE, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_CHE, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_CHE.csv")

lagged_exp_nat_count_CHE <- exposure(diffnet_CHE, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_CHE, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_CHE.csv")

diffnet_TUR <- diffnet_new[diffnet_new[["country"]] %in% countries[66]]
exp_nat_count_TUR <- exposure(diffnet_TUR, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_TUR, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_TUR.csv")

lagged_exp_nat_count_TUR <- exposure(diffnet_TUR, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_TUR, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_TUR.csv")

diffnet_ISR <- diffnet_new[diffnet_new[["country"]] %in% countries[33]]
exp_nat_count_ISR <- exposure(diffnet_ISR, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_ISR, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_ISR.csv")

lagged_exp_nat_count_ISR <- exposure(diffnet_ISR, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_ISR, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_ISR.csv")


# 9.5 Asia-Pacific (core):

diffnet_CHN <- diffnet_new[diffnet_new[["country"]] %in% countries[9]]
exp_nat_count_CHN <- exposure(diffnet_CHN, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_CHN, file = "Data\National exposure/Asia-Pacific (core)/Normal/nat_exp_count_CHN.csv")

lagged_exp_nat_count_CHN <- exposure(diffnet_CHN, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_CHN, file = "Data\National exposure/Asia-Pacific (core)/Lagged/lagged_nat_exp_count_CHN.csv")

diffnet_JPN <- diffnet_new[diffnet_new[["country"]] %in% countries[2]]
exp_nat_count_JPN <- exposure(diffnet_JPN, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_JPN, file = "Data\National exposure/Asia-Pacific (core)/Normal/nat_exp_count_JPN.csv")

lagged_exp_nat_count_JPN <- exposure(diffnet_JPN, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_JPN, file = "Data\National exposure/Asia-Pacific (core)/Lagged/lagged_nat_exp_count_JPN.csv")

diffnet_AUS <- diffnet_new[diffnet_new[["country"]] %in% countries[10]]
exp_nat_count_AUS <- exposure(diffnet_AUS, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_AUS, file = "Data\National exposure/Asia-Pacific (core)/Normal/nat_exp_count_AUS.csv")

lagged_exp_nat_count_AUS <- exposure(diffnet_AUS, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_AUS, file = "Data\National exposure/Asia-Pacific (core)/Lagged/lagged_nat_exp_count_AUS.csv")

diffnet_TWN <- diffnet_new[diffnet_new[["country"]] %in% countries[28]]
exp_nat_count_TWN <- exposure(diffnet_TWN, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_TWN, file = "Data\National exposure/Asia-Pacific (core)/Normal/nat_exp_count_TWN.csv")

lagged_exp_nat_count_TWN <- exposure(diffnet_TWN, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_TWN, file = "Data\National exposure/Asia-Pacific (core)/Lagged/lagged_nat_exp_count_TWN.csv")

diffnet_KOR <- diffnet_new[diffnet_new[["country"]] %in% countries[20]]
exp_nat_count_KOR <- exposure(diffnet_KOR, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_KOR, file = "Data\National exposure/Asia-Pacific (core)/Normal/nat_exp_count_KOR.csv")

lagged_exp_nat_count_KOR <- exposure(diffnet_KOR, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_KOR, file = "Data\National exposure/Asia-Pacific (core)/Lagged/lagged_nat_exp_count_KOR.csv")


# 9.6 Asia-Pacific (periphery):

diffnet_SGP <- diffnet_new[diffnet_new[["country"]] %in% countries[12]]
exp_nat_count_SGP <- exposure(diffnet_SGP, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_SGP, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_SGP.csv")

lagged_exp_nat_count_SGP <- exposure(diffnet_SGP, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_SGP, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_SGP.csv")

diffnet_IND <- diffnet_new[diffnet_new[["country"]] %in% countries[11]]
exp_nat_count_IND <- exposure(diffnet_IND, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_IND, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_IND.csv")

lagged_exp_nat_count_IND <- exposure(diffnet_IND, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_IND, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_IND.csv")

diffnet_HKG <- diffnet_new[diffnet_new[["country"]] %in% countries[31]]
exp_nat_count_HKG <- exposure(diffnet_HKG, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_HKG, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_HKG.csv")

lagged_exp_nat_count_HKG <- exposure(diffnet_HKG, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_HKG, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_HKG.csv")

diffnet_NZL <- diffnet_new[diffnet_new[["country"]] %in% countries[36]]
exp_nat_count_NZL <- exposure(diffnet_NZL, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_NZL, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_NZL.csv")

lagged_exp_nat_count_NZL <- exposure(diffnet_NZL, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_NZL, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_NZL.csv")

diffnet_SAU <- diffnet_new[diffnet_new[["country"]] %in% countries[24]]
exp_nat_count_SAU <- exposure(diffnet_SAU, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_SAU, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_SAU.csv")

lagged_exp_nat_count_SAU <- exposure(diffnet_SAU, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_SAU, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_SAU.csv")

diffnet_MYS <- diffnet_new[diffnet_new[["country"]] %in% countries[42]]
exp_nat_count_MYS <- exposure(diffnet_MYS, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_MYS, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_MYS.csv")

lagged_exp_nat_count_MYS <- exposure(diffnet_MYS, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_MYS, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_MYS.csv")

diffnet_MAC <- diffnet_new[diffnet_new[["country"]] %in% countries[65]]
exp_nat_count_MAC <- exposure(diffnet_MAC, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_MAC, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_MAC.csv")

lagged_exp_nat_count_MAC <- exposure(diffnet_MAC, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_MAC, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_MAC.csv")

diffnet_IRN <- diffnet_new[diffnet_new[["country"]] %in% countries[64]]
exp_nat_count_IRN <- exposure(diffnet_IRN, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_IRN, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_IRN.csv")

lagged_exp_nat_count_IRN <- exposure(diffnet_IRN, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_IRN, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_IRN.csv")

diffnet_THA <- diffnet_new[diffnet_new[["country"]] %in% countries[55]]
exp_nat_count_THA <- exposure(diffnet_THA, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_THA, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_THA.csv")

lagged_exp_nat_count_THA <- exposure(diffnet_THA, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_THA, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_THA.csv")

diffnet_PAK <- diffnet_new[diffnet_new[["country"]] %in% countries[1]]
exp_nat_count_PAK <- exposure(diffnet_PAK, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_PAK, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_PAK.csv")

lagged_exp_nat_count_PAK <- exposure(diffnet_PAK, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_PAK, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_PAK.csv")

diffnet_ARE <- diffnet_new[diffnet_new[["country"]] %in% countries[54]]
exp_nat_count_ARE <- exposure(diffnet_ARE, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_ARE, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_ARE.csv")

lagged_exp_nat_count_ARE <- exposure(diffnet_ARE, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_ARE, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_ARE.csv")

diffnet_PHL <- diffnet_new[diffnet_new[["country"]] %in% countries[74]]
exp_nat_count_PHL <- exposure(diffnet_PHL, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_PHL, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_PHL.csv")

lagged_exp_nat_count_PHL <- exposure(diffnet_PHL, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_PHL, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_PHL.csv") 

diffnet_BGD <- diffnet_new[diffnet_new[["country"]] %in% countries[81]]
exp_nat_count_BGD <- exposure(diffnet_BGD, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_BGD, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_BGD.csv")

lagged_exp_nat_count_BGD <- exposure(diffnet_BGD, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_BGD, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_BGD.csv")

diffnet_QAT <- diffnet_new[diffnet_new[["country"]] %in% countries[114]]
exp_nat_count_QAT <- exposure(diffnet_QAT, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_QAT, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_QAT.csv")

lagged_exp_nat_count_QAT <- exposure(diffnet_QAT, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_QAT, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_QAT.csv")

diffnet_VNM <- diffnet_new[diffnet_new[["country"]] %in% countries[19]]
exp_nat_count_VNM <- exposure(diffnet_VNM, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_VNM, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_VNM.csv")

lagged_exp_nat_count_VNM <- exposure(diffnet_VNM, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_VNM, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_VNM.csv")

diffnet_IDN <- diffnet_new[diffnet_new[["country"]] %in% countries[62]]
exp_nat_count_IDN <- exposure(diffnet_IDN, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_IDN, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_IDN.csv")

lagged_exp_nat_count_IDN <- exposure(diffnet_IDN, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_IDN, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_IDN.csv")

diffnet_JOR <- diffnet_new[diffnet_new[["country"]] %in% countries[93]]
exp_nat_count_JOR <- exposure(diffnet_JOR, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_JOR, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_JOR.csv")

lagged_exp_nat_count_JOR <- exposure(diffnet_JOR, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_JOR, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_JOR.csv")

diffnet_LBN <- diffnet_new[diffnet_new[["country"]] %in% countries[95]]
exp_nat_count_LBN <- exposure(diffnet_LBN, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_LBN, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_LBN.csv")

lagged_exp_nat_count_LBN <- exposure(diffnet_LBN, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_LBN, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_LBN.csv")

diffnet_KWT <- diffnet_new[diffnet_new[["country"]] %in% countries[89]]
exp_nat_count_KWT <- exposure(diffnet_KWT, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_KWT, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_KWT.csv")

lagged_exp_nat_count_KWT <- exposure(diffnet_KWT, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_KWT, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_KWT.csv")

diffnet_LKA <- diffnet_new[diffnet_new[["country"]] %in% countries[44]]
exp_nat_count_LKA <- exposure(diffnet_LKA, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_LKA, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_LKA.csv")

lagged_exp_nat_count_LKA <- exposure(diffnet_LKA, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_LKA, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_LKA.csv")

diffnet_OMN <- diffnet_new[diffnet_new[["country"]] %in% countries[101]]
exp_nat_count_OMN <- exposure(diffnet_OMN, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_OMN, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_OMN.csv")

lagged_exp_nat_count_OMN <- exposure(diffnet_OMN, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_OMN, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_OMN.csv")

diffnet_NPL <- diffnet_new[diffnet_new[["country"]] %in% countries[71]]
exp_nat_count_NPL <- exposure(diffnet_NPL, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_NPL, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_NPL.csv")

lagged_exp_nat_count_NPL<- exposure(diffnet_NPL, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_NPL, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_NPL.csv")

diffnet_FJI <- diffnet_new[diffnet_new[["country"]] %in% countries[104]]
exp_nat_count_FJI <- exposure(diffnet_FJI, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_FJI, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_FJI.csv")

lagged_exp_nat_count_FJI<- exposure(diffnet_FJI, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_FJI, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_FJI.csv")

diffnet_MNG <- diffnet_new[diffnet_new[["country"]] %in% countries[110]]
exp_nat_count_MNG <- exposure(diffnet_MNG, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_MNG, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_MNG.csv")

lagged_exp_nat_count_MNG<- exposure(diffnet_MNG, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_MNG, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_MNG.csv")

diffnet_BRN <- diffnet_new[diffnet_new[["country"]] %in% countries[107]]
exp_nat_count_BRN <- exposure(diffnet_BRN, valued=FALSE, normalize=FALSE)
write.csv(exp_nat_count_BRN, file = "Data/National exposure/Asia-Pacific (periphery)/Normal/nat_exp_count_BRN.csv")

lagged_exp_nat_count_BRN<- exposure(diffnet_BRN, valued=FALSE, normalize=FALSE, lags = 1L)
write.csv(lagged_exp_nat_count_BRN, file = "Data/National exposure/Asia-Pacific (periphery)/Lagged/lagged_nat_exp_count_BRN.csv")

# 11. Dataframes for modelling ####

# After extracting each network exposure scores (regional per community and country)
# Clean geographical network scores can be calculated by subtracting networks scores, using the bases files
# following the procedure laid out in Chapter section titled: "Computing exposure to geographically bounded networks".
# This can be done directly in Excel. 
# In the folder "Data/Bases" I provide the bases of lagged networks scores for each community. 
# These were manually extracted from the Network Exposure objects created in Step 5.
# The folder "Data/Clean scores" contains cleans scores (Normal and Lagged).

## North America ####

USACAN_int_exp <- read.csv("Data/Clean scores/Lagged/01.lagged_int_exp_North America.csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")
USACAN_reg_exp <- read.csv("Data/Clean scores/Lagged/02.lagged_reg_exp_North America.csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")
USACAN_nat_exp <- read.csv("Data/Clean scores/Lagged/03.lagged_nat_exp_North America.csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")

# Some transformations....
USACAN_int_exp <- data.frame(USACAN_int_exp[,-1], row.names = USACAN_int_exp[,1])
USACAN_reg_exp <- data.frame(USACAN_reg_exp[,-1], row.names = USACAN_reg_exp[,1])
USACAN_nat_exp <- data.frame(USACAN_nat_exp[,-1], row.names = USACAN_nat_exp[,1])

diffnet_NorthAmerica[["lagged_int_expo"]]   <- USACAN_int_exp
diffnet_NorthAmerica[["lagged_reg_expo"]]   <- USACAN_reg_exp
diffnet_NorthAmerica[["lagged_nat_expo"]]   <- USACAN_nat_exp
diffnet_NorthAmerica[["adopted"]]    <- toa_mat(diffnet_NorthAmerica)$cumadopt     

mydata_USACAN <- as.data.frame(diffnet_NorthAmerica)
mydata_USACAN$zebra_cited_per <- as.numeric(as.character(mydata_USACAN$zebra_cited_per))
mydata_USACAN$n_refs <- as.numeric(as.character(mydata_USACAN$n_refs))

# Write data as .csv
write.csv(mydata_USACAN, file = "Data/Clean scores/Lagged/mydata_USACAN_norm.csv")

## LAZEN ####

LAZEN_int_exp <- read.csv("Data/Clean scores/Lagged/04.lagged_int_exp_LAZEN.csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")
LAZEN_reg_exp <- read.csv("Data/Clean scores/Lagged/05.lagged_reg_exp_LAZEN.csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")
LAZEN_nat_exp <- read.csv("Data/Clean scores/Lagged/06.lagged_nat_exp_LAZEN.csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")

LAZEN_int_exp <- data.frame(LAZEN_int_exp[,-1], row.names = LAZEN_int_exp[,1])
LAZEN_reg_exp <- data.frame(LAZEN_reg_exp[,-1], row.names = LAZEN_reg_exp[,1])
LAZEN_nat_exp <- data.frame(LAZEN_nat_exp[,-1], row.names = LAZEN_nat_exp[,1])

diffnet_LAZEN[["lagged_int_expo"]]   <- LAZEN_int_exp
diffnet_LAZEN[["lagged_reg_expo"]]   <- LAZEN_reg_exp
diffnet_LAZEN[["lagged_nat_expo"]]   <- LAZEN_nat_exp
diffnet_LAZEN[["adopted"]]    <- toa_mat(diffnet_LAZEN)$cumadopt

mydata_LAZEN <- as.data.frame(diffnet_LAZEN)
mydata_LAZEN$pub_n <- as.numeric(as.character(mydata_LAZEN$pub_no))

write.csv(mydata_LAZEN, file = "Data/Clean scores/Lagged/mydata_LAZEN_norm.csv")

## EU-COST (core) ####

EU_core_int_exp <- read.csv("Data/Clean scores/Lagged/07.lagged_int_exp_EU-COST(core).csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")
EU_core_reg_exp <- read.csv("Data/Clean scores/Lagged/08.lagged_reg_exp_EU-COST(core).csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")
EU_core_nat_exp <- read.csv("Data/Clean scores/Lagged/09.lagged_nat_exp_EU-COST(core).csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")

EU_core_int_exp <- data.frame(EU_core_int_exp[,-1], row.names = EU_core_int_exp[,1])
EU_core_reg_exp <- data.frame(EU_core_reg_exp[,-1], row.names = EU_core_reg_exp[,1])
EU_core_nat_exp <- data.frame(EU_core_nat_exp[,-1], row.names = EU_core_nat_exp[,1])

diffnet_EU_core_COST[["lagged_int_expo"]]   <- EU_core_int_exp
diffnet_EU_core_COST[["lagged_reg_expo"]]   <- EU_core_reg_exp
diffnet_EU_core_COST[["lagged_nat_expo"]]   <- EU_core_nat_exp
diffnet_EU_core_COST[["adopted"]]    <- toa_mat(diffnet_EU_core_COST)$cumadopt

mydata_EU_core <- as.data.frame(diffnet_EU_core_COST)
write.csv(mydata_EU_core, file = "Data/Clean scores/Lagged/mydata_EU_core_norm.csv")

## EU-COST (periphery) ####

EU_peri_int_exp <- read.csv("Data/Clean scores/Lagged/10.lagged_int_exp_EU-COST(periphery).csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")
EU_peri_reg_exp <- read.csv("Data/Clean scores/Lagged/11.lagged_reg_exp_EU-COST(periphery).csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")
EU_peri_nat_exp <- read.csv("Data/Clean scores/Lagged/12.lagged_nat_exp_EU-COST(periphery).csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")

EU_peri_int_exp <- data.frame(EU_peri_int_exp[,-1], row.names = EU_peri_int_exp[,1])
EU_peri_reg_exp <- data.frame(EU_peri_reg_exp[,-1], row.names = EU_peri_reg_exp[,1])
EU_peri_nat_exp <- data.frame(EU_peri_nat_exp[,-1], row.names = EU_peri_nat_exp[,1])

diffnet_EU_peri_COST[["lagged_int_expo"]]   <- EU_peri_int_exp
diffnet_EU_peri_COST[["lagged_reg_expo"]]   <- EU_peri_reg_exp
diffnet_EU_peri_COST[["lagged_nat_expo"]]   <- EU_peri_nat_exp
diffnet_EU_peri_COST[["adopted"]]    <- toa_mat(diffnet_EU_peri_COST)$cumadopt

mydata_EU_peri <- as.data.frame(diffnet_EU_peri_COST)
write.csv(mydata_EU_peri, file = "Data/Clean scores/Lagged/mydata_EU_peri_norm.csv")

## Asia-Pacific (core) ####

ASOC_core_int_exp <- read.csv("Data/Clean scores/Lagged/13.lagged_int_exp_Asia-Pacific(core).csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")
ASOC_core_reg_exp <- read.csv("Data/Clean scores/Lagged/14.lagged_reg_exp_Asia-Pacific(core).csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")
ASOC_core_nat_exp <- read.csv("Data/Clean scores/Lagged/15.lagged_nat_exp_Asia-Pacific(core).csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")

ASOC_core_int_exp <- data.frame(ASOC_core_int_exp[,-1], row.names = ASOC_core_int_exp[,1])
ASOC_core_reg_exp <- data.frame(ASOC_core_reg_exp[,-1], row.names = ASOC_core_reg_exp[,1])
ASOC_core_nat_exp <- data.frame(ASOC_core_nat_exp[,-1], row.names = ASOC_core_nat_exp[,1])

diffnet_ASIA_OC_core[["lagged_int_expo"]]   <- ASOC_core_int_exp
diffnet_ASIA_OC_core[["lagged_reg_expo"]]   <- ASOC_core_reg_exp
diffnet_ASIA_OC_core[["lagged_nat_expo"]]   <- ASOC_core_nat_exp
diffnet_ASIA_OC_core[["adopted"]]    <- toa_mat(diffnet_ASIA_OC_core)$cumadopt

mydata_ASIA_OC_core <- as.data.frame(diffnet_ASIA_OC_core)
write.csv(mydata_ASIA_OC_core, file = "Data/Clean scores/Lagged/mydata_ASIA_OC_core_norm.csv")

## Asia-Pacific (periphery) ####

ASOC_core_peri_int_exp <- read.csv("Data/Clean scores/Lagged/16.lagged_int_exp_Asia-Pacific(periphery).csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")
ASOC_core_peri_reg_exp <- read.csv("Data/Clean scores/Lagged/17.lagged_reg_exp_Asia-Pacific(periphery).csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")
ASOC_core_peri_nat_exp <- read.csv("Data/Clean scores/Lagged/18.lagged_nat_Asia-Pacific(periphery).csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")

ASOC_core_peri_int_exp <- data.frame(ASOC_core_peri_int_exp[,-1], row.names = ASOC_core_peri_int_exp[,1])
ASOC_core_peri_reg_exp <- data.frame(ASOC_core_peri_reg_exp[,-1], row.names = ASOC_core_peri_reg_exp[,1])
ASOC_core_peri_nat_exp <- data.frame(ASOC_core_peri_nat_exp[,-1], row.names = ASOC_core_peri_nat_exp[,1])

diffnet_ASIA_OC_peri[["lagged_int_expo"]]   <- ASOC_core_peri_int_exp
diffnet_ASIA_OC_peri[["lagged_reg_expo"]]   <- ASOC_core_peri_reg_exp
diffnet_ASIA_OC_peri[["lagged_nat_expo"]]   <- ASOC_core_peri_nat_exp
diffnet_ASIA_OC_peri[["adopted"]]    <- toa_mat(diffnet_ASIA_OC_peri)$cumadopt

mydata_ASIA_OC_peri <- as.data.frame(diffnet_ASIA_OC_peri)
write.csv(mydata_ASIA_OC_peri, file = "Data/Clean scores/Lagged/mydata_ASIA_OC_peri_norm.csv")

## ALL ####

ALL_int_exp <- read.csv("Data/Clean scores/Lagged/00.lagged_int_exp_ALL.csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")
ALL_reg_exp <- read.csv("Data/Clean scores/Lagged/00.lagged_reg_exp_ALL.csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")
ALL_nat_exp <- read.csv("Data/Clean scores/Lagged/00.lagged_nat_exp_ALL.csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ",")

ALL_int_exp <- data.frame(ALL_int_exp[,-1], row.names = ALL_int_exp[,1])
ALL_reg_exp <- data.frame(ALL_reg_exp[,-1], row.names = ALL_reg_exp[,1])
ALL_nat_exp <- data.frame(ALL_nat_exp[,-1], row.names = ALL_nat_exp[,1])
ALL_int_exp <- data.frame(ALL_int_exp[,-23])
ALL_reg_exp <- data.frame(ALL_reg_exp[,-23])
ALL_nat_exp <- data.frame(ALL_nat_exp[,-23])

diffnet_ALL[["lagged_int_expo"]]   <- ALL_int_exp
diffnet_ALL[["lagged_reg_expo"]]   <- ALL_reg_exp
diffnet_ALL[["lagged_nat_expo"]]   <- ALL_nat_exp
diffnet_ALL[["adopted"]]    <- toa_mat(diffnet_ALL)$cumadopt

mydata_ALL <- as.data.frame(diffnet_ALL)
write.csv(mydata_ALL, file = "Data/Clean scores/Lagged/mydata_ALL_norm.csv")

# 12. Classification Tree ####

# mydata_all.csv is all communities, with all networks exposures togeter
mydata_ALL <- as.data.frame(read.csv("Clean scores/Lagged/mydata_ALL.csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ","))

mydata_ALL$adopted <- mydata_ALL$adopted > 0

mydata_ALL$adopted <- mydata_ALL$adopted == 1
mydata_ALL$community <- as.factor(mydata_ALL$community)

## Run regression for the whole dataset ####
deepmod_1 <- rpart(adopted ~ lagged_int_expo + lagged_reg_expo + lagged_nat_expo + mob_share + community + factor(per) + inst_size,
                   data    = mydata_ALL,
                   subset = (per <= toa) & (per < 2000), # For the early years of the diffusion process
                   control = rpart.control(cp = 0.01)
)

deepmod_1

# Classification distribution:
table(mydata_ALL_2$adopted)

# Assess how good the classifications are based on >0.5 threshold....
mean(ifelse((predict(deepmod_1, mydata_ALL_2)>0.5) == mydata_ALL_2$adopted, 
            TRUE, FALSE), na.rm = TRUE)

deepmod_2 <- rpart(adopted ~ lagged_int_expo + lagged_reg_expo + lagged_nat_expo + mob_share + community + factor(per) + inst_size,
                   data    = mydata_ALL_2,
                   subset = (per <= toa) & (per > 2000) & (per != 2017), # Subsequent years....
                   control = rpart.control(cp = 0.01)
)

deepmod_2

mean(ifelse((predict(deepmod_2, mydata_ALL_2)>0.5) == mydata_ALL_2$adopted, 
            TRUE, FALSE), na.rm = TRUE)

prp(deepmod_2)


deepmod_3 <- rpart(adopted ~ lagged_int_expo + lagged_reg_expo + lagged_nat_expo + mob_share + community + factor(per),
                   data    = mydata_ALL_2,
                   subset = (per <= toa) & (per < 2017), # The whole dataset
                   control = rpart.control(cp = 0.01)
)

deepmod_3

mean(ifelse((predict(deepmod_3, mydata_ALL_2)>0.5) == mydata_ALL_2$adopted, 
            TRUE, FALSE), na.rm = TRUE)

prp(deepmod_3)


## Run models for each Community ####

### North America ####
mydata_USACAN_2 <- as.data.frame(read.csv("Data/Clean scores/Lagged/mydata_USACAN_norm.csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ","))
mydata_USACAN_2 <- data.frame(mydata_USACAN_2[,-1], row.names = mydata_USACAN_2[,1])

mydata_USACAN_2$adopted <- mydata_USACAN_2$adopted > 0

###  LAZEN ####
mydata_LAZEN_2 <- as.data.frame(read.csv("Data/Clean scores/Lagged/mydata_LAZEN_norm.csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ","))
mydata_LAZEN_2 <- data.frame(mydata_LAZEN_2[,-1], row.names = mydata_LAZEN_2[,1])

mydata_LAZEN_2$adopted <- mydata_LAZEN_2$adopted > 0

### EU-COST (core) ####
mydata_EU_core_2 <- as.data.frame(read.csv("Data/Clean scores/Lagged/mydata_EU_core_norm.csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ","))
mydata_EU_core_2 <- data.frame(mydata_EU_core_2[,-1], row.names = mydata_EU_core_2[,1])

mydata_EU_core_2$adopted <- mydata_EU_core_2$adopted > 0

### EU-COST (periphery]) ####
mydata_EU_peri_2 <- as.data.frame(read.csv("Data/Clean scores/Lagged/mydata_EU_peri_norm.csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ","))
mydata_EU_peri_2 <- data.frame(mydata_EU_peri_2[,-1], row.names = mydata_EU_peri_2[,1])

mydata_EU_peri_2$adopted <- mydata_EU_peri_2$adopted > 0

### ASOC (core) ####
mydata_ASOC_core_2 <- as.data.frame(read.csv("Data/Clean scores/Lagged/mydata_ASIA_OC_core_norm.csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ","))
mydata_ASOC_core_2 <- data.frame(mydata_ASOC_core_2[,-1], row.names = mydata_ASOC_core_2[,1])

mydata_ASOC_core_2$adopted <- mydata_ASOC_core_2$adopted > 0

### ASOC (periphery)  ####
mydata_ASOC_per_2 <- as.data.frame(read.csv("Data/Clean scores/Lagged/mydata_ASIA_OC_peri_norm.csv",header=TRUE,stringsAsFactors=FALSE, check.names = FALSE, sep = ","))
mydata_ASOC_per_2 <- data.frame(mydata_ASOC_per_2[,-1], row.names = mydata_ASOC_per_2[,1])

mydata_ASOC_per_2$adopted <- mydata_ASOC_per_2$adopted > 0
colnames(mydata_ASOC_per_2)

### Model ####
deepmod <- rpart(adopted ~ lagged_int_expo + lagged_reg_expo + lagged_nat_expo + mob_share + factor(per),
                 data    = mydata_USACAN_2, # Replace the dataset with each community dataframe
                 subset =  (per <= toa) & per < 2017,
                 control = rpart.control(cp = 0.01)
)


deepmod

# Plot
rpart.plot(deepmod, type = 3, clip.right.labs = FALSE, branch = .3, under = TRUE)

# Get Tree rules
rules <- as.data.frame(rpart.rules(deepmod, cover = TRUE))
View(rules)

# Write rules as .csv
write.csv(rules, file = "Clean scores/Lagged/results_reg_tree_ASOC_peri.csv")

# Model assessment (i.e.goodness)
mean(ifelse((predict(deepmod, mydata_USACAN_2)>0.5) == mydata_USACAN_2$adopted, 
            TRUE, FALSE), na.rm = TRUE)

