#'@name 01_virtualspecies.R
#'@date 28.02.2025
#'@author Lisa Bald [bald@staff.uni-marburg.de]
#'@author 
#'@description train a model for each virtual species
#'@description code taken from the benchmark study of Valavi et al. 2023 DOI: https://doi.org/10.1111/geb.13639
#'@misc OSF repository of the study of Valavi et al. 2023: DOI: 10.17605/OSF.IO/G6DC3
#'@source https://osf.io/puk8v (of the code for the models)


library(virtualspecies)
library(tidyverse)

setwd("M:/user/bald/SDM/sdmPerformance/")

#https://www.sciencedirect.com/science/article/pii/S0304380020302659

if (!file.exists("data/variables.tif")){
  
  if(!dir.exists("data")) dir.create("data")
  bioclim=geodata::worldclim_country(country="Australia", path="data", var="bio", res=0.5)
  
  names(bioclim) <-  substr(names(bioclim), start=11, stop=16)
  
  bioclim=terra::subset(bioclim, c("bio_1", "bio_3", "bio_7", "bio_12"))
  
  
  
  border=geodata::gadm(country="Australia", path="data")
  border=sf::st_as_sf(border)
  border=border%>%dplyr::filter(NAME_1 %in% c("New South Wales", "Victoria", "Australian Capital Territory"))
  border=sf::st_transform(border, terra::crs(bioclim))
  
  bioclim=terra::crop(bioclim, border)
  bioclim=terra::mask(bioclim, border)
  bioclim=terra::project(bioclim, "epsg:3577")
  
  terra::writeRaster(bioclim, "data/variables.tif", overwrite=T)
  rm(border, bioclim)
}
# 2 - create virtual species ####
#-------------------------------#


r=terra::rast("data/variables.tif")
if(!dir.exists("data/virtualSpecies")) dir.create("data/virtualSpecies")


# create VS01 ####
#----------------#
paramsVS01 <- formatFunctions(bio_1 = c(fun = 'dnorm', mean = 18, sd = 4),
                              bio_3 = c(fun = 'dnorm', mean = 43, sd = 3),
                              bio_7 = c(fun = 'dnorm', mean = 30, sd = 5),
                              bio_12 = c(fun = 'dnorm', mean = 400, sd = 200))
VS01 <- virtualspecies::generateSpFromFun(raster.stack = r, parameters = paramsVS01, formula="bio_1 * bio_3 * bio_7 * bio_12")
terra::plot(VS01$suitab.raster)
saveRDS(VS01, paste0("data/virtualSpecies/", "VS01.RDS"))
terra::writeRaster(VS01$suitab.raster,  paste0("data/virtualSpecies/", "VS01.tif"));rm(VS01,paramsVS01)

# create VS02 ####
#----------------#

paramsVS02 <- formatFunctions(bio_1 = c(fun = 'dnorm', mean = 15, sd = 5),
                              bio_3 = c(fun = 'dnorm', mean = 45, sd = 5),
                              bio_7 = c(fun = 'dnorm', mean = 30, sd = 5),
                              bio_12 = c(fun = 'dnorm', mean = 500, sd = 200))
VS02 <- virtualspecies::generateSpFromFun(raster.stack = r, parameters = paramsVS02, formula="bio_1 * bio_3 * bio_7 * bio_12")
terra::plot(VS02$suitab.raster)
saveRDS(VS02, paste0("data/virtualSpecies/", "VS02.RDS"))
terra::writeRaster(VS02$suitab.raster,  paste0("data/virtualSpecies/", "VS02.tif"));rm(VS02,paramsVS02)

# create VS03 ####
#----------------#
paramsVS03 <- formatFunctions(bio_1 = c(fun = 'dnorm', mean = 18, sd = 4),
                              bio_3 = c(fun = 'dnorm', mean = 45, sd = 2),
                              bio_7 = c(fun = 'dnorm', mean = 30, sd = 5),
                              bio_12 = c(fun = 'dnorm', mean = 300, sd = 200))
VS03 <- virtualspecies::generateSpFromFun(raster.stack = r, parameters = paramsVS03, formula="bio_1 * bio_3 * bio_7 * bio_12")
terra::plot(VS03$suitab.raster)
saveRDS(VS03, paste0("data/virtualSpecies/", "VS03.RDS"))
terra::writeRaster(VS03$suitab.raster,  paste0("data/virtualSpecies/", "VS03.tif"));rm(VS03,paramsVS03)


# create VS04 ####
#----------------#
paramsVS04 <- formatFunctions(bio_1 = c(fun = 'dnorm', mean = 18, sd = 4),
                              bio_3 = c(fun = 'dnorm', mean = 45, sd = 4),
                              bio_7 = c(fun = 'dnorm', mean = 30, sd = 5),
                              bio_12 = c(fun = 'dnorm', mean = 700, sd = 300))
VS04 <- virtualspecies::generateSpFromFun(raster.stack = r, parameters = paramsVS04, formula="bio_1 * bio_3 * bio_7 * bio_12")
terra::plot(VS04$suitab.raster)
saveRDS(VS04, paste0("data/virtualSpecies/", "VS04.RDS"))
terra::writeRaster(VS04$suitab.raster,  paste0("data/virtualSpecies/", "VS04.tif"));rm(VS04,paramsVS04)


# create VS05 ####
#----------------#
paramsVS05 <- formatFunctions(bio_1 = c(fun = 'dnorm', mean = 13, sd = 7),
                              bio_3 = c(fun = 'dnorm', mean = 43, sd = 6),
                              bio_7 = c(fun = 'dnorm', mean = 25, sd = 7),
                              bio_12 = c(fun = 'dnorm', mean = 1000, sd = 600))
VS05 <- virtualspecies::generateSpFromFun(raster.stack = r, parameters = paramsVS05, formula="bio_1 * bio_3 * bio_7 * bio_12")
terra::plot(VS05$suitab.raster)
saveRDS(VS05, paste0("data/virtualSpecies/", "VS05.RDS"))
terra::writeRaster(VS05$suitab.raster,  paste0("data/virtualSpecies/", "VS05.tif"));rm(VS05,paramsVS05)


# create VS06 ####
#----------------#
paramsVS06 <- formatFunctions(bio_1 = c(fun = 'dnorm', mean = 18, sd = 3),
                              bio_3 = c(fun = 'dnorm', mean = 50, sd = 8),
                              bio_7 = c(fun = 'dnorm', mean = 27, sd = 10),
                              bio_12 = c(fun = 'dnorm', mean = 850, sd = 500))
VS06 <- virtualspecies::generateSpFromFun(raster.stack = r, parameters = paramsVS06, formula="bio_1 * bio_3 * bio_7 * bio_12")
terra::plot(VS06$suitab.raster)
saveRDS(VS06, paste0("data/virtualSpecies/", "VS06.RDS"))
terra::writeRaster(VS06$suitab.raster,  paste0("data/virtualSpecies/", "VS06.tif"));rm(VS06,paramsVS06)


# create VS07 ####
#----------------#
paramsVS07 <- formatFunctions(bio_1 = c(fun = 'dnorm', mean = 15, sd = 5),
                              bio_3 = c(fun = 'dnorm', mean = 49, sd = 5),
                              bio_7 = c(fun = 'dnorm', mean = 22, sd = 5),
                              bio_12 = c(fun = 'dnorm', mean = 850, sd = 300))
VS07 <- virtualspecies::generateSpFromFun(raster.stack = r, parameters = paramsVS07, formula="bio_1 * bio_3 * bio_7 * bio_12")
terra::plot(VS07$suitab.raster)
saveRDS(VS07, paste0("data/virtualSpecies/", "VS07.RDS"))
terra::writeRaster(VS07$suitab.raster,  paste0("data/virtualSpecies/", "VS07.tif"));rm(VS07,paramsVS07)


# create VS08 ####
#----------------#
paramsVS08 <- formatFunctions(bio_1 = c(fun = 'dnorm', mean = 15, sd = 2),
                              bio_3 = c(fun = 'dnorm', mean = 45, sd = 4),
                              bio_7 = c(fun = 'dnorm', mean = 29, sd = 5),
                              bio_12 = c(fun = 'dnorm', mean = 800, sd = 300))
VS08 <- virtualspecies::generateSpFromFun(raster.stack = r, parameters = paramsVS08, formula="bio_1 * bio_3 * bio_7 * bio_12")
terra::plot(VS08$suitab.raster)
saveRDS(VS08, paste0("data/virtualSpecies/", "VS08.RDS"))
terra::writeRaster(VS08$suitab.raster,  paste0("data/virtualSpecies/", "VS08.tif"));rm(VS08,paramsVS08)


# create VS09 ####
#----------------#
paramsVS09 <- formatFunctions(bio_1 = c(fun = 'dnorm', mean = 20, sd = 2),
                              bio_3 = c(fun = 'dnorm', mean = 45, sd = 2),
                              bio_7 = c(fun = 'dnorm', mean = 32, sd = 5),
                              bio_12 = c(fun = 'dnorm', mean = 350, sd = 100))
VS09 <- virtualspecies::generateSpFromFun(raster.stack = r, parameters = paramsVS09, formula="bio_1 * bio_3 * bio_7 * bio_12")
terra::plot(VS09$suitab.raster)
saveRDS(VS09, paste0("data/virtualSpecies/", "VS09.RDS"))
terra::writeRaster(VS09$suitab.raster,  paste0("data/virtualSpecies/", "VS09.tif"));rm(VS09,paramsVS09)



# create VS10 ####
#----------------#
paramsVS10 <- formatFunctions(bio_1 = c(fun = 'dnorm', mean = 14, sd = 4),
                              bio_3 = c(fun = 'dnorm', mean = 49, sd = 3),
                              bio_7 = c(fun = 'dnorm', mean = 22, sd = 4),
                              bio_12 = c(fun = 'dnorm', mean = 500, sd = 300))
VS10 <- virtualspecies::generateSpFromFun(raster.stack = r, parameters = paramsVS10, formula="bio_1 * bio_3 * bio_7 * bio_12")
terra::plot(VS10$suitab.raster)
saveRDS(VS10, paste0("data/virtualSpecies/", "VS10.RDS"))
terra::writeRaster(VS10$suitab.raster,  paste0("data/virtualSpecies/", "VS10.tif"));rm(VS10,paramsVS10)



# 3 - Conversion to PA data #####
#-------------------------------#

if(!dir.exists("data/paRaster/")) dir.create("data/paRaster")

vs=list.files("data/virtualSpecies/", pattern=".tif", full.names=T)

prevalence=data.frame(species=1:10, prevalence=c(0.35,0.34,0.33,0.29,0.26,0.21,0.15,0.12,0.11,0.05))

for(v in 1:length(vs)){
  # Conversion to presence-absence
  species=terra::rast(vs[v])
  PA <- virtualspecies::convertToPA(species, plot = F, species.prevalence = prevalence$prevalence[v])
  saveRDS(PA, paste0("data/paRaster/", gsub(".tif","",strsplit(vs[v], "/")[[1]][3]), ".RDS"))
  PA=terra::unwrap(PA$pa.raster)
  terra::writeRaster(PA, paste0("data/paRaster/", strsplit(vs[v], "/")[[1]][3]))
  rm(PA, species)
}

# 4 - Sample presence absence data ####
#-------------------------------------#

if(!dir.exists("data/PA")) dir.create("data/PA")
#sample pa data 20 times (replicates)
grid=expand.grid(vs=gsub(".tif","",list.files("data/paRaster/", pattern=".tif", full.names=F)) ,
            points=c(40,80,120,160,200,400),
            replicates=1:10)

set.seed(12345)
for (i in 1:nrow(grid)){
  PA_raster=terra::rast(paste0("data/paRaster/",grid$vs[i],".tif"))
  PA=virtualspecies::sampleOccurrences(x=PA_raster, 
                                       n=grid$points[i],
                                       type="presence-absence",
                                       sample.prevalence = 0.75)
  saveRDS(PA, paste0("data/PA/",grid$vs[i],"_",grid$points[i],"_",grid$replicates[i],".RDS"))
  PA=PA$sample.points
  PA=sf::st_as_sf(PA, coords=c("x","y"), crs=terra::crs(PA_raster))
  sf::write_sf(PA, paste0("data/PA/",grid$vs[i],"_",grid$points[i],"_",grid$replicates[i],".gpkg"))
  rm(PA,PA_raster)
}


