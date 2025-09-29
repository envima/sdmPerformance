
l=list.files("data/paRaster/",pattern="RDS", full.names = F)

for (i in l){
  vs=readRDS(paste0("data/paRaster/",i))
  r=terra::unwrap(vs$probability.of.occurrence)
  terra::writeRaster(r,  paste0("data/virtualSpecies/", gsub(".RDS","",i),".tif"));rm(vs)
}
