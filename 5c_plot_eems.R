

library(rEEMSplots)
library(raster)
library(rgdal)

mcmcpath = "./output"
plotpath = "./output/plots"
sr <- "+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

projection_none <- "+proj=longlat +datum=WGS84"
projection <- sr
usa <- readOGR("gadm36_USA_shp", "gadm36_USA_1")
usa <- subset(usa, usa@data$NAME_1 %in% c("Texas"))


eemsResults <- rEEMSplots::eems.plots(mcmcpath, plotpath, longlat = T,out.png = FALSE, projection.in = projection_none, projection.out = projection_mercator)

eems.plots(mcmcpath, plotpath, longlat = F, out.png = FALSE,add.grid=T, add.demes=T, projection.in = projection_none, projection.out = projection, m.plot.xy = { plot(usa, col = NA, add = TRUE) }, q.plot.xy = { plot(usa, col = NA, add = TRUE) })

eems.plots(mcmcpath, plotpath, longlat = F, out.png = FALSE,add.grid=F, add.demes=T, projection.in = projection_none, projection.out = projection, m.plot.xy = { plot(usa, col = NA, add = TRUE) }, q.plot.xy = { plot(usa, col = NA, add = TRUE) })

