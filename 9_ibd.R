library(gmt)

# no scientific notation
options(scipen=999)

# read in fst stats
fst <- read.table("aoudad_stats.txt", header=T)
fst <- fst[fst$stat == "Fst",]

# read in localities
loc <- read.table("aoudad_localities_texas.txt", header=T)

# subset fst comparisons to those in texas (and not the wildlife center)
fst <- fst[fst$pop1 %in% loc$Locality,]
fst <- fst[fst$pop2 %in% loc$Locality,]

# calculate geographic distances between Texas localities
distance_km <- c()
for(a in 1:nrow(fst)) {
	distance_km <- c(distance_km, geodist(loc$Lat[loc$Locality == fst$pop1[a]], loc$Long[loc$Locality == fst$pop1[a]], loc$Lat[loc$Locality == fst$pop2[a]], loc$Long[loc$Locality == fst$pop2[a]], units="km"))
}

fst <- cbind(fst, distance_km)

plot(fst$distance_km, fst$calculated_stat, pch=19, cex=0.8, xlab="Geographic Distance (km)", ylab="FST", ylim=c(0,0.2), xlim=c(0,700))
abline(lm(fst$calculated_stat ~ fst$distance_km))
summary(lm(fst$calculated_stat ~ fst$distance_km))