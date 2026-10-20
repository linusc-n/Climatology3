"""
Linus Nützi
Climatology 3, exercise 1, ENSO-Composites
19.10.26
"""

setwd("C:/Users/linus/OneDrive/Desktop/Uni/HS26/Geographie/Klima 3/Data/")

library(ncdf4)
library(maps)
library(abind)


####### READ
### Read netcdf files and check
f1 <- nc_open("ModE-RA_ensmean_temp2_abs_1420-2009.nc")
f2 <- nc_open("ModE-RA_ensmean_geopoth_50000_abs_1420-2009.nc")
f1$dim
f2$dim
length(f1$dim)
length(f2$dim)
f1$dim[[1]]
f2$dim[[1]]
f1$dim[[1]]$name
f2$dim[[1]]$name
f1$dim[[1]]$units
f2$dim[[1]]$units
f1$dim[[2]]$units
f2$dim[[2]]$units


####### PREPROCESS
### set dimensions
time <- f1$dim[[1]]$vals ## note that unit is "months since 1421-01-15 00:00:00"
lon <- f1$dim[[2]]$vals
lat <- f1$dim[[3]]$vals

## read entire variable
temp <- ncvar_get(f1,varid="temp2")
z500 <- ncvar_get(f2,varid="geopotential_height")
dim(temp)
dim(z500)

### generate array: lon, lat, year, mon
temp_2 <- array(temp,dim=c(length(lon),length(lat),12,length(time)/12))
dim(temp_2)
z500_2 <- array(z500,dim=c(length(lon),length(lat),12,length(time)/12))
dim(z500_2)

temp_3 <- aperm(temp_2,c(1,2,4,3))
dim(temp_3)
z500_3 <- aperm(z500_2,c(1,2,4,3))
dim(z500_3)

### create boreal winter (Jan-April) mean values
temp.jfma <- apply(temp_3[,,,1:4],c(1:3),mean)
dim(temp.jfma)
z500.jfma <- apply(z500_3[,,,1:4],c(1:3),mean)
dim(z500.jfma) 

### new time axes (yr)
yr.jfma <- c(1421:2008)
length(yr.jfma)


######## ANALYSIS
### Lists of years with El Nino and years with La Nina
elnino <- scan("ElNinoYears.txt", skip=1)
lanina <- scan("LaNinaYears.txt", skip=1)

### get years with events (year 0) and following years (year 1)
selelnino_0 <- match(elnino,yr.jfma)
selelnino_0
selelnino_1 <- match(elnino+1,yr.jfma)

sellanina_0 <- match(lanina,yr.jfma)
sellanina_1 <- match(lanina+1,yr.jfma)

### composites averaging over selected years for both El Nino and La Nina in T and Z500
dim(temp.jfma[,,selelnino_0])
comp_temp_elnino_0 <- apply(temp.jfma[,,selelnino_0],c(1:2),mean)
dim(z500.jfma[,,selelnino_0])
comp_z500_elnino_0 <- apply(z500.jfma[,,selelnino_0],c(1:2),mean)

dim(temp.jfma[,,selelnino_1])
comp_temp_elnino_1 <- apply(temp.jfma[,,selelnino_1],c(1:2),mean)
dim(z500.jfma[,,selelnino_1])
comp_z500_elnino_1 <- apply(z500.jfma[,,selelnino_1],c(1:2),mean)


dim(temp.jfma[,,sellanina_0])
comp_temp_lanina_0 <- apply(temp.jfma[,,sellanina_0],c(1:2),mean)
dim(z500.jfma[,,sellanina_0])
comp_z500_lanina_0 <- apply(z500.jfma[,,sellanina_0],c(1:2),mean)

dim(temp.jfma[,,sellanina_1])
comp_temp_lanina_1 <- apply(temp.jfma[,,sellanina_1],c(1:2),mean)
dim(z500.jfma[,,sellanina_1])
comp_z500_lanina_1 <- apply(z500.jfma[,,sellanina_1],c(1:2),mean)

### reference period: Mean of all years; 5 years prior events
ref_temp <- apply(temp.jfma,c(1:2),mean)
ref_z500 <- apply(z500.jfma,c(1:2),mean)

ref5years_elnino <- c(elnino-1,elnino-2,elnino-3,elnino-4,elnino-5)
selref_elnino <- match(ref5years_elnino,yr.jfma)
selref_elnino

ref5years_lanina <- c(lanina-1,lanina-2,lanina-3,lanina-4,lanina-5)
selref_lanina <- match(ref5years_lanina,yr.jfma)
selref_lanina

### testing for NA-Indizes; if >0, NA-Index in data
sum(is.na(selref_elnino))
sum(is.na(selref_lanina))
sum(is.na(selelnino_1))
sum(is.na(sellanina_1))

### correction and application
selref_elnino <- match(ref5years_elnino, yr.jfma) 
selref_elnino <- unique(selref_elnino[!is.na(selref_elnino)])

ref5years_temp_elnino <- apply(temp.jfma[,,selref_elnino],c(1,2),mean)
ref5years_z500_elnino <- apply(z500.jfma[,,selref_elnino],c(1,2),mean)

ref5years_temp_lanina <- apply(temp.jfma[,,selref_lanina],c(1,2),mean)
ref5years_z500_lanina <- apply(z500.jfma[,,selref_lanina],c(1,2),mean)

### differences to get anomalies
diff0_temp_elnino_0 <- comp_temp_elnino_0-ref_temp
diff0_z500_elnino_0 <- comp_z500_elnino_0-ref_z500

diff1_temp_elnino_1 <- comp_temp_elnino_1-ref5years_temp_elnino
diff1_z500_elnino_1 <- comp_z500_elnino_1-ref5years_z500_elnino


diff0_temp_lanina_0 <- comp_temp_lanina_0-ref_temp
diff0_z500_lanina_0 <- comp_z500_lanina_0-ref_z500

diff1_temp_lanina_1 <- comp_temp_lanina_1-ref5years_temp_lanina
diff1_z500_lanina_1 <- comp_z500_lanina_1-ref5years_z500_lanina


diff0_temp_enln_0 <- comp_temp_elnino_0-comp_temp_lanina_0
diff0_z500_enln_0 <- comp_z500_elnino_0-comp_z500_lanina_0

diff1_temp_enln_1 <- comp_temp_elnino_1-comp_temp_lanina_1
diff1_z500_enln_1 <- comp_z500_elnino_1-comp_z500_lanina_1

######## PLOT
### revert lat for y axis and colourscale
x <- lon
lat
y <- rev(lat)

mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")

### preparation for for-loop
vars <- list(
  diff0_temp_elnino_0 = diff0_temp_elnino_0,
  diff0_z500_elnino_0 = diff0_z500_elnino_0,
  diff1_temp_elnino_1 = diff1_temp_elnino_1,
  diff1_z500_elnino_1 = diff1_z500_elnino_1,
  diff0_temp_lanina_0 = diff0_temp_lanina_0,
  diff0_z500_lanina_0 = diff0_z500_lanina_0,
  diff1_temp_lanina_1 = diff1_temp_lanina_1,
  diff1_z500_lanina_1 = diff1_z500_lanina_1,
  diff0_temp_enln_0   = diff0_temp_enln_0,
  diff0_z500_enln_0   = diff0_z500_enln_0,
  diff1_temp_enln_1   = diff1_temp_enln_1,
  diff1_z500_enln_1   = diff1_z500_enln_1
)

titel <- c(
  "Temperature Anomaly: Years with El Nino",
  "Z500 Anomaly: Years with El Nino",
  "Temperature Anomaly: Years prior to El Nino",
  "Z500 Anomaly: Years prior to El Nino",
  "Temperature Anomaly: Years with La Nina",
  "Z500 Anomaly: Years with La Nina",
  "Temperature Anomaly: Years prior to La Nina",
  "Z500 Anomaly: Years prior to La Nina",
  "Temperature Anomaly: Comparison Years with El Nino vs La Nina",
  "Z500 Anomaly: Comparison Years with El Nino vs La Nina",
  "Temperature Anomaly: Comparison Years before El Nino vs La Nina",
  "Z500 Anomaly: Comparison Years before El Nino vs La Nina"
)

###plot all at once to avoid manual inputs
library(maps)

for (i in seq_along(vars)) {
  name   <- names(vars)[i]
  mydata <- vars[[i]][, rev(1:length(lat))]
  
  maxabs <- max(max(mydata, na.rm = TRUE), abs(min(mydata, na.rm = TRUE)))
  mylevs <- maxabs * (c(0:10) - 5) / 5
  
  png(paste0(name, ".png"), width = 1200, height = 800, res = 150)
  
  filled.contour(
    x, y, mydata,
    levels = mylevs,
    col = mycol,
    plot.axes = {
      axis(1)
      axis(2)
      map("world", interior = FALSE, add = TRUE)
    },
    xlab = "Longitude",
    ylab = "Latitude",
    main = titel[i]
  )
  
  dev.off()
}



# Close the NetCDF connection
nc_close(f1)
