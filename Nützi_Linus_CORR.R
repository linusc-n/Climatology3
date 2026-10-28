"""
Linus Nützi
Climatology 3, exercise 2, Correlations
26.10.26
"""

setwd("C:/Users/linus/OneDrive/Desktop/Uni/HS26/Geographie/Klima 3/Data/")


library(ncdf4)
library(maps)
library(abind)


####### READ
### Read netcdf files 
f1 <- nc_open("ModE-RA_ensmean_slp_abs_1420-2009.nc")
f2 <- nc_open("ModE-RA_ensmean_temp2_abs_1420-2009.nc")
f3 <- nc_open("ModE-RA_ensmean_totprec_abs_1420-2009.nc")


####### PREPROCESS
### read dimensions 
time <- f1$dim[[1]]$vals
lon <- f1$dim[[2]]$vals
lat <- f1$dim[[3]]$vals

### read entire variables
slp <- ncvar_get(f1,varid="slp")
temp <- ncvar_get(f2, varid="temp2")
totprec <- ncvar_get(f3, varid="totprec")

### generate array: lon, lat, year, mon
slp_2 <- array(slp,dim=c(length(lon),length(lat),12,length(time)/12))
temp_2 <- array(temp,dim=c(length(lon),length(lat),12,length(time)/12))
totprec_2 <- array(totprec,dim=c(length(lon),length(lat),12,length(time)/12))

slp_3 <- aperm(slp_2,c(1,2,4,3))
temp_3 <- aperm(temp_2,c(1,2,4,3))
totprec_3 <- aperm(totprec_2,c(1,2,4,3))

### create winter (Dec-Feb) mean values
slp_4 <- abind(slp_3[,,-(length(time)/12),12],slp_3[,,-1,1:2],along=4)
slp.djf <- apply(slp_4[,,,1:3],c(1:3),mean)

### create summer (June-Aug) mean values for temperature and precipitation
temp.jja <- apply(temp3[,,,6:8],c(1:3),mean)
prec.jja <- apply(prec3[,,,6:8],c(1:3),mean)

### do new time axes (yr)
yr.djf <- c(1422:2008) # because of the shifting with Dec
yr.jja <- c(1421:2008)

selyrs_slp <- which(yr.djf >= 1709 & yr.djf <= 2008)
selyrs_temp <- which(yr.jja >= 1870 & yr.jja <= 2008)

### Coordinates gridcell: Reykjavik: 338.5E, 64.1N; Lisbon: 350.8, 38.7N; Lugano: 9E, 46N; Northpole: 0E, 90N; Central Pacific Ocean: 145W = -145E, 0N)
clon <- -145     #change manually
clon <- ifelse(clon > 180, clon -360, clon)
clat <- 0     #change manually

lon_index=which.min(abs(lon-clon))
lat_index=which.min(abs(lat-clat))

### selecting gridcell: 
location_slp <- slp.djf[lon_index,lat_index,selyrs_slp]
plotfield_slp <- apply(slp.djf[,,selyrs_slp],c(1,2),cor,location_slp)

location_temp <- temp.jja[lon_index,lat_index,selyrs_temp]
plotfield_temp <- apply(prec.jja[,,selyrs_temp],c(1,2),cor,location_temp) # same argument is used for every grid cell


###### Plotting
### Prepare plotting
x <- lon
y <- rev(lat)
#z <- plotfield_slp[,rev(1:length(lat))]
z <- plotfield_temp[,rev(1:length(lat))]

### plot
mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- max(max(z),abs(min(z)))*(c(0:10)-5)/5
#filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, plot.title = title(main = "SLP Correlation of the Northpole for the period 1709-2008 (Dec-Feb)"))
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, plot.title = title(main = "Temperature Correlation of the Central Pacific for the period 1870-2008 (Jun-Aug)"))

