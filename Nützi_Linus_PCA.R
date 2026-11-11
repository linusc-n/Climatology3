"""
Linus Nützi
Climatology 3, exercise 2, Correlations
10.11.26
"""

setwd("C:/Users/linus/OneDrive/Desktop/Uni/HS26/Geographie/Klima 3/Data/")


library(ncdf4)
library(maps)
library(abind)


####### READ
### Read netcdf files 
f1 <- nc_open("ModE-RA_ensmean_slp_abs_1420-2009.nc")
f2 <- nc_open("ModE-RA_ensmean_temp2_abs_1420-2009.nc")

lon <- f1$dim[[2]]$vals
lat <- f1$dim[[3]]$vals
time <- f1$dim[[1]]$vals

### Read and prepare variables
slp <- ncvar_get(f1,varid="slp")
temp <- ncvar_get(f2,varid="temp2")

slp_2 <- array(slp,dim=c(length(lon),length(lat),12,length(time)/12))
slp_3 <- aperm(slp_2,c(1,2,4,3))

temp_2 <- array(temp,dim=c(length(lon),length(lat),12,length(time)/12))
temp_3 <- aperm(temp_2,c(1,2,4,3))


### extended boreal winter Nov-April
slp_4 <- abind(slp_3[,,-(length(time)/12),11:12],slp_3[,,-1,1:4],along=4)
slp.ndjfma <- apply(slp_4,c(1:3),mean)

temp_4 <- abind(temp_3[,,-(length(time)/12),11:12],temp_3[,,-1,1:4],along=4)
temp.ndjfma <- apply(temp_4,c(1:3),mean)

yr.slp <- c(1422:2008)
yr.temp <- c(1422:2008)

selyrs_slp <- which(yr.slp >= 1752 & yr.slp <= 1901)
selyrs_temp <- which(yr.temp >= 1781 & yr.temp <= 2001)

#years to label timeseries x-axis
selyrs_slp_nonidx <- yr.slp[selyrs_slp]
selyrs_temp_nonidx <- yr.temp[selyrs_temp]


### Adjust for regions
idx_lat_nord <- which(lat >= 20)
lat_nord <- lat[idx_lat_nord]

idx_lat_trop <- which(lat >= -20 & lat <= 20)
lat_trop <- lat[idx_lat_trop]

slp_nord <- slp.ndjfma[ , idx_lat_nord, selyrs_slp]
temp_trop <- temp.ndjfma[ , idx_lat_trop, selyrs_temp]
 

#### Prepare PCA
y_slp.scaled <- apply(slp_nord,c(1,2),scale)
y_temp.scaled <- apply(temp_trop,c(1,2),scale)

y_slp.weighted <- apply(y_slp.scaled,c(1,2),"*",sqrt(cos(3.14159*lat_nord/180))) 
y_temp.weighted <- apply(y_temp.scaled,c(1,2),"*",sqrt(cos(3.14159*lat_trop/180))) 

y_slp.weighted[1,1,1]
y_slp.perm <- aperm(y_slp.weighted,c(2,3,1))
y_slp.mat <- matrix(y_slp.perm,nrow=length(selyrs_slp),ncol=length(lon)*length(lat_nord))

y_temp.weighted[1,1,1]
y_temp.perm <- aperm(y_temp.weighted,c(2,3,1))
y_temp.mat <- matrix(y_temp.perm,nrow=length(selyrs_temp),ncol=length(lon)*length(lat_trop))

### PCA
#inverted SLP-data to get negative values over iceland and positive values over azores
pc_slp <- prcomp(y_slp.mat,retx=T,scale.=F)
plotpc_slp <- -matrix(pc_slp$rotation[,1],nrow=length(lon),ncol=length(lat_nord))
plottime_slp <- -pc_slp$x[,1]

pc_temp <- prcomp(y_temp.mat,retx=T,scale.=F)
plotpc_temp <- matrix(pc_temp$rotation[,1],nrow=length(lon),ncol=length(lat_trop))
plottime_temp <- pc_temp$x[,1]

######## Plotting
### time series
plot(selyrs_slp_nonidx,plottime_slp,type="l", main = "PCA Timeseries Nov-April Mean Sea level pressure (1751-1901)", xlab="year", ylab="PC1-Score")
plot(selyrs_temp_nonidx,plottime_temp,type="l", main = "PCA Timeseries Nov-April Mean Temperature (1780-2001)", xlab="year", ylab="PC1-Score")


### field
# revert lat
x <- lon
y_slp <- rev(lat_nord)
y_temp <- rev(lat_trop)
z_slp <- plotpc_slp[,rev(1:length(lat_nord))]
z_temp <- plotpc_temp[,rev(1:length(lat_trop))]

# plot
mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")

mylevs_slp <- max(max(z_slp),abs(min(z_slp)))*(c(0:10)-5)/5
filled.contour(x, y_slp, z_slp, levels = mylevs_slp, col = mycol, plot.title = title(main = "PCA Nov-April Mean Sea level pressure (1751-1901)", xlab = "Longitude", ylab = "Latitude"),
               plot.axes = {
                 axis(1, at = c(-180, -120, -60, 0, 60, 120, 180))
                 axis(2, at = c(20, 30, 40, 50, 60, 70, 80, 90))
                 map("world", interior = FALSE, add = TRUE)
               }
)

mylevs_temp <- max(max(z_temp),abs(min(z_temp)))*(c(0:10)-5)/5
filled.contour(x, y_temp, z_temp, levels = mylevs_temp, col = mycol, plot.title = title(main = "PCA Nov-April Mean Temperature (1780-2001)", xlab = "Longitude", ylab = "Latitude"),
  plot.axes = {
    axis(1, at = c(-180, -120, -60, 0, 60, 120, 180))
    axis(2, at = c(-20, -10, 0, 10, 20))
    map("world", interior = FALSE, add = TRUE)
  }
)
