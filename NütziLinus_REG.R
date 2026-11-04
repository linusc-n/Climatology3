"""
Linus Nützi
Climatology 3, exercise 3, Regressions
2.11.26
"""

setwd("C:/Users/linus/OneDrive/Desktop/Uni/HS26/Geographie/Klima 3/Data/")

library(ncdf4)
library(maps)

### regression
# read the data
f1 <- nc_open("ModE-RA_ensmean_geopoth_50000_abs_1420-2009.nc")
f2 <- nc_open("ModE-RA_ensmean_temp2_abs_1420-2009.nc")
f3 <- nc_open("ModE-RA_ensmean_totprec_abs_1420-2009.nc")

####### PREPROCESS
### read dimensions 
time <- f1$dim[[1]]$vals
lon <- f1$dim[[2]]$vals
lat <- f1$dim[[3]]$vals

### read entire variables
z500 <- ncvar_get(f1,varid="geopotential_height")
temp <- ncvar_get(f2, varid="temp2")
totprec <- ncvar_get(f3, varid="totprec")

### generate array: lon, lat, year, mon
z500_2 <- array(z500,dim=c(length(lon),length(lat),12,length(time)/12))
temp_2 <- array(temp,dim=c(length(lon),length(lat),12,length(time)/12))
totprec_2 <- array(totprec,dim=c(length(lon),length(lat),12,length(time)/12))

z500_3 <- aperm(z500_2,c(1,2,4,3))
temp_3 <- aperm(temp_2,c(1,2,4,3))
totprec_3 <- aperm(totprec_2,c(1,2,4,3))

### create annual  mean values
z500.yr <- apply(z500_3[,,,1:12],c(1:3),mean)
temp.yr <- apply(temp_3[,,,1:12],c(1:3),mean)
totprec.yr <- apply(totprec_3[,,,1:12],c(1:3),mean)

### create winter (January - March) mean values
z500.jfm <- apply(z500_3[,,,1:3],c(1:3),mean)
temp.jfm <- apply(temp_3[,,,1:3],c(1:3),mean)
totprec.jfm <- apply(totprec_3[,,,1:3],c(1:3),mean)

### create summer (July - September) mean values
z500.jas <- apply(z500_3[,,,7:9],c(1:3),mean)
temp.jas <- apply(temp_3[,,,7:9],c(1:3),mean)
totprec.jas <- apply(totprec_3[,,,7:9],c(1:3),mean)

### do new time axes (yr)
yr <- c(1600:2000)

### apply to predictors
predictors <- read.table("exercise3_predictors.txt",header=TRUE)
selyrs_predictors <- predictors[predictors$YR %in% yr, 2:4]


### response:

# Add a column of ones to represent the intercept.
x <- as.matrix(cbind(1,selyrs_predictors))

#adjust dims
idx     <- which(yr >= 1600 & yr <= 2000)

y <- z500.yr[,,idx]

dim(z500.yr)
dim(y)
dim(x)

xm <- solve(t(x)%*%x)%*%t(x)
dim(xm)


#### regression
b <- apply(y, c(1:2), function(fy,fx) lm(fy~fx)$coef,as.matrix(selyrs_predictors))

b_tsi  <- b[2, , ]
b_co2  <- b[3, , ]
b_volc <- b[4, , ]

lat <- rev(lat) 
b <- b[ , , dim(b)[3]:1]

regmap <- function(field, idx, xm) {
  apply(field[ , , idx], c(1,2), "%*%", t(xm))
}
b_z500_yr   <- regmap(z500.yr,   idx, xm)
b_temp_yr   <- regmap(temp.yr,   idx, xm)
b_temp_jfm  <- regmap(temp.jfm,  idx, xm)
b_temp_jas  <- regmap(temp.jas,  idx, xm)
b_prec_jfm  <- regmap(totprec.jfm, idx, xm)
b_prec_jas  <- regmap(totprec.jas, idx, xm)


### prepare plot
lon2 <- ifelse(lon > 180, lon - 360, lon)
ord  <- order(lon2)          # reorder longitude to rising values
lon2 <- lon2[ord]

plot_map <- function(b, k, main, lon, lat, ord, pal = "RdBu", rev = TRUE) {
  z <- b[k, ord, ]
  zmax   <- max(abs(z), na.rm = TRUE)
  mylevs <- seq(-zmax, zmax, length.out = 22)
  mycol  <- hcl.colors(21, pal, rev = rev)
  
  filled.contour(
    x = lon, y = lat, z = z,
    levels = mylevs, col = mycol,
    plot.title = title(main = main, xlab = "Longitude", ylab = "Latitude"),
    plot.axes  = {
      axis(1); axis(2)                      
      map("world", interior = FALSE, add = TRUE)
    }
  )
}

# a) annual means: TSI (k=2) und CO2 (k=3) on Z500 and temp
par(mfrow = c(2, 2))
#plot_map(b_z500_yr, 2, "Z500 ~ TSI (ANN)",  lon, lat)
#plot_map(b_z500_yr, 3, "Z500 ~ CO2 (ANN)",  lon, lat)
#plot_map(b_temp_yr, 2, "Temp ~ TSI (ANN)",  lon, lat)
#plot_map(b_temp_yr, 3, "Temp ~ CO2 (ANN)",  lon, lat)

# b) Volcanic aerosoles (k=4) on totprec and temp, JFM/JAS
par(mfrow = c(2, 2))
#plot_map(b_prec_jfm, 4, "Prec ~ VOLC (JFM)", lon, lat, ord, "BrBG", rev = FALSE)
#plot_map(b_prec_jas, 4, "Prec ~ VOLC (JAS)", lon, lat, ord, "BrBG", rev = FALSE)
#plot_map(b_temp_jfm, 4, "Temp ~ VOLC (JFM)", lon, lat)
plot_map(b_temp_jas, 4, "Temp ~ VOLC (JAS)", lon, lat)
