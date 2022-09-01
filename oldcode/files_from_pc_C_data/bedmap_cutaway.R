bed0 <- "\\\\aad.gov.au\\files\\AADC\\Scientific_Data\\Data\\gridded_new\\data\\secure.antarctica.ac.uk\\data\\bedmap2\\bedmap2_tiff\\bedmap2_bed.tif"
srf0 <- "\\\\aad.gov.au\\files\\AADC\\Scientific_Data\\Data\\gridded_new\\data\\secure.antarctica.ac.uk\\data\\bedmap2\\bedmap2_tiff\\bedmap2_surface.tif"

library(anglr)
library(lazyraster)
library(raster)
lbed <- lazyraster(bed0)
lsrf <- lazyraster(srf0)

ex <- extent(lbed)

#p <- drawPoly()

ndim <- dim(raster::raster(bed0))[1:2]
prj <- projection(raster::raster(bed0))
deci <- 2
line <- reproj::reproj(cbind(c(-174, 53), -54, -56), source = 4326, target = prj)[,1:2]
poly <- rbind(line, 
      c(line[2, 1] - 7e6, line[2, 2]), 
      c(line[2, 1] - 7e6, line[1, 2]), 
      line[1, ])


## 5087, 3678 native
#bed <- trim(mask(as_raster(lbed, c(5087, 3678)%/%4), p))
bed <- as_raster(lbed, ndim %/% deci)
#bed[bed < -2000] <- NA

bed <- trim(bed)
srf <- trim(mask(as_raster(lsrf, ndim %/% deci), spPolygons(poly, crs = prj), inverse = F))
spl <- spLines(line)
cell <- cellFromLine(bed, spl)[[1]]
zbed <- raster::extract(bed, as.integer(cell))
bad <- is.na(zbed)
cell <- cell[!bad]
zbed <- zbed[!bad]
zsrf <- extract(as_raster(lsrf, ndim %/% deci), xyFromCell(bed, cell))
bad <- is.na(zsrf)
cell <- cell[!bad]
zbed <- zbed[!bad]
zsrf <- zsrf[!bad]

xyzz <- cbind(xFromCell(bed, cell), yFromCell(bed, cell), zbed, zsrf)


# poly_xyz <- rbind(bed_xyz, srf_xyz)
path2quad <- function(x) {
  cbind(head(x, -1L), tail(x, -1))[,c(1, 2, 2, 1)]
}
mkquad <- function(xx) {
 idx <- path2quad(seq_len(nrow(xx)))
 xy <- xyzz[as.vector(t(idx)), 1:2]
 z <- as.vector(rbind(matrix(xx[as.vector(t(idx[,1:2])), 3], 2), matrix(xx[as.vector(t(idx[,1:2])), 4], 2)[2:1, ]))
 cbind(xy, z)
}

msrf <- as.mesh3d(srf, triangles = FALSE)
mbed <- as.mesh3d(bed, triangles = FALSE) #TRUE, max_triangles = 2e6)

## globurate
exag <- 10
msrf$vb[1:3, ] <- t(quadmesh::llh2xyz(reproj::reproj(t(msrf$vb[1:3, ]), source = prj, target = 4326), exag = exag))
mbed$vb[1:3, ] <- t(quadmesh::llh2xyz(reproj::reproj(t(mbed$vb[1:3, ]), source = prj, target = 4326), exag = exag))

plot3d(mbed, color = "grey", xlab = "", ylab = "", zlab = "")
wire3d(msrf, color = "white")

quds <- mkquad(xyzz)
quds[,1:2] <- reproj::reproj(quds[,1:2], source = prj, target = 4326)[,1:2]


quds <- quadmesh::llh2xyz(quds, exag = exag)
library(rgl)
quads3d(quds, color = "white")
clear3d(type = "bboxdec")
bg3d("lightskyblue4")

aspect3d("iso")

g <- graticule::graticule(seq(-180, 180, by = 10), seq(-85, -60, by = 5), xlim = c(-180, 180), tiles = TRUE)
plot3d(globe(silicate::SC0(g)), add = TRUE, color = "black", gproj = "+proj=geocent +a=4308137")
cl <- silicate::SC0(SOmap::SOmap_data$continent[2, ])
cl$vertex$z_ <- 500
plot3d(globe(cl), color = "black", add = TRUE, lwd = 2)
ice <- raadtools::readice()
play3d(spin3d(rpm = 2))
#plot3d(globe(silicate::TRI0(g)), add = TRUE)
