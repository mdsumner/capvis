
bed0 <- raadtools::topofile("ibcso_bed")
srf0 <- raadtools::topofile("ibcso_is")
basename(bed0)
#"IBCSO_v2_bed.tif"
basename(srf0)
#"IBCSO_v2_ice-surface.tif"


library(anglr)
library(whatarelief)
library(vaster)
library(ximage)
info1 <- vapour::vapour_raster_info(bed0)
info2 <- vapour::vapour_raster_info(srf0)

info1 <- vapour::vapour_raster_info(bed0)
dm <- dev.size("px")
lbed <- elevation(extent = info1$extent, dimension = dm, projection = "EPSG:3031",
                  source = vapour::vapour_vrt(bed0, projection = "EPSG:3031"))
lsrf <- elevation(extent = info1$extent, dimension = dm, projection = "EPSG:3031",
                  source = vapour::vapour_vrt(srf0, projection = "EPSG:3031"))

## trim, we can't work with ximage-ready matrix, needs matrix() orientation
rmat <- t((lsrf - lbed)[nrow(lsrf):1,])
rmat[rmat <=0] <- NA
rowgo <- colSums(is.na(rmat)) == dm[1]
colgo <- rowSums(is.na(rmat)) == dm[2]

## trimmed extent
ex <-  c(sort(x_from_col(info1$extent, dm, range(which(!colgo)))),
              sort(y_from_row(info1$extent, dm, range(which(!rowgo)))))

## now we have dm, and ex that we want
## (the extent of the surface valid data, a smallish dimension to match our requirements)
lbed <- elevation(extent = ex, dimension = dm, projection = "EPSG:3031",
                  source = vapour::vapour_vrt(bed0, projection = "EPSG:3031"))
lsrf <- elevation(extent = ex, dimension = dm, projection = "EPSG:3031",
                  source = vapour::vapour_vrt(srf0, projection = "EPSG:3031"))
par(mar = rep(0, 4))
ximage(lsrf -lbed,  asp = 1, extent = ex)

# # ## for ggplot2
# xyz <- vaster::vaster_long(ex, dimension = dm, data = lsrf - lbed)
# xyz <- xyz[!is.na(xyz[,3]), ]
# xyz <- tibble::as_tibble(xyz)
# ximage(lsrf - lbed, extent = info1$extent, asp = 1)
# library(ggplot2)
# ggplot(xyz) + geom_raster(aes(x, y, fill = z)) + coord_equal() + ggtitle("Bedmap surface - bed")
#
#


#p <- drawPoly()

library(raster)
ndim <- dim(lbed)

line <- rgdal::project(cbind(c(-174, 53), -54, -56), "EPSG:3031")
poly <- rbind(line,
              c(line[2, 1] - 7e6, line[2, 2]),
              c(line[2, 1] - 7e6, line[1, 2]),
              line[1, ])


ximage(lbed, extent = ex, asp = 1)
lines(line)
#
# ti <- tile_index(grout::grout(bed0, blocksize = info1$block))
# windows <- split(t(as.matrix(ti[,2:5])), rep(1:nrow(ti), each = 4))
# ## scan until we get values
# for (i in seq_along(windows)) {
#   d <- vapour::vapour_read_raster_int(bed0, window = windows[[i]])
#
# }

## 5087, 3678 native
#bed <- trim(mask(as_raster(lbed, c(5087, 3678)%/%4), p))
#bed <- as_raster(lbed, ndim %/% deci)
#bed[bed < -2000] <- NA



# bed <- trim(bed)
spl <- spLines(line)
bed <- raster(lbed)
extent(bed) <- extent(ex)
crs(bed) <- "EPSG:3031"
srf <- raster(lsrf)
extent(srf) <- extent(ex)
crs(srf) <- "EPSG:3031"
srf <- trim(mask(srf, spPolygons(poly, crs = prj), inverse = F))

cell <- cellFromLine(bed, spl)[[1]]
zbed <- raster::extract(bed, as.integer(cell))
bad <- is.na(zbed)
cell <- cell[!bad]
zbed <- zbed[!bad]
zsrf <- extract(srf, xyFromCell(bed, cell))
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
prj <- "+proj=stere +lat_ts=-70 +lat_0=-90 +datum=WGS84"
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
