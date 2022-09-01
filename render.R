
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
#bed0 <- raadtools::topofile("ibcso_bed")
#srf0 <- raadtools::topofile("ibcso_is")
#"/vsicurl/https://download.pangaea.de/dataset/937574/files/IBCSO_v2_bed.nc"

#bed0 <- "/rdsi/PUBLIC/raad/data/download.pangaea.de/dataset/937574/files/IBCSO_v2_bed.tif"
#srf0 <- "/rdsi/PUBLIC/raad/data/download.pangaea.de/dataset/937574/files/IBCSO_v2_ice-surface.tif"
#img0 <- "/rdsi/PUBLIC/raad/data/download.pangaea.de/dataset/937574/files/IBCSO_v2_ice-surface_RGB.tif"
img0 <- "/vsicurl/https://download.pangaea.de/dataset/937574/files/IBCSO_v2_ice-surface_RGB.tif"

bed0 <- "IBCSO_v2_bed.tif"
srf0 <- "IBCSO_v2_ice-surface.tif"
prj <-  "+proj=stere +lat_0=-90 +lat_ts=-65 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84"


g <- graticule::graticule(seq(-180, 180, by = 10), seq(-85, -60, by = 5), xlim = c(-180, 180), tiles = TRUE)

bed0 <- vapour::vapour_vrt(bed0, projection = prj)
srf0 <- vapour::vapour_vrt(srf0, projection = prj)
#img0 <- vapour::vapour_vrt(img0, projection = prj)
info1 <- vapour::vapour_raster_info(bed0)
info2 <- vapour::vapour_raster_info(srf0)
#info3 <- vapour::vapour_raster_info(img0)

library(raster)
library(whatarelief)
library(ximage)
library(anglr)


#dm <- dev.size("px")
dm <- c(2048, 2048)
ex <- c(-2512500, 2737500, -2212500, 2137500)
ex <- info1$extent
isrf <- elevation(extent = ex, dimension = dm, projection = prj,
                  source = srf0)
ibed <- elevation(extent = ex, dimension = dm, projection = prj,
                  source = bed0)
# irgb <- imagery(extent = ex, dimension = dm * 4, projection = prj,
#                 source = img0)

#load("capvis.Rdata")

#p <- drawPoly()


line <- rgdal::project(cbind(c(-174, 53), -54, -56), info1$projection)
poly <- rbind(line,
              c(line[2, 1] - 7e6, line[2, 2]),
              c(line[2, 1] - 7e6, line[1, 2]),
              line[1, ])

spl <- spLines(line)
bed <- raster(ibed)
extent(bed) <- extent(ex)
crs(bed) <- info1$projection
srf <- raster(isrf)
extent(srf) <- extent(ex)
crs(srf) <- info1$projection
icemask0 <- terra::as.polygons(terra::rast(srf >0))[2, ]
icemask <- terra::intersect(icemask0, terra::vect(poly, type = "polygon", crs = info1$projection))
icemask <- sf::as_Spatial(sf::st_as_sf(icemask))
#icemask0 <- sf::as_Spatial(sf::st_as_sf(icemask))
srf2 <- srf
srf <- mask(srf, icemask, inverse = F)
#srf <- trim(mask(srf, spPolygons(poly, crs = info1$projection), inverse = F))

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
quds <- mkquad(xyzz)
quds[,1:2] <- sf::sf_project(pts = quds[,1:2], from = info1$projection, to = "OGC:CRS84")[,1:2]
exag <- 10

quds <- quadmesh::llh2xyz(quds, exag = exag)

msrf2 <- as.mesh3d(srf2, triangles = F)
msrf <- as.mesh3d(srf, triangles = F)
mbed <- as.mesh3d(bed, triangles = F) #TRUE, max_triangles = 2e6)

mbed$meshColor <- "vertices"
msrf$meshColor <- "vertices"
mbed$material$color <- rep(scales::alpha("#243D57", 1), ncol(mbed$vb))
mbed$material$color[na.omit(which(msrf2$vb[3, ] > 0))] <- scales::alpha("#654933", 1)
#plot3d(mbed)


bedllmeshpts <- sf::sf_project(pts = t(mbed$vb[1:3, ]), from = info1$projstring, to = "OGC:CRS84")
bedllmeshpts[bedllmeshpts[,2] > -60, ] <- NA
srfllmeshpts <- sf::sf_project(pts = t(msrf$vb[1:3, ]), from = info1$projstring, to = "OGC:CRS84")
srfllmeshpts[srfllmeshpts[,2] > -60, ] <- NA

msrf$vb[1:3, ] <- t(quadmesh::llh2xyz(srfllmeshpts, exag = exag))
mbed$vb[1:3, ] <- t(quadmesh::llh2xyz(bedllmeshpts, exag = exag))


plot3d(mbed, xlab = "", ylab = "", zlab = "", specular = "black")
shade3d(msrf, color = "white")


cpt <- cbind(122.5918498657864, -75.28014648303909)
cz <- extract(srf, rgdal::project(cpt, projection(srf)))
core <- quadmesh::llh2xyz(cbind(cpt[c(1, 1), ], c(-1000, cz * 3)),
                          exag = exag)

# l <- quadmesh::llh2xyz(cbind(c(-180, 180), c(-90, -40), 0),
#                        exag = exag)
#
# lines3d(l)


library(rgl)
lines3d(core, lwd = 3, col = "red")

quads3d(quds, color = "white")
clear3d(type = "bboxdec")
bg3d("lightskyblue4")

aspect3d("iso")
bg3d("grey60")
g <- silicate::SC0(g)
g$vertex <- setNames(tibble::as_tibble(quadmesh::llh2xyz(cbind(as.matrix(g$vertex), z_ = 0))),
                     c("x_", "y_", "z_"))
g$meta$proj <- "+proj=geocent +datum=WGS84"
plot3d(g, add = TRUE, color = "black", gproj = "+proj=geocent +a=4308137")
# cl <- silicate::SC0(SOmap::SOmap_data$continent[2, ])
# cl$vertex$z_ <- 500
# plot3d(globe(cl), color = "black", add = TRUE, lwd = 2)
#ice <- raadtools::readice()
#play3d(spin3d(rpm = 2))
#plot3d(globe(silicate::TRI0(g)), add = TRUE)
#rgl::snapshot3d("capvis3.png", webshot = FALSE)
#save(mbed, msrf, core, quds, g, file = "capvis_render.Rdata")
