
#bed0 <- raadtools::topofile("ibcso_bed")
#srf0 <- raadtools::topofile("ibcso_is")
#"/vsicurl/https://download.pangaea.de/dataset/937574/files/IBCSO_v2_bed.nc"

bed0 <- "/rdsi/PUBLIC/raad/data/download.pangaea.de/dataset/937574/files/IBCSO_v2_bed.tif"
srf0 <- "/rdsi/PUBLIC/raad/data/download.pangaea.de/dataset/937574/files/IBCSO_v2_ice-surface.tif"
#img0 <- "/rdsi/PUBLIC/raad/data/download.pangaea.de/dataset/937574/files/IBCSO_v2_ice-surface_RGB.tif"
img0 <- "/vsicurl/https://download.pangaea.de/dataset/937574/files/IBCSO_v2_ice-surface_RGB.tif"

prj <-  "+proj=stere +lat_0=-90 +lat_ts=-65 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84"
bed0 <- vapour::vapour_vrt(bed0, projection = prj)
srf0 <- vapour::vapour_vrt(srf0, projection = prj)
img0 <- vapour::vapour_vrt(img0, projection = prj)
info1 <- vapour::vapour_raster_info(bed0)
info2 <- vapour::vapour_raster_info(srf0)
info3 <- vapour::vapour_raster_info(img0)

library(whatarelief)
library(anglr)

library(vaster)
library(ximage)


#dm <- dev.size("px")
dm <- c(2048, 2048)
ex <- c(-2512500, 2737500, -2212500, 2137500)
isrf <- elevation(extent = ex, dimension = dm, projection = prj,
                  source = srf0)
ibed <- elevation(extent = ex, dimension = dm, projection = prj,
                  source = bed0)
# irgb <- imagery(extent = ex, dimension = dm * 4, projection = prj,
#                 source = img0)


#
#
# mesh <- as.mesh3d(setExtent(raster(isrf), extent(ex)),
#                   image_texture = setExtent(brick(irgb), extent(ex)))
#
#
# plot3d(mesh)



#
# basename(bed0)
# #"IBCSO_v2_bed.tif"
# basename(srf0)
# #"IBCSO_v2_ice-surface.tif"
#
#


#dm <-dev.size("px")
# lbed <- elevation(extent = info1$extent, dimension = dm, projection = info1$projection,
#                   source = vapour::vapour_vrt(bed0, projection = info1$projection))
# lsrf <- elevation(extent = info1$extent, dimension = dm, projection = info1$projection,
#                   source = vapour::vapour_vrt(srf0, projection = info1$projection))
#
# ## trim, we can't work with ximage-ready matrix, needs matrix() orientation
# rmat <- t((lsrf - lbed)[nrow(lsrf):1,])
# rmat[rmat <=0] <- NA
# rowgo <- colSums(is.na(rmat)) == dm[1]
# colgo <- rowSums(is.na(rmat)) == dm[2]
#
# ## trimmed extent
# ex <-  c(sort(x_from_col(info1$extent, dm, range(which(!colgo)))),
#               sort(y_from_row(info1$extent, dm, range(which(!rowgo)))))


## now we have dm, and ex that we want
## (the extent of the surface valid data, a smallish dimension to match our requirements)
# lbed <- elevation(extent = ex, dimension = dm, projection = info1$projection,
#                   source = vapour::vapour_vrt(bed0, projection = info1$projection))
# lsrf <- elevation(extent = ex, dimension = dm, projection = info1$projection,
#                   source = vapour::vapour_vrt(srf0, projection = info1$projection))
par(mar = rep(0, 4))
ximage(isrf -ibed,  asp = 1, extent = ex)

save(info1, ibed, isrf, dm, ex, file = "capvis.Rdata")

# # ## for ggplot2
# xyz <- vaster::vaster_long(ex, dimension = dm, data = lsrf - lbed)
# xyz <- xyz[!is.na(xyz[,3]), ]
# xyz <- tibble::as_tibble(xyz)
# ximage(lsrf - lbed, extent = info1$extent, asp = 1)
# library(ggplot2)
# ggplot(xyz) + geom_raster(aes(x, y, fill = z)) + coord_equal() + ggtitle("Bedmap surface - bed")
#
#
