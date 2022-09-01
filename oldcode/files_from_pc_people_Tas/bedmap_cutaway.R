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

denser <- 4
n_cells <- length(cellFromLine(bed, spLines(line))[[1]]) * denser
idx <- 1:n_cells
xy <- cbind(approxfun(range(idx), line[,1])(idx), 
            approxfun(range(idx), line[,2])(idx))
xyz <- cbind(xy, extract(bed, xy, method = "bilinear"))
xyz2 <- cbind(xy, extract(srf, xy, method = "bilinear"))

bad <- is.na(xyz2[,3])
xyz <- xyz[!bad, ]
xyz2 <- xyz2[!bad, ]
l <- rbind(xyz, xyz2[nrow(xyz2):1, ], xyz[1, ])

## http://natverse.org/news/2019/02/nat-1-8-13/
uv <- Morpho::tangentPlane(l)
DotProduct=function(a,b) sum(a*b)
# find coordinates w.r.t. to our two basis vectors
xy_tri <- cbind(u=apply(l, 1, DotProduct, uv[[1]]),
                v=apply(l, 1, DotProduct, uv[[2]]))
plot(xy_tri)
p2s <- function(x) cbind(.vx0 = utils::head(x, -1L), .vx1 = utils::tail(x, -1))

#S <- p2s(seq_len(nrow(xy_tri)-1))
#S <- rbind(S, c(S[length(S)], 1L))
#ps <- RTriangle::pslg(P = xy_tri[-nrow(xy_tri), ], S = S)
#tri <- RTriangle::triangulate(ps, Y= TRUE, j = TRUE, S = 0)

tri <- decido::earcut(xy_tri)







# # poly_xyz <- rbind(bed_xyz, srf_xyz)
# path2quad <- function(x) {
#   cbind(head(x, -1L), tail(x, -1))[,c(1, 2, 2, 1)]
# }
# mkquad <- function(xx) {
#  idx <- path2quad(seq_len(nrow(xx)))
#  xy <- xyzz[as.vector(t(idx)), 1:2]
#  z <- as.vector(rbind(matrix(xx[as.vector(t(idx[,1:2])), 3], 2), matrix(xx[as.vector(t(idx[,1:2])), 4], 2)[2:1, ]))
#  cbind(xy, z)
# }
# quds <- mkquad(xyzz)
# quds[,1:2] <- reproj::reproj(quds[,1:2], source = prj, target = 4326)[,1:2]
# 
# 
# quds <- quadmesh::llh2xyz(quds, exag = exag)
# quads3d(quds, color = "white")

library(rgl)

msrf <- as.mesh3d(srf, triangles = FALSE)
mbed <- as.mesh3d(bed, triangles = FALSE) #TRUE, max_triangles = 2e6)

## globurate
exag <- 10
msrf$vb[1:3, ] <- t(quadmesh::llh2xyz(reproj::reproj(t(msrf$vb[1:3, ]), source = prj, target = 4326), exag = exag))
mbed$vb[1:3, ] <- t(quadmesh::llh2xyz(reproj::reproj(t(mbed$vb[1:3, ]), source = prj, target = 4326), exag = exag))

mbed$meshColor <- "vertices"
mbed$material$color <- scales::alpha(ifelse(mbed$vb[3, ] > 0, "grey", "lightblue4"))
plot3d(mbed, xlab = "", ylab = "", zlab = "")
wire3d(msrf, color = "white")

clear3d(type = "bboxdec")
bg3d("lightskyblue4")

aspect3d("iso")
ll <- l
ll[,1:2] <- reproj::reproj(ll[,1:2], target = 4326, source = prj)[,1:2]
triangles3d(quadmesh::llh2xyz(ll[tri, ], exag = exag), color = "white")

g <- graticule::graticule(seq(-180, 180, by = 10), seq(-85, -60, by = 5), xlim = c(-180, 180), tiles = TRUE)
plot3d(globe(silicate::SC0(g)), add = TRUE, color = "black", gproj = "+proj=geocent +a=4308137")
cl <- silicate::SC0(SOmap::SOmap_data$continent[2, ])
cl$vertex$z_ <- 500
plot3d(globe(cl), color = "lightblue", add = TRUE, lwd = 1)

play3d(spin3d(rpm = 2))
#plot3d(globe(silicate::TRI0(g)), add = TRUE)
