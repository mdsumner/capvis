library(lazyraster)
library(anglr)
library(rgl)
bed <- lazyraster("C:\\data\\BEDMAP\\bedmap2_bed.tif")
surf <- lazyraster("C:\\data\\BEDMAP\\bedmap2_surface.tif")

ex <- new("Extent", xmin = 1924027.2075697, xmax = 2825048.26982457, 
          ymin = -1352918.47272506, ymax = 504571.101769608)

ex <- new("Extent", xmin = -487936.866774132, xmax = 2866633.85731326, 
          ymin = -2267801.39747617, ymax = 1003598.15163385)

b1 <- as_raster(crop(bed,ex), c(1024, 1024))
s1 <- as_raster(crop(surf,ex), c(1024, 1024))
b1 <- as_raster(bed, c(1024, 1024))
s1 <- as_raster(surf, c(1024, 1024))
b1[b1 < -1500] <- NA
s1[s1 < -1500] <- NA
libary(sp)
pp <- new("SpatialPolygons", polygons = list(new("Polygons", Polygons = list(
  new("Polygon", labpt = c(1800806.92216858, -1169708.16422678
  ), area = 6347715944639.04, hole = FALSE, ringDir = 1L, coords = structure(c(2076507.69502821, 
                                                                               2118093.2825169, 2908219.44480195, 3032976.20726801, 163570.670548626, 
                                                                               343774.882999602, 2076507.69502821, 2076507.69502821, 657051.589228126, 
                                                                               615466.00173944, 407538.064296006, -2364834.43494977, -2448005.60992715, 
                                                                               -1269747.29774769, 657051.589228126, 657051.589228126), .Dim = c(8L, 
                                                                                                                                                2L)))), plotOrder = 1L, labpt = c(1800806.92216858, -1169708.16422678
                                                                                                                                                ), ID = "1", area = 6347715944639.04)), plotOrder = 1L, bbox = structure(c(163570.670548626, 
                                                                                                                                                                                                                           -2448005.60992715, 3032976.20726801, 657051.589228126), .Dim = c(2L, 
                                                                                                                                                                                                                                                                                            2L), .Dimnames = list(c("x", "y"), c("min", "max"))), proj4string = new("CRS", 
                                                                                                                                                                                                                                                                                                                                                                    projargs = NA_character_))


m1 <- as.mesh3d(b1, triangles = FALSE)
m2 <- as.mesh3d(s1, max_triangles = 2e5)

#m2 <- as.mesh3d(raster::mask(s1, pp, inverse = TRUE), max_triangles = 2e5)

prj <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
m1$vb[1:3, ] <- t(quadmesh::llh2xyz(reproj::reproj(t(m1$vb[1:3, ]), source = prj, target = 4326), exag = 20))
m2$vb[1:3, ] <- t(quadmesh::llh2xyz(reproj::reproj(t(m2$vb[1:3, ]), source = prj, target = 4326), exag = 20))

plot3d(m1, color = "white")
wire3d(m2, col = "green", size = 1, alpha = 0.5)
aspect3d("iso")




