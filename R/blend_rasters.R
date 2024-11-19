#' Blend two rasters together by weighted mean
#' 
#' The function blends two overlapping rasters together by computing a weighted
#' mean between the two rasters. Weights are derived by distance from linestrings
#' associated with each raster.  
#' 
#' @param r1 A SpatRaster
#' @param r2 A second spatRaster 
#' @param l1 A line associated with r1
#' @param l2 A line associated with r2
#' @param blenddist the distance where blending will be made
#' 
#' @returns a blended raster
#' 
#' @export
blend_rasters <- function(r1, r2, 
                          l1, l2, 
                          blenddist) {
    
    #prepare dem
    fulldem <- terra::rast(stars::st_mosaic(stars::st_as_stars(r1),
                                     stars::st_as_stars(r2)))
    # fulldem <- mean(c(resample(r1, fulldem),
    #                   resample(r2, fulldem)),
    #                 na.rm=TRUE)
    resolution <- terra::res(fulldem)
    
    # ensure both lines extend to the blend buffer border
    # This is necessary when the l1 or l2 do not extent to the border
    # of the blend buffer to avoid artefacts.
    l1 <- extend_line(l1, resolution, blenddist)
    l2 <- extend_line(l2, resolution, blenddist)
    
    
    # get buffers from blending distance
    l1b <- sf::st_buffer(l1, blenddist)
    l2b <- sf::st_buffer(l2, blenddist)
    int <- sf::st_intersection(l1b, l2b)
    
    # resample r2 to match r1
    l1dem <- terra::crop(r1, terra::vect(int)) 
    l2dem <- terra::crop(r2, terra::vect(int)) %>% 
        terra::resample(l1dem)
    
    invisible(capture.output(w <- derive_weights(l1, 
                                                 l2,
                                                 l1dem,
                                                 l2dem,
                                                 blenddist)
    ))
    # compute weighted mean, resample, and replace values in the merged fulldem
    r3 <- terra::weighted.mean(c(l1dem, l2dem),
                        w, na.rm=TRUE)
    r3 <- terra::resample(r3, fulldem, method = "near")
    fulldem[!is.na(r3)] <- r3[!is.na(r3)]
    
    return(fulldem)
} 


extend_line <- function(l, resolution, blenddist) {
    crs <- sf::st_crs(l)
    l <- sf::st_coordinates(l)[,1:2]
    bear <- bearing_2d(l[c(nrow(l), 1),])
    move <- move_coords(bear, seq(from = min(resolution), 
                                  blenddist, 
                                  by = min(resolution)*2))
    move[,1] <- move[,1] + l[1,1]
    move[,2] <- move[,2] + l[1,2]
    l <- rbind(move, l)
    l <- sf::st_linestring(l) %>% 
        sf::st_sfc() %>% 
        sf::st_sf() %>% 
        sf::st_set_crs(crs) 
    
    return(l)
}

derive_weights <- function(l1, 
                           l2,
                           l1dem,
                           l2dem,
                           blenddist) {
    d1 <- terra::distance(l1dem, terra::vect(l1)) %>% 
        terra::resample(l1dem, method = "near")
    d2 <- terra::distance(l2dem, terra::vect(l2)) %>% 
        terra::resample(l2dem, method = "near")
    
    # handle edges
    dif <- d2-d1
    maxval <-  unlist(blenddist/terra::global(dif, max, na.rm=TRUE))
    minval <- unlist(blenddist/abs(terra::global(dif, min, na.rm=TRUE)))
    ind <- dif > 0
    dif[ind] <- dif[ind] * maxval
    dif[!ind] <- dif[!ind] * minval
    ind <- dif < 1 & dif >= 0
    dif[ind] <- 1
    ind <- dif > -1 & dif <= 0
    dif[ind] <- -1
    
    lw1 <- abs(dif + blenddist)
    lw2 <- abs(dif - blenddist)
    
    # final weights
    dw1 <- d1
    dw1 <- -dw1 + (terra::global(dw1, min)[[1]])+blenddist
    dw1[dw1 < 0] <- 0
    dw2 <- d2
    dw2 <- -dw2 + (terra::global(dw2, min)[[1]])+blenddist
    dw2[dw2 < 0] <- 0
    
    return(c(dw1, dw2))
}

