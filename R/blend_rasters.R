#' @export
blend_rasters <- function(r1, r2, 
                          l1, l2, 
                          blenddist) {
    
    #prepare dem
    # fulldem <- sprc(r1, r2) %>% merge()
    fulldem <- terra::rast(stars::st_mosaic(stars::st_as_stars(r1),
                                     stars::st_as_stars(r2)))
    # fulldem <- mean(c(resample(r1, fulldem),
    #                   resample(r2, fulldem)),
    #                 na.rm=TRUE)
    resolution <- terra::res(fulldem)
    
    # ensure both lines extend to the blend buffer border
    # This is necessary when the l1 or l2 do not extent to the border
    # of the blend buffer to avoid artefacts.
    crs <- sf::st_crs(l1)
    l1 <- sf::st_coordinates(l1)[,1:2]
    bear <- bearing_2d(l1[c(nrow(l1), 1),])
    move <- move_coords(bear, seq(from = min(resolution), 
                                  blenddist, 
                                  by = min(resolution)*2))
    move[,1] <- move[,1] + l1[1,1]
    move[,2] <- move[,2] + l1[1,2]
    l1 <- rbind(move, l1)
    l1 <- sf::st_linestring(l1) %>% 
        sf::st_sfc() %>% 
        sf::st_sf() %>% 
        sf::st_set_crs(crs) 
    
    crs <- sf::st_crs(l2)
    l2 <- sf::st_coordinates(l2)[,1:2]
    bear <- bearing_2d(l2[c(1, nrow(l2)),])
    move <- move_coords(bear, seq(from = min(resolution), 
                                  blenddist, 
                                  by = min(resolution)*2))
    move[,1] <- move[,1] + l2[nrow(l2),1]
    move[,2] <- move[,2] + l2[nrow(l2),2]
    l2 <- rbind(l2, move)
    l2 <- sf::st_linestring(l2) %>% 
        sf::st_sfc %>% 
        sf::st_sf() %>% 
        sf::st_set_crs(crs)
    
    
    # l1 <- rln[[1]]
    # l2 <- rln[[2]]
    l1b <- sf::st_buffer(l1, blenddist)
    l2b <- sf::st_buffer(l2, blenddist)
    int <- sf::st_intersection(l1b, l2b) #%>% 
    # st_union()
    
    # l1 <- st_intersection(l1, st_buffer(int, 25000)) # change as needed
    # int <- st_buffer(int, 10000)
    
    l1dem <- terra::crop(r1, terra::vect(int)) #%>% 
    # mask(vect(int))
    l2dem <- terra::crop(r2, terra::vect(int)) %>% 
        # mask(vect(int)) %>% 
        terra::resample(l1dem)
    
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
    
    dw1 <- d1
    dw1 <- -dw1 + (terra::global(dw1, min)[[1]])+blenddist
    dw1[dw1 < 0] <- 0
    
    
    dw2 <- d2
    dw2 <- -dw2 + (terra::global(dw2, min)[[1]])+blenddist
    dw2[dw2 < 0] <- 0
    
    
    r3 <- terra::weighted.mean(c(l1dem, l2dem),
                        c(dw1, dw2), na.rm=TRUE)
    r3 <- terra::resample(r3, fulldem, method = "near")
    # r2 <- resample(r2, fulldem, method = "near")
    
    fulldem[!is.na(r3)] <- r3[!is.na(r3)]
    # fulldem[is.na(fulldem)] <- r2[is.na(fulldem)]
    
    return(fulldem)
} 