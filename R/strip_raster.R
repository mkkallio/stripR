#' @export
strip_raster <- function(line,
                         r,
                         shape, 
                         tolerance,
                         buffer_size,
                         # start_point = NULL,
                         verbose = FALSE) {
    
    L1 <- NULL
    . <- NULL
    
    if(verbose) message("Preparing rotation..")
    
    # create ríbbon from line input
    rib <- straighten_line(line, tolerance, 0)
    skel <- rib$Skeleton
    rib <- rib$Ribbon
    
    # figure out splitpoints
    angles <- shape$angles
    lengths <- shape$lengths
    n_segments <- length(angles)
    
    coords <- sf::st_coordinates(rib)
    miny <- min(coords[,"Y"])
    maxy <- max(coords[,"Y"])
    dy <- (maxy-miny) * cumsum(lengths / sum(lengths))[-n_segments]
    split_at <- miny + dy
    
    # nodes at which splitting takes place
    ind <- sapply(split_at, \(y) {
        which.min(abs(coords[,"Y"] - y))
    })
    
    
    
    # since the skeleton needs to be rotated at more points than the 
    # straight ribbon, each split will have a slightly different 
    # rotation angle. 
    skel_coords <- sf::st_coordinates(skel)
    rib_coords <- sf::st_coordinates(rib)
    
    
    # PROCESS EACH SPLIT
    if(verbose) {
        npb <- length(ind)+1
        message("Processing ", npb, " segments...")
        pb <- utils::txtProgressBar(0, npb, style = 3)
    } 
    
    segment_dems <- list()
    segment_lines <- list()
    for(i in 1:(length(ind)+1)) {
        
        # create lines between splitpoints
        test <- i == 1
        test2 <- i == (length(ind)+1)
        if(test) {
            start <- 1
            end <- ind[i]
        } else if(test2) {
            start <- ind[length(ind)]
            end <- nrow(skel_coords)
        } else {
            start <- ind[i-1]
            end <- ind[i]
        }
        lines <- skel_coords[start:end,] %>% 
            dplyr::as_tibble() %>% 
            dplyr::group_by(L1) %>% 
            dplyr::group_split() %>% 
            lapply(\(x) {sf::st_linestring(as.matrix(x[,1:2])) })  %>%
            do.call(c, .) %>% 
            sf::st_sfc() %>% 
            sf::st_set_crs(sf::st_crs(skel)) %>% 
            sf::st_cast("LINESTRING") 
        lines_buf <- sf::st_buffer(lines, buffer_size)
        
        target_angle <- (angles[i] + 360) %% 360
        
        # split angle
        rib_line <- rib_coords[start:end,] %>% 
            dplyr::as_tibble() %>% 
            dplyr::group_by(L1) %>% 
            dplyr::group_split() %>% 
            lapply(\(x) {sf::st_linestring(as.matrix(x[,1:2])) }) %>% 
            do.call(c, .) %>% 
            sf::st_sfc() %>% 
            sf::st_set_crs(sf::st_crs(skel)) %>% 
            sf::st_cast("LINESTRING") 
        
        # ----------------------
        # crop rasters for each line segment betweem splitpoints, and rotate
        rr <- list()
        rln <- list()
        for(ii in seq_along(lines)) {
            ln <- lines[ii] 
            r_crop <- terra::crop(r, terra::ext(sf::st_sf(lines_buf[ii]))) %>% 
                terra::mask(terra::vect(lines_buf[ii]))
            start <- lwgeom::st_startpoint(ln) 
            end <- lwgeom::st_endpoint(ln) 
            
            # figure out the rotation angle
            ln_coords <- sf::st_coordinates(c(start, end))
            bear <- bearing_2d(ln_coords)
            rotate_angle <- (target_angle + 360 - bear) %% 360
            
            #update rotate_angle based on the ribbon
            rb_start <- lwgeom::st_startpoint(rib_line[ii]) 
            rb_end <- lwgeom::st_endpoint(rib_line[ii]) 
            rb_coords <- sf::st_coordinates(c(rb_start, rb_end))
            bear_rb <- bearing_2d(rb_coords)
            rotate_angle <- (rotate_angle + bear_rb) %% 360
            
            # rotate raster and the line
            r_crop <- r_crop %>% 
                rotate_raster(start, rotate_angle, trim = TRUE) 
            ln <- (sf::st_geometry(ln) - start) * rotation(rotate_angle) + start
            r_crop <- terra::mask(r_crop, terra::vect(sf::st_buffer(ln, buffer_size)))
            
            
            # move rotated features so that they connect with the previous
            # split
            test <- ii > 1 || i > 1
            if(test) {
                move <- last_end - start
                ln <- ln + move
                
                move <- sf::st_coordinates(move)
                r_crop <- terra::shift(r_crop, dx = move[1], dy = move[2])
            }
            
            # record to the list
            ln <- sf::st_set_crs(ln, sf::st_crs(lines))
            last_end <- lwgeom::st_endpoint(ln)
            rr[[ii]] <- r_crop
            rln[[ii]] <- ln
        }
        
        
        
        # BLEND ROTATED AND SHIFTED RASTERS
        n <- length(rr)
        test <- n > 1
        if(test) {
            for(ii in 2:length(rr)) {
                blended <- blend_rasters(rr[[ii-1]], rr[[ii]],
                                         rln[[ii-1]], rln[[ii]],
                                         blenddist = buffer_size)
                
                if(ii == n) {
                    segment_dems[[i]] <- blended
                } else {
                    rr[[ii]] <- blended
                }
            } 
        } else {
            segment_dems[[i]] <- rr[[1]]
        }
        
        segment_lines[[i]] <- do.call(c, rln)
        
        if(verbose) utils::setTxtProgressBar(pb, i)
    }
    if(verbose) {
        close(pb)
        npb <- length(segment_dems)
        message("Blending ", npb, " segments...")
        pb <- utils::txtProgressBar(0, npb, style = 3)
    } 

    # combine all previous blended rasters
    # TODO: this does not need to be a separate phase, but can be 
    # integrated in the above loop
    blend <- segment_dems[[1]]
    for(ii in 2:length(segment_dems)) {
        blend <- blend_rasters(r1 = blend, 
                               r2 = segment_dems[[ii]],
                               l1 = segment_lines[[ii-1]], 
                               l2 = segment_lines[[ii]],
                               blenddist = buffer_size)
        if(verbose) utils::setTxtProgressBar(pb, ii)
        message(ii)
    } 
    if(verbose) close(pb)
    
    segments <- do.call(c, segment_lines)
    
    return(list(rotated_line = segments,
                rotated_r = blend))
}