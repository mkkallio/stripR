#' @export
strip_vector <-  function(line,
                          v,
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
    
    segment_vect <- list()
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
        rv <- list()
        rln <- list()
        for(ii in seq_along(lines)) {
            ln <- lines[ii] 
            v_crop <- sf::st_intersection(v, lines_buf[ii]) 
            
            start <- lwgeom::st_startpoint(ln) 
            end <- lwgeom::st_endpoint(ln) 
            
            # figure out the rotation angle
            ln_coords <- sf::st_coordinates(c(start, end))
            bear <- bearing_2d(ln_coords)
            rotate_angle <- (target_angle + 360 - bear) %% 360
            
            #update rotate_angle
            rb_start <- lwgeom::st_startpoint(rib_line[ii]) 
            rb_end <- lwgeom::st_endpoint(rib_line[ii]) 
            rb_coords <- sf::st_coordinates(c(rb_start, rb_end))
            bear_rb <- bearing_2d(rb_coords)
            rotate_angle <- (rotate_angle + bear_rb) %% 360
            
            # rotate raster and the line
            v_crop_moved <-  (sf::st_geometry(v_crop) - start) * rotation(rotate_angle) + start
            ln <- (sf::st_geometry(ln) - start) * rotation(rotate_angle) + start
            
            # move rotated features so that they connect with the previous
            # split
            test <- ii > 1 || i > 1
            if(test) {
                move <- last_end - start
                ln <- ln + move
                v_crop_moved <- v_crop_moved + move
                
            }
            
            # record to the list
            ln <- sf::st_set_crs(ln, sf::st_crs(lines))
            v_crop_moved <- sf::st_set_crs(v_crop_moved, sf::st_crs(v))
            last_end <- lwgeom::st_endpoint(ln)
            
            v_crop <- sf::st_set_geometry(v_crop, v_crop_moved) %>% 
                dplyr::mutate(major_rotate = i,
                       minor_rotate = ii,
                       rotation_angle = rotate_angle)
            rv[[ii]] <- v_crop
            rln[[ii]] <- ln
        }
        
        segment_lines[[i]] <- do.call(c, rln)
        segment_vect[[i]] <- do.call(rbind, rv)
        
        if(verbose) utils::setTxtProgressBar(pb, i)
    }
    if(verbose) close(pb)
    
    
    segments <- do.call(c, segment_lines)
    vectors <- do.call(rbind, segment_vect)
    
    return(list(rotated_line = segments,
                rotated_v = vectors))
}