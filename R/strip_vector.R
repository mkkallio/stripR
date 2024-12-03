#' Create a strip map with vector data to a specified shape. EXPERIMENTAL
#' 
#' The function creates a strip of a specified shape, based on a linestring to 
#' be bent to shape. The vectors are rotated, shifted and blended in sequence
#' from the start of the line. The vectors and linestrings are placed so that
#' they approximately follow the target shape. 
#' 
#' @param line An sf linestring
#' @param v sf object with the vectors to be strip'd.
#' @param tolerance Tolerance for the Douglas-Peucker simplification algorithm.
#' See details.
#' @param buffer_size Buffer size applied to the linestring and used to mask 
#' the vectors.
#' @param shape A list of angles and lengths for each segment of the shape. 
#' Potentially obtained using create_shape(). If left NULL, the shape is a
#' straight line pointing upward (angle = 0). See details.
#' @param verbose Whether or not to print progress information.
#' 
#' @returns a list with two elements: the rotated linestring and vectors.
#' 
#' @export
strip_vector <-  function(line,
                          v,
                          tolerance,
                          buffer_size,
                          shape = NULL, 
                          remove_overlap = TRUE,
                          split_features = FALSE,
                          verbose = FALSE) {
    
    L1 <- NULL
    . <- NULL
    
    if(verbose) message("Preparing rotation..")
    
    # create ríbbon from line input
    rib <- straighten_line(line, tolerance, 0)
    skel <- rib$Skeleton
    rib <- rib$Ribbon
    
    # figure out split points
    coords <- sf::st_coordinates(rib)
    test <- is.null(shape)
    if(test) {
        angles <- 0
        lengths <- 1
        n_segments <- 1
        
        split_nodes <- nrow(coords)
        inds <- 1
    } else {
        test <- utils::hasName(shape, c("angles", "lengths"))
        if(any(!test)) stop("Shape not valid. Use create_shape()-function first.")
        angles <- shape$angles
        lengths <- shape$lengths
        n_segments <- length(angles)
        
        miny <- min(coords[,"Y"])
        maxy <- max(coords[,"Y"])
        dy <- (maxy-miny) * cumsum(lengths / sum(lengths))[-n_segments]
        split_at <- miny + dy
        
        # nodes at which splitting takes place
        split_nodes <- sapply(split_at, \(y) {
            which.min(abs(coords[,"Y"] - y))
        })
        inds <- 1:(length(split_nodes) +1)
    }
    
    
    
    # since the skeleton needs to be rotated at more points than the 
    # straight ribbon, each split will have a slightly different 
    # rotation angle. 
    skel_coords <- sf::st_coordinates(skel)
    rib_coords <- sf::st_coordinates(rib)
    
    
    # PROCESS EACH SPLIT
    if(verbose) {
        npb <- n_segments
        message("Processing ", npb, " segments...")
        pb <- utils::txtProgressBar(0, npb, style = 3)
    } 
    
    # segment_vect <- list()
    # segment_lines <- list()
    rv <- list()
    rln <- list()
    for(i in inds) {
        
        node_ind <- get_node_indices(i, split_nodes, skel_coords, n_segments)
        start <- node_ind[1]
        end <- node_ind[2]
        
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
            v_crop_moved <-  (sf::st_geometry(v_crop) - start) * rotation_matrix(rotate_angle) + start
            ln <- (sf::st_geometry(ln) - start) * rotation_matrix(rotate_angle) + start
            
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
            ln <- sf::st_sf(ln,
                            major_rotate = i,
                            minor_rotate = ii,
                            rotation_angle = rotate_angle)
            
            rv <- append(rv, list(v_crop))
            rln <- append(rln, list(ln))
        }
        
        # segment_lines[[i]] <- do.call(rbind, rln)
        # segment_vect[[i]] <- do.call(rbind, rv)
        
        if(verbose) utils::setTxtProgressBar(pb, i)
    }
    if(verbose) close(pb)
    
    if(remove_overlap) {
        rv <- remove_overlap(rln, rv, buffer_size, split_features)
    }
    
    
    segments <- do.call(rbind, rln)
    vectors <- do.call(rbind, rv)
    
    return(list(rotated_line = segments,
                rotated_v = vectors))
}


remove_overlap <- function(rln, rv, buffer_size, split_features) {
    
    crs <- st_crs(rv[[1]])
    
    # REMOVE OVERLAPS FROM VECTORS
    for(i in 2:length(rv)) {
        
        # create transect and buffer zones
        l1 <- rln[[i-1]]
        l2 <- rln[[i]]
        
        l1c <- sf::st_coordinates(l1)
        l2c <- sf::st_coordinates(l2)
        
        nl1c <- nrow(l1c)
        transect <- rbind(l1c[(nl1c-1):nl1c,1:2],
                          l2c[2,1:2]) %>%
            sf::st_linestring() %>% 
            sf::st_sfc() %>% 
            create_transects(width = buffer_size*2)
        transect <- transect[2]
        
        t_left <- sf::st_buffer(transect, buffer_size, singleSide = TRUE) 
        t_right <- sf::st_buffer(transect, -buffer_size, singleSide = TRUE) 
        
        v1 <- rv[[i-1]] %>% sf::st_set_crs(NA)
        v2 <- rv[[i]] %>% sf::st_set_crs(NA)
        
        # handle different geometry types
        if(split_features) {
            v1 <- st_intersection(t_left, v1) %>% 
                st_union() %>% 
                st_buffer(0.1) %>% 
                st_difference(v1, .) %>% 
                sf::st_set_crs(crs)
            
            v2 <- st_intersection(t_right, v2) %>% 
                st_union() %>% 
                st_buffer(0.1) %>% 
                st_difference(v2, .) %>% 
                sf::st_set_crs(crs)

        } else {
           predicate <- sf::st_intersects
           
           ind <- v1 %>% 
               predicate(t_left, ., sparse = FALSE) %>% 
               as.vector()    
           v1 <- dplyr::filter(v1, !ind) %>% 
               sf::st_set_crs(crs)
           
           ind <- v2 %>% 
               predicate(t_right, ., sparse = FALSE) %>% 
               as.vector()    
           v2 <- dplyr::filter(v2, !ind) %>% 
               sf::st_set_crs(crs)
        }
       
      
       
        rv[[i-1]] <- v1
        rv[[i]] <- v2
    }
    
    return(rv)
}
