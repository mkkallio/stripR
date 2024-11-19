#' Create a strip map with raster data to a specified shape
#' 
#' The function creates a strip of a specified shape, based on a linestring to 
#' be bent to shape. The raster is rotated, shifted and blended in sequence
#' from the start of the line. The rasters and linestrings are placed so that
#' they approximately follow the target shape. 
#' 
#' @param line An sf linestring
#' @param r A SpatRaster 
#' @param tolerance Tolerance for the Douglas-Peucker simplification algorithm.
#' See details.
#' @param buffer_size Buffer size applied to the linestring and used to mask 
#' the raster.
#' @param shape A list of angles and lengths for each segment of the shape. 
#' Potentially obtained using create_shape(). If left NULL, the shape is a
#' straight line pointing upward (angle = 0). See details.
#' @param verbose Whether or not to print progress information.
#' 
#' @returns a list with two elements: the rotated linestring and raster.
#' 
#' @export
strip_raster <- function(line,
                         r,
                         tolerance,
                         buffer_size,
                         shape = NULL, 
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
    
    segment_dems <- list()
    segment_lines <- list()
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
        
        # rib_lines also for adjusting split angle
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
        # crop rasters for each line segment between split points, and rotate
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
                rotate_raster(rotate_angle, start) 
            ln <- (sf::st_geometry(ln) - start) * rotation_matrix(rotate_angle) + start
            r_crop <- terra::mask(r_crop, 
                                  terra::vect(sf::st_buffer(ln, buffer_size)))
            
            
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
    } 
    
    # combine all previous blended rasters
    # TODO: this does not need to be a separate phase, but can be 
    # integrated in the above loop
    # TODO: blending does not currently work for shapes which join one another
    # in the end. This requires work on how the lines and their extensions are
    # handled.
    test <- n_segments > 1
    if(test) {
        if(verbose) {
            npb <- length(segment_dems)
            message("Blending ", npb, " segments...")
            pb <- utils::txtProgressBar(0, npb, style = 3)
        } 
        blend <- segment_dems[[1]]
        for(ii in 2:length(segment_dems)) {
            blend <- blend_rasters(r1 = blend, 
                                   r2 = segment_dems[[ii]],
                                   l1 = segment_lines[[ii-1]], 
                                   l2 = segment_lines[[ii]],
                                   blenddist = buffer_size)
            if(verbose) utils::setTxtProgressBar(pb, ii)
        } 
        if(verbose) close(pb)
    } else {
        blend <- segment_dems[[1]]
    }
   
    
    segments <- do.call(c, segment_lines)
    
    return(list(rotated_line = segments,
                rotated_r = blend))
}



 get_node_indices <- function(i, split_nodes, skel_coords, n_segments) {
     
     test <- n_segments == 1
     if(test) {
         start <- 1
         end <- nrow(skel_coords)
     } else {
         # create lines between splitpoints
         test <- i == 1
         test2 <- i == (length(split_nodes)+1)
         if(test) {
             start <- 1
             end <- split_nodes[i]
         } else if(test2) {
             start <- split_nodes[length(split_nodes)]
             end <- nrow(skel_coords)
         } else {
             start <- split_nodes[i-1]
             end <- split_nodes[i]
         }
     }
     
    return(c(start = start, end = end))
     
 }
 