#' Straighten a linestring to create a strip
#' 
#' The function creates a straight strip from the linestring provided. Douglas-
#' Peucker simplification algorithm is used to find nodes which are used to 
#' split the linestring in to segments. The segments are then rotated to
#' the angle provided, and stitched together. The result is a straight line
#' with details of the original linestring. The width of the strip (and amount
#' of detail in the linestring) is controlled by the tolerance argument.
#' 
#' @param sf An sf linestring
#' @param tolerance tolerance input for the Douglas-Peucker algorithm.
#' @param angle The direction of of the straightened line. Default = 0 (north).
#' 
#' @returns a list with 1. the straightened line, and 2. the original line with
#' segments defined by the simplification.
#' 
#' @export
straighten_line <- function(sf, 
                            tolerance, 
                            angle = 0) {
    
    
    main_coords <- sf %>% sf::st_coordinates()
    main_coords <- main_coords[!duplicated(main_coords),]
    main <- sf::st_linestring(main_coords[,1:2])
    
    skel <- sf::st_simplify(main, dTolerance = tolerance) #Douglas-Peucker tolerance in map units 
    skel_coords <- sf::st_coordinates(skel)
    
    xi <- main_coords[,1] %in% skel_coords[,1]
    yi <- main_coords[,2] %in% skel_coords[,2]
    inds <- which(xi&yi)
    breakpoints <- main_coords[inds,1:2]
    keep <- which(!duplicated(breakpoints))
    inds <- inds[keep]
    
    for (i in 2:length(inds)) {
        lcoords <- main_coords[inds[i-1]:inds[i],1:2]
        
        start <- sf::st_point(lcoords[1,1:2])
        bearcoords <- lcoords[c(1,nrow(lcoords)),1:2]
        # bear <- bearing(lcoords[1,1:2], lcoords[nrow(lcoords),])
        bear <- bearing_2d(bearcoords)
        rotate <- 360-bear
        # rotate <- ifelse(bear < 0, 180+abs(bear), 180-bear)
        
        
        line <- sf::st_linestring(lcoords[,1:2]) 
        rotated_line <- (line - start) * rotation_matrix(rotate) + start
        rotated_line <- (rotated_line - start) * rotation_matrix(angle) + start
        
        
        if(i == 2) {
            ribbon <- rotated_line
            last_line <- rotated_line
            skeleton <- line
        } else {
            last_end <- sf::st_coordinates(last_line)
            new_start <- sf::st_coordinates(rotated_line)
            
            diff <- sf::st_point(last_end[nrow(last_end),1:2]) -
                sf::st_point(new_start[1,1:2])
            
            rotated_line <- rotated_line + diff
            
            ribbon <- c(ribbon, rotated_line)
            
            last_line <- rotated_line
            
            skeleton <- c(skeleton, line)
        }
    }
    
    ribbon <- ribbon %>% 
        sf::st_sfc() %>% 
        sf::st_sf() %>% 
        sf::st_cast("LINESTRING") %>% 
        sf::st_set_crs(sf::st_crs(sf))
    skeleton <- skeleton %>% 
        sf::st_sfc() %>% 
        sf::st_sf() %>% 
        sf::st_cast("LINESTRING") %>% 
        sf::st_set_crs(sf::st_crs(sf))
    ribbon$ID <- 1:nrow(ribbon)
    skeleton$ID <- 1:nrow(skeleton)
    
    return(list(Ribbon = ribbon,
                Skeleton = skeleton))
}
