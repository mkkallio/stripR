#' Create a strip map with raster data using transects. EXPERIMENTAL
#' 
#' The function creates a strip of the specified shape with raster data 
#' using transects derived from the linestring. Raster data is then sampled
#' from the raster with the linestring. The values are then translated to
#' transects computed from the target shape, and gaps are filled with 
#' using IDW interpolation. EXPERIMENTAL AND SUBJECT TO CHANGE!
#' 
#' 
#' @param line An sf linestring
#' @param r sf object with the vectors to be strip'd.
#' @param shape a sf linestring object showing the target shape.
#' @param n_transects The number of transects to create from the line
#' @param width Width of the transect
#' @param radius radius for inverse distance weighting-based interpolation to
#' fill gaps between transects.
#' @param ... Further arguments passed to terra::interpIDW()
#' 
#' @returns a list with two elements: the rotated linestring and vectors.
#' 
#' @export
strip_transect <- function(line, 
                           r, 
                           shape,
                           n_transects, 
                           width,
                           radius,
                           ...) {
    
    
    line_transects <- line %>%
        sf::st_line_sample(n_transects) %>% 
        create_transects(window_size = 3, width = width) %>% 
        sf::st_cast("LINESTRING")
    
    
    shape_transects <- shape %>%
        sf::st_line_sample(n_transects) %>% 
        create_transects(window_size = 3, width = width) %>% 
        sf::st_cast("LINESTRING")
    
    
    r_values <- terra::extract(r, terra::vect(line_transects))
    
    points <- list()
    for(i in unique(r_values$ID)) {
        
        v <- r_values[r_values[,"ID"] == i,]
        n <- nrow(v)
        
        pts <- shape_transects[i] %>% 
            sf::st_line_sample(n) %>% 
            sf::st_cast("POINT") %>% 
            sf::st_sf(v = v[,2])
        points[[i]] <- pts
    }
    points <- do.call(rbind, points) %>% 
        sf::st_set_crs(sf::st_crs(line))
    
    template <- terra::rast(extent = terra::ext(points), 
                            resolution = terra::res(r))
    
    r2 <- terra::interpIDW(template, 
                           terra::vect(points), 
                           field = "v", 
                           radius = radius, 
                           ...)
    
    return(r2)
}