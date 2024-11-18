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