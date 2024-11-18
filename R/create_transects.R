#' @export
create_transects <- function(sf,
                             window_size = 3, 
                             width = 500) {
    transect <- list()
    pts <- sf::st_coordinates(sf)
    n_pts <- nrow(pts)
    ref_dir <- c(-1, 0)
    
    for (i in 1:n_pts) {
        start <- max(1, i - floor(window_size / 2))
        end <- min(n_pts, i + floor(window_size / 2))
        window_pts <- pts[start:end, ]
        
        mean_x <- mean(window_pts[, 1])
        mean_y <- mean(window_pts[, 2])
        
        # Calculate covariance matrix components
        covXX <- sum((window_pts[, 1] - mean_x) ^ 2) / nrow(window_pts)
        covXY <- sum((window_pts[, 1] - mean_x) * (window_pts[, 2] - mean_y)) / nrow(window_pts)
        covYY <- sum((window_pts[, 2] - mean_y) ^ 2) / nrow(window_pts)
        
        # Calculate the angle of the heading and perpendicular components
        theta <- 0.5 * atan2(2 * covXY, covXX - covYY)
        perpX <- -sin(theta)
        perpY <- cos(theta)
        
        # Ensure consistent orientation
        dot_product <- perpX * ref_dir[1] + perpY * ref_dir[2]
        if (dot_product < 0) {
            # Flip the perpendicular vector
            perpX <- -perpX
            perpY <- -perpY
        }
        
        # Calculate endpoints of the perpendicular vector
        x0 <- pts[i, 1]
        y0 <- pts[i, 2]
        coords <- rbind(c(x = x0 + width * perpX,
                          y = y0 + width * perpY),
                        c(x = x0 - width * perpX,
                          y = y0 - width * perpY))
        transect[[i]] <- sf::st_linestring(coords)
        
    }
    transect <- do.call(c, transect) %>% 
        sf::st_sfc() %>% 
        sf::st_set_crs(sf::st_crs(sf))
    # st_set_crs(transect) <- crs(points)
    
    return(transect)
}