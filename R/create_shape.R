#' @export
create_shape <- function(shape, 
                         tolerance) {
    
    . <- NULL
    
    l <- sf::st_simplify(shape, 
                         preserveTopology = TRUE,
                         dTolerance = tolerance) %>% 
        sf::st_geometry() %>% 
        sf::st_sfc() 
    
    p <- sf::st_cast(l, "POINT")
    
    l <- lwgeom::st_split(l, p) %>% 
        sf::st_collection_extract("LINESTRING")
    
    # casting orders features by latitude. need to preserve the correct order.
    order <- l %>% 
        lwgeom::st_startpoint() %>% 
        sf::st_equals(p) %>% 
        unlist() %>% 
        match(1:length(l), .)
    l <- l[order]
    
    lengths <- units::drop_units(sf::st_length(l))
    angles <- sapply(seq_along(l), \(i) {
        l[i] |>
            sf::st_coordinates() |>
            bearing_2d()
    })
    n <- length(l)
    names(lengths) <- 1:n
    names(angles) <- 1:n
    
    return(list(shape = sf::st_sf(l, id = 1:n),
                lengths = lengths,
                angles = angles))
}