#' Creates the list of lengths and angles for use in strip_*()-functions
#' 
#' The function takes an input of a linestring, which is simplified based on
#' given tolerance. The lengths and angles between nodes are then recorded.
#' 
#' 
#' @param line a sf linestring object showing the target line.
#' @param tolerance The number of transects to create from the line
#
#' 
#' @returns a list with three elements: 1. the simplified shape, and its 
#' associated 2. lengths, and 3. angles. 
#' 
#' @export
create_shape <- function(line, 
                         tolerance) {
    
    . <- NULL
    
    l <- sf::st_simplify(line, 
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