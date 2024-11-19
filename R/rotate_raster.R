#' Rotate a raster around a specified point
#' 
#' The function rotates a raster around a specified point.Inspired by Damien
#' Caillaud's answer in Stackoverflow: https://gis.stackexchange.com/a/479974. 
#' 
#' @param r A SpatRaster to be rotated
#' @param angle Angle in degrees, with positive angle in clock-wise direction
#' @param point A sf POINT indicating the point around which rotation should
#' be done. if NULL, centroid of r will be used.
#' 
#' @returns a rotated SpatRaster
#' 
#' @export
rotate_raster <- function(r, 
                          angle,
                          point = NULL) {
    
    # rotation (affine transformation) using stars
    r <- stars::st_as_stars(r)
    
    # get the the resolution for x and y dimensions
    scale_x <- stars::st_dimensions(r)$x$delta
    scale_y <- stars::st_dimensions(r)$y$delta
    
    # get the translation from the lower left corner of the bbox
    trans_x <- sf::st_bbox(r)[1] |> as.numeric()
    trans_y <- sf::st_bbox(r)[2] |> as.numeric()
    
    # bbox of original raster
    bbox <- sf::st_bbox(r)
    sfc_bbox <- sf::st_as_sfc(bbox)
    
    test <- is.null(point)
    if(test) point <- sf::st_centroid(sfc_bbox)
    
    # convert to radians:
    rad <- angle * (pi / 180)
    
    # define affine parameters for rotation around centroid
    a <- scale_x * -sin(rad)
    b <- scale_x * cos(rad)     
    c <- scale_y * cos(rad)   
    d <- scale_y * sin(rad)
    e <- trans_x 
    f <- trans_y 
    
    # apply transformation
    r_rot <- r
    stars::st_geotransform(r_rot) = c(e, b, a, f, d, c)
    
    # convert to rectilinear grid
    r_rot <- stars::st_warp(r_rot, crs = sf::st_crs(r))

    # set bounding box using terra and ensure the original and rotated
    # have same crs and resolution
    bbox_rot <- (sfc_bbox - point) * rotation_matrix(angle) + point
    ext_rot <- terra::ext(sf::st_as_sf(bbox_rot))
    r_rot <- terra::rast(r_rot)
    terra::ext(r_rot) <- terra::ext(ext_rot)
    r <- terra::extend(terra::rast(r), ext_rot)
    r_rot <- terra::resample(r_rot, r)
    r_rot <- terra::trim(r_rot)
    
    return(r_rot)
}
