#' @export
rotate_raster <- function(r, 
                          centroid, 
                          angle,
                          trim = FALSE) {
    
    # inputs r, angle, centroid around which to rotate
    # image magick
    maxval <- unlist(terra::global(r, max, na.rm=TRUE))
    minval <- unlist(terra::global(r, min, na.rm=TRUE))
    img <- as.array(r) / maxval
    img <- abind::abind(img, img, along = 3)
    img <- abind::abind(img, img[,,1], along = 3)
    img <- magick::image_read(img) 
    img <- magick::image_rotate(img, angle)
    r_rot <- raster::as.raster(img) |> as.matrix() |> terra::rast()
    
    # new extent
    # dims <- dim(r_rot) - dim(r)
    extent <- terra::ext(r)
    ext_rot <- extent |> 
        terra::as.polygons() |> 
        sf::st_as_sf() |>
        sf::st_set_crs(sf::st_crs(r)) |> 
        sf::st_geometry()
    ext_rot <- (ext_rot - centroid) * rotation(angle) + centroid
    ext_rot <- ext_rot |> 
        terra::vect() 
    
    # set extent and resample to make sure rasters align and crop artefacts
    r_ext <- terra::extend(r, terra::ext(ext_rot), snap = "out")
    terra::ext(r_rot) <- terra::ext(ext_rot)
    resolution <- max(terra::res(r_rot))*2
    r_rot <- r_rot |> 
        terra::resample(r_ext, method = "near") 
    
    # trim values from 
    if(trim) {
        r_rot <- terra::mask(r_rot, 
                             terra::vect(sf::st_buffer(sf::st_as_sf(ext_rot), -resolution)))
    } else {
        r_rot <- terra::mask(r_rot,ext_rot)
    }
    
    # remove colour table 
    terra::coltab(r_rot) <- NULL
    # r_rot <- (r_rot + (minval-1)) * (maxval / 255)
    r_rot <- terra::stretch(r_rot, minv = minval, maxv = maxval)
    
    return(r_rot)
}