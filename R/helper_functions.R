#' Create rotation matrix for rotation of coordinates
#' 
#' 
#' @param a Angle in degrees, with positive angles in clockwise direction.
#' 
#' @returns a rotated SpatRaster
#' 
#' @export
rotation_matrix <- function(a) {
    r = a * pi / 180 #degrees to radians
    matrix(c(cos(r), sin(r), -sin(r), cos(r)), nrow = 2, ncol = 2)
} 


#' Compute the bearing between two points
#' 
#' Computes a 2d bearing (angle) of a line between two points.
#' 
#' 
#' @param mat A 2x2 matrix with X and Y coordinate in columns, and start and end 
#'   point as rows. 
#' 
#' @returns The bearing (angle) in degrees. 0 indicates North.
#' 
#' @export
bearing_2d <- function(mat) {
    theta <- atan2(mat[2,1] - mat[1,1], mat[2,2] - mat[1,2])
    rad <- ifelse(theta < 0, 2*pi + theta, theta) 
    deg <- rad2deg(rad)
    return(deg)
}

#' Convert radians to degrees
#' 
#' @param rad A vector of radians to convert 
#' 
#' @returns Vector of degrees
#' 
#' @export
rad2deg <- function(rad) {
    (rad * 180) / (pi)
}

#' Compute displacement in X and Y coordinates
#' 
#' Compute displacement in X and Y coordinates based on bearing angle and 
#' distance along the bearing angle.
#' 
#' @param bear Angle of bearing
#' @param dist distance 
#' 
#' @returns A vector of two values: displacement in X and Y coordinate.
#' 
#' @export
move_coords <- function(bear, dist) {
    movex <- dist * sin(bear*pi/180)
    movey <- dist * cos(bear*pi/180)
    return(cbind(x=movex, y=movey))
}

