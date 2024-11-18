#' @export
rotation = function(a){
    r = a * pi / 180 #degrees to radians
    matrix(c(cos(r), sin(r), -sin(r), cos(r)), nrow = 2, ncol = 2)
} 

#' @export
ang <- function(x,y) { 
    z <- x + 1i * y
    res <- 90 - Arg(z) / pi * 180
    res %% 360
}

#' @export
bearing_2d <- function(mat) {
    theta <- atan2(mat[2,1] - mat[1,1], mat[2,2] - mat[1,2])
    rad <- ifelse(theta < 0, 2*pi + theta, theta) 
    deg <- rad2deg(rad)
    return(deg)
}

#' @export
rad2deg <- function(rad) {
    (rad * 180) / (pi)
}

#' @export
move_coords <- function(bear, dist) {
    movex <- dist * sin(bear*pi/180)
    movey <- dist * cos(bear*pi/180)
    return(cbind(x=movex, y=movey))
}

