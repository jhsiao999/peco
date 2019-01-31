#' Pairwise distance between two circular variables
#'
#' We define distance between two angles: the minimum of the differences
#' in both clockwise and counterclockwise directions.
#'
#' @param y1 A vector of angles
#' @param y2 A vector of angles
#'
#' @return A vector of distances between angles
#'
#' @author Joyce Hsiao, Matthew Stephens
#'
#' @export
circ_dist <- function(y1,y2){
  pmin(abs(y2-y1), abs(2*pi-(abs(y2-y1))))
}


#' Rotate circular variable shift_var to minimize distance between ref_var and shift_var
#'
#' Because the origin of the cell cycle phases is arbitrary, we transform
#' the angles prior to computing the distance (rotation and shifting)
#' to minimize the distance between two vectors. After this, one can apply
#' circ_dist to compute the distance between the output value and ref_var.
#'
#' @param ref_var A vector of reference angles
#' @param shift_var A vector of angles to be compared to ref_var
#'
#' @return  The transformed values of shift_var after rotation and shifting
#'
#' @author Matthew Stephens
#'
#' @export
rotation <- function(ref_var,shift_var){

  df <- data.frame(flip=rep(c(1,-1), each=length(shift_var)),
                   shift = rep(shift_var, 2))

  for (i in 1:nrow(df)) {
    shift_var_tmp <- df$flip[i]*((shift_var+df$shift[i])%%(2*pi))
    df$dis[i] <- mean(circ_dist(ref_var, shift_var_tmp))
  }
  which_cutoff <- which.min(df$dis)

  shift_var_new <- df$flip[which_cutoff]*((shift_var+df$shift[which_cutoff])%%(2*pi))

  return(shift_var_new)
}




