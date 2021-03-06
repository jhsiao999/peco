#' @title Pairwise Distance Between Two Circular Variables
#'
#' @description We define distance between two angles as the minimum of
#'   the differences in both clockwise and counterclockwise directions.
#'
#' @param y1 A vector of angles.
#' 
#' @param y2 A vector of angles.
#'
#' @return A vector of distances between the angles.
#'
#' @examples
#' # Vector of angles.
#' theta_ref <- seq(0,2*pi, length.out=100)
#'
#' # Shift the origin of theta_ref to pi.
#' theta_compare <- shift_origin(theta_ref, origin = pi)
#' mean(circ_dist(theta_ref, theta_compare))
#'
#' # After rotation of angles, difference is zero between the original
#' # and the shifted angles.
#' theta_compare_rotated <- rotation(ref_var=theta_ref,
#'     shift_var=theta_compare)
#' mean(circ_dist(theta_ref, theta_compare_rotated))
#'
#' @author Joyce Hsiao, Matthew Stephens
#' 
#' @export
#' 
circ_dist <- function(y1,y2)
    pmin(abs(y2 - y1), abs(2*pi - abs(y2 - y1)))

#' @title Rotate Circular Variable to Minimize Distance Between
#'   Variables.
#'
#' @description Rotate circular variable \code{shift_var} to minimize
#'   the distance between \code{ref_var} and \code{shift_var}. Because
#'   the origin of the cell cycle phases is arbitrary, we transform the
#'   angles prior to computing the distance (rotation and
#'   shifting). After this, one can apply \code{circ_dist} to compute
#'   the distance between the output value and \code{ref_var}.
#'
#' @param ref_var A vector of reference angles.
#' 
#' @param shift_var A vector of angles to be compared to ref_var.
#'
#' @return The transformed value of \code{shift_var} after rotation and
#'   shifting.
#'
#' @examples
#'
#' # Create a vector of angles.
#' theta_ref <- seq(0,2*pi, length.out=100)
#' 
#' # Shift the origin of theta_ref to pi.
#' theta_compare <- shift_origin(theta_ref, origin = pi)
#' 
#' # Rotate theta_compare in a such a way that the distance between
#' # theta_ref and theta_compare is minimized.
#' theta_compare_rotated <-
#'   rotation(ref_var=theta_ref, shift_var=theta_compare)
#' 
#' par(mfrow = c(1,2))
#' plot(x=theta_ref, y = theta_compare)
#' plot(x=theta_ref, y = theta_compare_rotated)
#'
#' @author Matthew Stephens
#' @export
#'
rotation <- function(ref_var,shift_var) {
    df <- data.frame(flip=rep(c(1,-1), each=length(shift_var)),
                     shift = rep(shift_var, 2))

    for (i in seq_len(nrow(df))) {
        shift_var_tmp <- df$flip[i] * ((shift_var + df$shift[i]) %% (2*pi))
        df$dis[i] <- mean(circ_dist(ref_var, shift_var_tmp))
    }
    which_cutoff <- which.min(df$dis)

    shift_var_new <-
      df$flip[which_cutoff] * ((shift_var + df$shift[which_cutoff]) %% (2*pi))

    return(shift_var_new)
}
