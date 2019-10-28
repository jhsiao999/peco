#' @name circ_dist
#'
#' @title Pairwise distance between two circular variables
#'
#' @description We define distance between two angles: the minimum of
#' the differences in both clockwise and counterclockwise directions.
#'
#' @param y1 A vector of angles.
#' @param y2 A vector of angles.
#'
#' @return A vector of distances between angles.
#'
#' @examples
#' # a vector of angles
#' theta_ref <- seq(0,2*pi, length.out=100)
#'
#' # shift the origin of theta_ref to pi
#' theta_compare <- shift_origin(theta_ref, origin = pi)
#' mean(circ_dist(theta_ref, theta_compare))
#'
#' # after rotation of angles, difference is 0 between the original
#' # and the shifted angles
#' theta_compare_rotated <- rotation(ref_var=theta_ref,
#'     shift_var=theta_compare)
#' mean(circ_dist(theta_ref, theta_compare_rotated))
#'
#'
#' @author Joyce Hsiao, Matthew Stephens
#' @export
circ_dist <- function(y1,y2) {
    pmin(abs(y2-y1), abs(2*pi-(abs(y2-y1))))
}

#' @name rotation
#'
#' @title Rotate circular variable shift_var to minimize distance
#' between ref_var and shift_var
#'
#' @description Because the origin of the cell cycle phases is
#' arbitrary, we transform the angles prior to computing the distance
#' (rotation and shifting) to minimize the distance between two
#' vectors. After this, one can apply circ_dist to compute the
#' distance between the output value and ref_var.
#'
#' @param ref_var A vector of reference angles.
#' @param shift_var A vector of angles to be compared to ref_var.
#'
#' @return The transformed values of shift_var after rotation and
#' shifting.
#'
#' @examples
#' # a vector of angles
#' theta_ref <- seq(0,2*pi, length.out=100)
#'
#' # shift the origin of theta_ref to pi
#' theta_compare <- shift_origin(theta_ref, origin = pi)
#'
#' # rotate theta_compare in a such a way that the distance
#' # between theta_ref and thet_compare is minimized
#' theta_compare_rotated <- rotation(ref_var=theta_ref, shift_var=theta_compare)
#'
#' par(mofrow=c(1,2))
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
        shift_var_tmp <- df$flip[i]*((shift_var+df$shift[i])%%(2*pi))
        df$dis[i] <- mean(circ_dist(ref_var, shift_var_tmp))
    }
    which_cutoff <- which.min(df$dis)

    shift_var_new <-
        df$flip[which_cutoff]*((shift_var+df$shift[which_cutoff])%%(2*pi))

    return(shift_var_new)
}
