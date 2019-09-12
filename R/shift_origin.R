#' @title Shift origin of the angles
#'
#' @description Shift origin of the angles for visualization
#'
#' @param phase A vector of angles (in radians).
#'
#' @param origin the new origin of the angles.
#'
#' @return A vector of angles shifted to the new origin.
#'
#' @author Joyce Hsiao
#'
#' @export
#'
shift_origin <- function(phase, origin) {
  phase2 <- phase
  if (origin > 0) {
    phase2[phase>=origin] <- phase[phase>=origin]-origin
    phase2[phase<origin] <- (phase[phase<origin]-origin+2*pi)
  }
  if (origin < 0) {
    a <- phase[phase >= 0 & phase< (2*pi-abs(origin))] + abs(origin)
    b <- -1*(2*pi - (phase[phase <= 2*pi & phase >= (2*pi - abs(origin))])) + abs(origin)
    phase2 <- c(a,b)
  }
  phase2
}
