# l_infty projection
proj_linf_ball <- function(x, r) {
  pmax(pmin(x, r), -r)
}