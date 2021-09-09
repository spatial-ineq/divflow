#' get.dist.weighted.composite
#'
#' Given an area geoid, a set of attribute variables for all areas, and a list-column
#' of spatial weights between areas, calculate a "spatial composite" that
#' incorporates neighbors' attributes, spatial weights, and population or other
#' weights.
#'
#' @param i a given tract geoid
#' @param x full ct data to subset from
#' @param value.col Attribute column name in `x` to calculate weighted avg of
#'   adjacent tracts for.
#' @param weight.col column name in `x` for non spatial weights. Probably population
#'   or households
#' @param dist.decay.fcn a distance decay or proximity fcn to convert distances to
#'   spatial weights
#' @param dist.col string for column name in `spws` that pairs geoids with distances
#'   from i
#' @param spws spatial weight matrices with a geoid column and other list columns for
#'  spatial weights
#' @param ... passed on to `dst.decay.fcn`
#'
#' @export get.dist.weighted.composite
get.dist.weighted.composite <- function(i,  x
                                        , value.col = 'value'
                                        , weight.col = 'weight'
                                        , dist.decay.fcn
                                        ,dist.col = 'dists'
                                        ,spatial.weights = spws
                                        ,...) {

  # all neighbors within dists of i, based on supplied spatial weights
  nbs <- spatial.weights %>% filter(geoid %in% i) %>% pull(dist.col)
  # organize as tibble
  nbs <- tibble(geoid = names(nbs[[1]])
                ,dist.from.i = nbs[[1]] )

  # combine distances and attributes; filter out loops (where i==j)
  js <- x[x$geoid %in% nbs$geoid, ] %>%
    left_join(nbs, by = 'geoid') %>%
    filter(geoid != i)

  # convert distance to weight w/ proximity fcn
  js$spatial.weight <-
    dist.decay.fcn(js$dist.from.i
                   , ...)

  # apply both spatial weight and pop or hh weight to get spatial mean
  spu <- stats::weighted.mean(js$value
                              ,js$weight * js$spatial.weight
                              ,na.rm = T)
  return(spu)
}


#' dsts.below.threshold
#'
#' Quick fcn with hardcoded items, including assumption of meters from crs length
#' unit
#'
#' @param sfx an sf object to selectively get distances from i, probably
#'   containing point geometries.
#' @param i row # index of given tract/centroid in `sfx`
#' @param j.colm list-column of geoids in `sfx` to get distance from
#'
#' @examples
#' \dontrun{
#' ctrs <- cts %>% st_centroid()
#' dst.nbs <-
#'  st_is_within_distance(ctrs
#'                          , dist =
#'                          units::set_units(20,'miles'))
#' cts$nbs <- map(dst.nbs, ~tmp[., ]$geoid)
#' ctrs %>%
#' mutate(dists =
#'          map(1:nrow(ctrs)
#'           ,~dsts.below.threshold(ctrs, ., 'nbs')))
#'}
#'
dsts.below.threshold <- function(sfx, i, j.colm = 'nbs') {

  #  browser()
  # get sf corresponding to geoids j
  js <- pull(sfx, j.colm)[[i]]
  geos <- sfx[sfx$geoid %in% js, ]$geometry

  # get pairwise distance between i and j
  as.vector(st_distance(sfx[i, ]
                        , geos)) %>%
    set_units('m') %>%
    set_units('km') %>%
    round(digits = 2) %>%
    setNames(js)

}
