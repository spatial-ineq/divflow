
# sp-weighted composites -------------------------------------------------------



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



# distances within cutoff ------------------------------------------------------


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




# flow weights -----------------------------------------------------------------



#' get.flow.weights.within.distance
#'
#' Get "flow weights," i.e., proportion of visitor/visited from other tracts within
#' distance threshold. Relies on precalculated distance matrix. Rounds to 4 digits
#' for sparseness.
#'
#' @inheritParams get.dist.weighted.composite
#' @param flow.type Whether to get weights based on all visitors to i from other
#'   proximate areas (default); or based on all proximiate areas visited by people in
#'   area i.
#' @param drop.loops Whether to drop loops (where origin==destination)
#' @param prox.col name for list-column in `spws` with all areas within distance
#'   cutoff from each other
#' @param flow.counts data.frame with origin/destination/estimated visits (n).
#' @param weight.floor flow weight floor. Minimum percent of incoming/visited trips
#'   between tracts to be included. 0.1% by default.
#'
#' @export get.flow.weights.within.distance
get.flow.weights.within.distance <- function(i
                                             ,flow.type = c('incoming', 'visited')
                                             ,drop.loops = T
                                             ,prox.col =  'below.cutoff'
                                             ,flow.counts = sfg
                                             ,spatial.weights = spws
                                             ,weight.floor = .001) {

  flow.type <- flow.type[[1]]

  # get geoids for js within distance threshold
  j.ids <- spatial.weights  %>% filter(geoid %in% i) %>% pull(prox.col) %>% unlist()
  # drop loops if appropriate
  if(drop.loops)
    flow.counts <- flow.counts %>% filter(origin != dest)

  # Get either all trips from Js to I (incoming) or from I to all Js (visited)
  if(flow.type == 'incoming') {

    js <- flow.counts %>%
      filter(dest %in% i &
               origin %in% j.ids) %>%
      mutate(j = origin)

  } else if(flow.type == 'visited') {

    js <- flow.counts %>%
      filter(origin %in% i &
               dest %in% j.ids) %>%
      mutate(j = dest)
  } else
    stop("need 'incoming' or 'visited' as flow.type")

  # get flow weights
  flwws <- js %>%
    mutate(flow.weight =
             n / sum(n)) %>%
    select(j, flow.weight) %>%
    filter(flow.weight >=
             weight.floor)

  # return as named vector
  flwws$flow.weight %>%
    as.vector() %>%
    round(digits = 4) %>%
    setNames(flwws$j)
}

#' Della.wrapper_flow.weights
#'
#' Wrapper function to generate all flow weights. In addition to parameters, add a
#' psatial weights object `spws` with a proximity list, as calculated with the
#' `creating SpWs` scripts. This needs the `below.cutoff` list-column with
#' tract/blockgroup geoids within the distance maximum.
#'
#' @inheritParams get.flow.weights.within.distance
#'
Della.wrapper_flow.weights <- function(cz
                                       ,agg2tracts = T
                                       ,drop.loops = T
                                       ,weight.floor = .001
                                       ,sfg.dir
                                       ,year = '2019'
                                       ,save.dir = NULL
) {

  #browser()

  # filter spws to area within cz, and within distance range of those
  local.spws <- spws %>% geox::geosubset(cz = cz) %>% select(geoid, below.cutoff)
  in.range <- unlist(local.spws$below.cutoff) %>% unique()
  # to counties
  cos.in.range <- in.range %>% substr(1,5) %>% unique()

  # sfg2load
  czs2load <- geox::rx %>%
    filter(countyfp %in%
             unique(substr(in.range
                           ,1,5))
    ) %>%
    pull(cz) %>% unique()

  sfg <- sfg.seg::read.sfg.CZs(czs2load,
                               sfg.dir = sfg.dir,
                               year = year)

  if(agg2tracts)
    sfg <- sfg.seg::cbg.flows2tracts(sfg)

  # get inc and visited weights
  flwws <- local.spws %>%
    mutate(inc.flow.weights =
             map(geoid
                 ,~get.flow.weights.within.distance(i = .x
                                                    ,'incoming'
                                                    ,flow.counts = sfg
                                                    ,spatial.weights = local.spws
                                                    ,weight.floor = 0.001
                                                    ,drop.loops = drop.loops
                 )))
  flwws <- flwws %>%
    mutate(visited.flow.weights =
             map(geoid
                 ,~get.flow.weights.within.distance(i = .x
                                                    ,'visited'
                                                    ,flow.counts = sfg
                                                    ,spatial.weights = local.spws
                                                    ,weight.floor = 0.001
                                                    ,drop.loops = drop.loops
                 )))

  flwws <- flwws %>% select(geoid, matches('flow.weights'))


  if(!is.null(save.dir)) {
    if(!exists(save.dir))
      dir.create(save.dir)

    write_rds(flwws
              ,file = paste0(save.dir
                            ,as.character(cz), '-', as.character(year)
                            ,'-flow-weights.rds'))
    }

  return(flwws)
}

