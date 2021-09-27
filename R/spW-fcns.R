
# organizational note -----------------------------------------------------

#' most fcns that would be here are just moved to the sending Della job scripts
#' that generate the weights or composites.
#'
#' In R/creating SpWs/sending della jobs/


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
#' Get "flow weights," i.e., proportion of visitor/visited from other tracts
#' within distance threshold. Will use a precalculated distance matrix to filter
#' by distance if `spatial.weights` is supplied. Does not filter by distace if
#' it's left as null. Rounds to 4 digits for sparseness.
#'
#' @inheritParams get.dist.weighted.composite
#' @param flow.type Whether to get weights based on all visitors to i from other
#'   proximate areas (default); or based on all proximiate areas visited by
#'   people in area i.
#' @param drop.loops Whether to drop loops (where origin==destination)
#' @param prox.col name for list-column in `spws` with all areas within distance
#'   cutoff from each other
#' @param flow.counts data.frame with origin/destination/estimated visits (n).
#' @param weight.floor flow weight floor. Minimum percent of incoming/visited
#'   trips between tracts to be included. 0.1% by default.
#'
#' @export get.flow.weights.within.distance
get.flow.weights.within.distance <- function(i
                                             ,flow.type = c('incoming', 'visited')
                                             ,drop.loops = T
                                             ,prox.col =  'below.cutoff'
                                             ,flow.counts = sfg
                                             ,spatial.weights = NULL
                                             ,weight.floor = .001) {

  flow.type <- flow.type[[1]]

  # get geoids for js within distance threshold or within same region
  if(!is.null(spatial.weights))
    j.ids <- spatial.weights  %>% filter(geoid %in% i) %>% pull(prox.col) %>% unlist()
  else {
    j.ids <- c( unique(flow.counts$origin)
                ,unique(flow.counts$dest) )
  }

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

#' Della.wrapper_flow.weights_dst.cutoff
#'
#' Wrapper function to generate all flow weights. In addition to parameters,
#' keep a spatial weights object `spws` in global env. This is a a proximity
#' list as calculated with the `creating SpWs` scripts and needs the
#' `below.cutoff` list-column with tract/blockgroup geoids within the distance
#' maximum.
#'
#' @inheritParams get.flow.weights.within.distance
#'
Della.wrapper_flow.weights_dst.cutoff <- function(cz
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

  # drop loops if appropriate
  if(drop.loops)
    sfg <- sfg %>% filter(origin != dest)

  # agg2 tracts
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
                 )))
  flwws <- flwws %>%
    mutate(visited.flow.weights =
             map(geoid
                 ,~get.flow.weights.within.distance(i = .x
                                                    ,'visited'
                                                    ,flow.counts = sfg
                                                    ,spatial.weights = local.spws
                                                    ,weight.floor = 0.001
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


#' Della.wrapper_flow.weights_by.rt
#'
#' Wrapper function to generate all flow weights by either CZ or CBSA.
#'
#' @inheritParams get.flow.weights.within.distance
#' @param list.colms whether to return list columns (or tibble colms)
#'
Della.wrapper_flow.weights_by.rt <- function(
   cz_id = NULL
  ,cbsa_id = NULL
  ,agg2tracts = T
  ,drop.loops = T
  ,weight.floor = .001
  ,sfg.dir
  ,year = '2019'
  ,save.dir = NULL
  #,list.colms = T
) {

  if(is.null(cbsa_id))
    czs2load <- cz_id
  else
    czs2load <- geox::rx %>%
      filter(cbsa %in% cbsa_id) %>%
      pull(cz) %>% unique()

  # sfg load
  sfg <- sfg.seg::read.sfg.CZs(czs2load,
                               sfg.dir = sfg.dir,
                               year = year)

  # subset to trips w/in region (cz/cbsa)
  sfg <- geox::geosubset(sfg
                         ,subset.cols = c('origin', 'dest')
                         ,cz = cz_id
                         ,cbsa = cbsa_id)

  # agg2 tracts
  if(agg2tracts)
    sfg <- sfg.seg::cbg.flows2tracts(sfg)

  # drop loops if appropriate
  if(drop.loops)
    sfg <- sfg %>% filter(origin != dest)


  # get inc and visited weights
  inc.flwws <- sfg %>%
    group_by(dest) %>%
    mutate(flww = n / sum(n)) %>%
    ungroup() %>%
    filter(flww >= weight.floor) %>%
    select(geoid = dest, j = origin, flww) %>%
    nest(inc.flwws = c(j, flww))

  vis.flwws <- sfg %>%
    group_by(origin) %>%
    mutate(flww = n / sum(n)) %>%
    ungroup() %>%
    filter(flww >= weight.floor) %>%
    select(geoid = origin, j = dest, flww) %>%
    nest(vis.flwws = c(j, flww))

  flwws <- full_join(inc.flwws, vis.flwws)

  # browser()

  if(!is.null(save.dir)) {
    if(!exists(save.dir))
      dir.create(save.dir, recursive = T)

    rids <- region2identifiers( region.type = region.type
                                ,region.id = region.id)
    write_rds(flwws
              ,file = paste0(save.dir
                             ,rids
                             ,'-flow-weights.rds'))
  }

  return(flwws)
}


# helpers -----------------------------------------------------------------

#' region2identifiers
#'
#' Turns region.type and region.id into a cleaned & concatenated set of
#' identifiers, which I use for partial saves.
#'
#' @export region2identifiers
region2identifiers <- function(region.type, region.id) {

  requireNamespace('geox')
  rids <- geox::add.rns(tibble(rt = region.type, rid = region.id ))
  rids <- rids %>% paste0(collapse = '-')

  # transform '/'s in identifier (sometimes appear in CBSA names)
  rids <- rids %>% gsub('\\/', '-', .)
  return(rids)
}
