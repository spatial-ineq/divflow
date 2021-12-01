

# setting up graphs -------------------------------------------------------------

#' rid2sfg
#'
#' @inheritDotParams load.sfg
#'
rid2sfg <- function(cz = NULL, cbsa = NULL
                    , ...) {

  require(tidyverse)
  requireNamespace('sfg.seg')
  requireNamespace('geox')

  .countyfps <- geox::x2cos(cz=cz, cbsa = cbsa)
  .czs <- geox::rx %>% filter(countyfp %in% .countyfps) %>% pull(cz) %>% unique()

  sfg <- load.sfg(.czs, ...)

  return(sfg)
}



#' sfx2sfg
#'
#' From a given `sf` object, get correspondence with commuting zones, and load all
#' safegraph data for trips with origins at those commuting zones. For now, just done
#' for annual, but we could make this monthly. If you want to start with just CZs,
#' just use sfg.seg::read.sfg.CZs.
#'
#' @param sfx an sf object within bounds of which to load all flows.
#' @inheritDotParams load.sfg
#'
#' @export sfx2sfg
sfx2sfg <- function( sfx
                     , year = 2019
                     , ...) {

  require(tidyverse)
  require(sf)
  requireNamespace('sfg.seg')
  requireNamespace('geox')

  # get county overlap with bounding area, and use that to load sfg data
  cos <- tigris::counties(year = year)
  cos <- cos %>% st_transform(st_crs(sfx))

  ovcos <- st_crop(cos
                   ,sfx)

  .czs <-
    geox::rx %>%
    filter(countyfp %in% ovcos$GEOID) %>%
    pull(cz) %>% unique()

  sfg <- load.sfg(.czs
                  ,year = year
                  ,...)

  return(sfg)
}


#' load.sfg
#'
#' @param min.flows flow floor to immediately filter to. Default 10.
#' @param trim.loops whether to trim loops, where origin==dest. Default true.
#' @inheritParams spatialize.graph
#' @param year year to load for. Default 2019.
#' @param ddir,sfg.dir location of safegraph data, saved by cz.
#'
load.sfg <- function(.czs
                     ,min.flows = 10
                     ,tracts.or.groups = c('ct', 'bg')
                     ,trim.loops = T
                     ,year = '2019'
                     ,ddir = Sys.getenv('drop_dir')
                     ,sfg.dir = 'sfg-processed/orig_dest_annual/'
                     ,...) {


  sfg <- sfg.seg::read.sfg.CZs(.czs
                               ,sfg.dir =
                                 paste0( ddir
                                         ,sfg.dir)
                               ,year = as.character(year)
  )

  sfg <- sfg %>%
    filter(n >= min.flows)

  # filter to trips starting and ending in overlapping CZs (not just starting there)
  sfg <- sfg %>%
    geox::geosubset(c('origin', 'dest')
                    ,cz = .czs)

  # aggregate to tracts if appropriate
  # note this has to be done before trimming loops
  if(tracts.or.groups[1] == 'ct')
    sfg <- sfg.seg::cbg.flows2tracts(sfg)

  # trim loops
  if(trim.loops)
    sfg <- sfg %>%
    filter(origin != dest)

  return(sfg)

}


#' sfg2gh
#'
#' From a origin-destination data.frame representing safegraph data, turn it into a
#' tidygraph object. Assumes `origin`, `dest`, and `n` columns. Also attaches
#' geometry and calculates "normalized connectedness," which will be different if
#' it's made undirected or left directed.
#'
#' @param sfg OD dataframe
#' @param directed whether to leave directed or make undirected (which aggregates
#'   trips to/from the same tracts, rather than duplicating the O-D pair for each
#'   direction.)
#'
#' @export sfg2gh
sfg2gh <- function(  sfg
                   , directed = F
                   , ...) {

  require(tidygraph)

  # get flow str in dif ways based on directedness
  if(directed) {

    # get straightforward %n_ij === n_ij/N_i
    sfg <- sfg %>%
      group_by(origin) %>%
      mutate(total.from = sum(n)) %>%
      group_by(dest) %>%
      mutate(total.inc = sum(n)) %>%
      ungroup() %>%
      mutate(  perc.to.dest = n / total.inc
               ,perc.from.origin = n / total.from)

    gh <- tidygraph::tbl_graph(
      edges = sfg,
      directed = directed
    )

  } else if(!directed) {

    # combine directed edges with igraph
    gh <- igraph::graph.data.frame(sfg, directed = T) # directed is true because it describes initial data.frame
    gh <- igraph::as.undirected(gh,
                                mode = "collapse",
                                edge.attr.comb = "sum")

    gh <- tidygraph::as_tbl_graph(gh)
    el <- gh %>%
      activate('edges') %>%
      as_tibble()

    # gets flows between the two tracts
    gh <- gh %>%
      activate('edges') %>%
      mutate(tstr =
               map2_dbl(from, to,
                        ~get.normalized.undirected.connectedness(.x, .y,
                                                                 el = el))
      )
  }

  return(gh)
}


#' spatialize.graph
#'
#' Spatializes a graph based on the node `name`s. Names should align with geoids for
#' either block groups or tracts (default).
#'
#' @param gh a graph, as setup by `sfg2gh`
#' @param frame.sf an `sf` object to crop to and use CRS from. Nodes outside the bbox of
#'   this area will be trimmed. Defaults to Lambert conformal conic projection if
#'   left null. (Why do I have all this complexity instead of just using coord_sf
#'   when mappi ng? Think I need to revise.. I think b/c i felt it allowed for more
#'   controlled trimming but..)
#' @inheritParams sfg2gh
#' @param tracts.or.groups Aggregate to census tract (ct) or leave at block group (bg)
#'
#' @export spatialize.graph
spatialize.graph <- function(gh
                             , frame.sf = NULL
                             , tracts.or.groups = c('ct', 'bg')
                             , directed = F
                             , year = 2019
                             , ...) {

  requireNamespace('sfnetworks')

  # get geometries based on type
  if(tracts.or.groups[1] == 'ct')
    tigris.call <- tigris::tracts
  else
    tigris.call <- tigris::block_groups

  eligible.ids <-
    gh %>% activate('nodes') %>% as_tibble()

  cbsf <- eligible.ids %>%
    rename(geoid = name) %>%
    geox::attach.geos(query.fcn = tigris.call
                      ,year = year)

  if(!is.null(frame.sf))
    cbsf <- cbsf %>%
    st_transform(st_crs(frame.sf)) %>%
    st_crop(frame.sf)
  else
    cbsf <- cbsf %>%
    #st_transform(4326)
    st_transform(crs = '+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45')

  # get centroids
  ctrs <- st_centroid(cbsf)

  # inner join to filter nodes out of crop distance
  gh <- gh %>%
    activate('nodes') %>%
    inner_join(ctrs
               ,by = c('name' = 'geoid'))

  # spatialize edges
  gh <- gh %>%
    activate('edges') %>%
    sfnetworks::as_sfnetwork(
       directed = directed
      ,edges_as_lines = T)

  # get lengths
  gh <- gh %>%
    activate('edges') %>%
    mutate(dst = st_length(geometry))

  return(gh)
}



#' apply.flow.filters
#'
#' Apply flow filters to a graph, as generated from `sfg2gh` and spatialized with
#' `spatialize.graph`. Crops to bounds.sf, filters by distance, total trips, and flow
#' strength based on parameters. Note that cropping to bounds sf is done to retain
#' edges that start outside of map area but ends within. This is done by cropping by
#' edges, which also creates a number of unnamed/NA-named nodes, which are artifacts
#' of nodes outside of cropping area when edges are partially inside.
#'
#' @inheritParams spatialize.graph
#' @inheritParams sfg2gh
#' @param min.tie.str If not null, floor for tie strength
#' @param min.flows minimum flows
#' @param tie.str.drop.deciles If not null, number of deciles of tie str to keep after all
#'   other trims
#'
#' @export apply.flow.filters
apply.flow.filters <- function(gh
                               ,tie.str.drop.deciles = 5
                               ,directed = F
                               ,frame.sf = NULL
                               ,min.tie.str = NULL
                               ,min.flows = 10
                               ,max.dst = units::set_units(10, 'miles')
                               ,...) {

  require(tidygraph)

  # browser()

  sfe <- gh %>% activate('edges') %>% as_tibble()

  if(! 'dst' %in% colnames(sfe) )
    gh <- spatialize.graph(gh)

  # crop edges that are not within bounds
  if(!is.null(frame.sf))
    gh <- gh %>%
    activate('edges') %>%
    st_crop(frame.sf)

  # filter by flows & distance
  if(!is.null(max.dst))
  gh <- gh %>%
    activate('edges') %>%
    filter(dst <= max.dst)

  if(!is.null(min.flows))
    gh <- gh %>%
    activate('edges') %>%
    filter(n >= min.flows)

  # filter by tie strength
  if(!is.null(min.tie.str))
    gh <- gh %>%
    activate('edges') %>%
    filter_at( vars(matches('^perc|^tstr'))
               ,any_vars(. >= min.tie.str))

  if(!is.null(tie.str.drop.deciles))
    gh <- gh %>%
    activate('edges') %>%
    filter_at( vars(matches('^perc|^tstr'))
               ,any_vars(. >= # drop n quartiles
                           quantile(.,
                                    seq(0,1,.1))[tie.str.drop.deciles + 1])
               )
  return(gh)
}


#' setup.gh.wrapper
#'
#' From a cz/cbsa identifier or an `sf` object, sets up a graph as all nodes surrounding
#' the area. Built to create a graph within a bounding box.
#'
#' @inheritDotParams sfx2sfg
#' @inheritDotParams sfg2gh
#' @inheritDotParams spatialize.graph
#' @inheritDotParams apply.flow.filters
#'
#' @export setup.gh.wrapper
setup.gh.wrapper <- function( sfx = NULL
                              , cz = NULL, cbsa = NULL
                              ,...) {

  if(!is.null(sfx)){
    sfg <- sfx2sfg(sfx = sfx
                   ,year = 2019
                   , ...)

    frame.sf <- st_bbox(sfx)
  } else {
    sfg <- rid2sfg(cz = cz, cbsa = cbsa,
                   ...)
    frame.sf <- NULL
    }

  gh <- sfg2gh(sfg = sfg
               ,...)

  gh <- spatialize.graph( gh = gh
                         , frame.sf = frame.sf
                         , ...)

  gh <- apply.flow.filters(gh
                           ,frame.sf = frame.sf
                           ,...)

  return(gh)
}



# setup for running ergm -------------------------------------------------------


#' get.normalized.undirected.connectedness
#'
#' Gets a measure of "tie strength" of total visits between tracts, scaled by total
#' flows into and out of each tract. Takes individual node/tract identifiers and an
#' edgelist containing all flows in between them.
#'
#'
#' @param x,y node ids
#' @param el edge list
#' @param flow.colm column containing edge weights/flows between nodes.
#'
get.normalized.undirected.connectedness <- function(x, y,
                                                    el, flow.colm = "n") {

  # subset to ods containing either of the two tracts
  total.trips <-
    el[with(el,
            from %in% c(x, y) |
              to %in% c(x, y) ), ]

  # vs flow between just the two tracts
  xy.link <-
    el[with(el,
            (from == x & to == y |
               from == y & to == x)), ]

  # get ratio
  flow.btwn <- pull(xy.link, flow.colm) %>% sum()
  total.trips <- pull(total.trips, flow.colm) %>% sum()

  flow.btwn / total.trips
}



#' ergm.clean.n.convert
#'
#' the ergms packages can't handle some types of information, like factors and NAs.
#' This applies final cleans neccessary for an ERGM to run without error and then
#' coerces to `network` object.
#'
#' @param population.floor Applies a population floor. Will assume a node-based `pop`
#'   column to filter by
#'
#' @return a `network` object for use with statnet libraries.
#'
ergm.clean.n.convert <- function(gh, pop.floor = 0) {

  #browser()

  # factors -> ints
  gh <- gh %>%
    activate("nodes") %>%
    mutate(across(where(is.factor),
                  as.integer)) %>%
    activate("edges") %>%
    mutate(across(where(is.factor),
                  as.integer))

  # apply population floor
  if(is.numeric(pop.floor))
    gh <- gh %>%
      activate("nodes") %>%
      filter(pop > pop.floor)

  # no NAs
  gh <- gh %>%
    activate("nodes") %>%
    mutate(across(where(is.numeric),
                  ~ifelse(is.na(.x),
                          0, .x)))

  # coerce to statnet network object.
  # (network object can't deal with
  ergh <- intergraph::asNetwork(gh)

  return(ergh)
}



#' is.in.city.center
#'
#' Given a tidy graph with tract centroid geometries, get each nodes' status as
#' within or outside the city center, defined as the largest Place by population
#' in the commuting zone. References a dataset bundled with `divM` that already
#' pairs these population centers with commuting zones.
#'
#' @param sfx a tidy graph with node centroids in a geometry column
#' @param cz a cz corresponding with the graph; if left `NULL`, fcn will assume
#'   a column with the cz id exists for nodes already
#' @param filter2cc Whether or not to filter or just add column
#'
#' @return a tidy graph with new `in.cc` column for nodes, based on whether or
#'   not they are in a population center. If filter2cc is `TRUE`, only keeps
#'   nodes within population center
#'
#' @export is.in.city.center
is.in.city.center <- function(sfx, cz = NULL, filter2cc = FALSE, ...) {

  requireNamespace("divM")
  require(sf)

  if(is.null(cz))
    cz <- unique(sfx$cz)

  # ensure in sf
  sfx <- st_sf(sfx)

  ccenter <-
    divM::largest.plc.in.cz %>%
    filter(cz.id %in% cz) %>%
    st_transform(st_crs(sfx))

  sbgp <-
    st_intersects(
      st_point_on_surface(sfx)
      , ccenter)

  sfx <- sfx %>%
    mutate(in.cc =
             lengths(sbgp) > 0)

  if(filter2cc)
    sfx <- sfx %>%
    filter(in.cc)

  return(sfx)
}



# convenience fcns -------------------------------------------------------------

#' gh2nodes
#'
#' convenience fcn to get nodes as tibble from tidygraph
#'
gh2nodes <- function(gh) {
  require(tidygraph)
  gh %>%
    activate('nodes') %>%
    as_tibble()
}

#' gh2edges
#'
#' convenience fcn to get edges as tibble from tidygraph
#'
gh2edges <- function(gh) {
  require(tidygraph)
  gh %>%
    activate('edges') %>%
    as_tibble()
}

gh2coords <- function(gh) {
  require(sf)
  require(tidygraph)
  gh %>%
    gh2nodes() %>%
    st_sf() %>%
    st_coordinates()


}


# resources and refernces -----------------------------------------------------------


#' http://statnet.org/Workshops/valued.html#32_modeling_ordinal_relational_data_using_ergmrank
#'
#' https://cran.r-project.org/web/packages/ergm/vignettes/ergm-term-crossRef.html
