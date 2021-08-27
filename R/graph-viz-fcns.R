

# setting up graph -------------------------------------------------------------


#' sfx2sfg
#'
#' From a given `sf` object, get correspondence with commuting zones, and load all
#' safegraph data for trips with origins at those commuting zones. For now, just done
#' for annual, but we could make this monthly. If you want to start with just CZs,
#' just use sfg.seg::read.sfg.CZs.
#'
#' @param eligible.sf an sf object within bounds of which to load all flows.
#'
#' @export sfx2sfg
sfx2sfg <- function( eligible.sf
                    ,min.flows = 10
                    ,base.dir = Sys.getenv('drop_dir')
                    ,sfg.dir = 'sfg-processed/orig_dest_annual/'
                    ,year= 2019
                    ,trim.loops = T) {

  requireNamespace('sfg.seg')
  requireNamespace('geox')

  cos <- tigris::counties(year = year)
  cos <- cos %>% st_transform(st_crs(eligible.sf))

  ovcos <- st_crop(cos
                   ,eligible.sf)

  czs2load <-
    geox::rx %>%
    filter(countyfp %in% ovcos$GEOID)

  czs2load <- unique(czs2load$cz)

  sfg <- sfg.seg::read.sfg.CZs(czs2load
                               ,sfg.dir =
                                 paste0( base.dir
                                         ,sfg.dir)
                               )
  sfg <- sfg %>%
    filter(n >= min.flows)

  # filter to trips starting and ending in overlapping CZs (not just starting there)
  sfg <- sfg %>%
    geox::geosubset(c('origin', 'dest')
                    ,cz = czs2load)

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
#' @param tracts.or.groups Aggregate to census tract (ct) or leave at block group (bg)
#' @param directed whether to leave directed or make undirected (which aggregates
#'   trips to/from the same tracts, rather than duplicating the O-D pair for each
#'   direction.)
#'
#' @export sfg2gh
sfg2gh <- function(sfg,
                   tracts.or.groups = c('ct', 'bg'),
                   directed = F) {

  requireNamespace('sfg.seg')
  require(tidygraph)

  # aggregate to tracts if appropriate
  if(tracts.or.groups[1] == 'ct')
    sfg <- sfg.seg::cbg.flows2tracts(sfg)

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
    gh <- igraph::graph.data.frame(od, directed = T) # directed is true because it describes initial data.frame
    gh <- igraph::as.undirected(dgh,
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
#' @param sfx an `sf` object to crop to and use CRS from. Nodes outside the bbox of
#'   this area will be trimmed.
#' @inheritParams sfg2gh
#'
#' @export spatialize.graph
spatialize.graph <- function( gh
                             ,frame.sf
                             ,tracts.or.groups = c('ct', 'bg')
                             ,directed = F
                             ,year = 2019) {

  # get geometries based on type
  if(tracts.or.groups[1] == 'ct')
    tigris.call <- tigris::tracts
  else
    tigris.call <- tigris::block_groups


  ellids <- c(sfg$origin, sfg$dest) %>% unique()
  full_cbsf <- tibble(geoid = ellcbs) %>%
    geox::attach.geos(query_fcn = tigris.call
                       ,year = year) %>%
    st_transform(st_crs(frame.sf))

  cbsf <- st_crop(full_cbsf,
                  frame.sf)

  # get centroids
  ctrs <- st_centroid(cbsf)

  # inner join to filter out of crop distance
  # (although we'll also filter by distance later)
  gh <- gh %>%
    activate('nodes') %>%
    inner_join(ctrs
               ,by = c('name' = 'geoid'))


  # spatialize edges
  gh <- gh %>%
    activate('edges') %>%
    sfnetworks::as_sfnetwork(
      .
      , directed = directed
      , edges_as_lines = T)

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
#' @param min.str If null, filters to top quartile after other filters applied
#'
#' @export apply.flow.filters
apply.flow.filters <- function(gh,
                               bounds.sf = NULL,
                               directed = F,
                               min.str = NULL,
                               min.flows = 10,
                               max.dst = units::set_units(10, 'miles')
) {

  require(tidygraph)

  sfe <- gh %>% activate('edges') %>% as_tibble()
  if(! 'dst' %in% colnames(sfe))
    gh <- spatialize.graph(gh)

  # crop edges that are not within bounds
  if(!is.null(bounds.sf))
    gh <- gh %>%
    activate('edges') %>%
    st_crop(bounds.sf)

  # filter by flows & distance
  gh <- gh %>%
    activate('edges') %>%
    filter(dst <= max.dst,
           n >= min.flows)

  # filter by tie strength
  if(!is.null(min.str))
    gh <- gh %>%
    activate('edges') %>%
    filter_at( vars(matches('^perc|^tstr'))
               ,any_vars(. > min.str))
  else
    gh <- gh %>%
      activate('edges') %>%
      filter_at( vars(matches('^perc|^tstr'))
                 ,any_vars(. > # drop n quartiles
                             quantile(.,
                                      seq(0,1,.25))[4]))
  return(gh)
}


#' setup.gh.wrapper
#'
#' From a `sf` object that will center graph, a buffer distance to define area around
#' center, and a set of trimming parameters, sets up a graph as all nodes surrounding
#' the area
#'
#' @inheritParams sfx2sfg
#' @inheritParams sfg2gh
#' @inheritParams spatialize.graph
#' @inheritParams apply.flow.filters
#'
#' @export setup.gh.wrapper
setup.gh.wrapper <- function( sfx
                             ,map.buffer = units::set_units(5, 'miles')
                             ,max.dst = units::set_units(5, 'miles')
                             ,min.flows = 10
                             ,min.str = NULL
                             ,directed = F
                             ,tracts.or.groups = c('ct', 'bg')
                             ,year = 2019
                             ,base.dir = Sys.getenv('drop_dir')
                             ,sfg.dir = 'sfg-processed/orig_dest_annual/'
                             ,crs = '+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45'
                             ) {

  ctr <- sfx %>%
    st_transform(crs) %>%
    st_centroid()

  eligible.area <- ctr %>%
    st_buffer(map.buffer + max.dst)

  bounds.area <- ctr %>%
    st_buffer(map.buffer)

  sfg <- sfx2sfg(eligible.area
                 ,min.flows = min.flows
                 ,base.dir = base.dir
                 ,sfg.dir = sfg.dir
                 ,year = year )

  gh <- sfg2gh(sfg,
               tracts.or.groups = tracts.or.groups,
               directed = directed)

  gh <- spatialize.graph(gh
                         ,sfx = eligible.area
                         ,tracts.or.groups = tracts.or.groups
                         ,directed = directed
                         ,year = year)

  gh <- apply.flow.filters(gh
                           ,bounds.sf = bounds.area
                           ,directed = directed
                           ,min.str = min.str
                           ,min.flows = min.flows
                           ,max.dst = max.dst)

  return(gh)
}





# making ggplot layers ----------------------------------------------------------------

flow.map.base <- function() {

  require(ggplot2)
  require(ggraph)
  list(
    scale_edge_alpha_continuous(guide = 'none'
                                ,range = c(.05,1))
    ,scale_edge_width_continuous(guide = 'none'
                                ,range = c(.05,2))
    ,scale_color_discrete(guide = 'none')
    ,theme_void()
    )
}


#' get.div.layers
#'
#' Wraps visaux helpers with very specific settings for combining everything.
#'
#' @export get.div.layers
get.div.layers <- function( bounds.sf
                           ,map.buffer = units::set_units(5, 'miles')
                           ,dropbox.dir = Sys.getenv("drop_dir")
                           ,subdir = "shapefiles/nhpn/") {

  div.lyrs <- visaux::add.map.layers(bounds.sf
                                     ,add.counties = NULL
                                     ,add.places = "#583799"
                                     ,lwd = .5
                                     ,spatial.trim = st_crop)
  nhpn <- visaux::get.NHPN(bounds.sf)
  hwys <- nhpn %>% filter(signt1 %in% c("I", "U"))

  div.lyrs$ints <-
    geom_sf(data = filter(hwys,
                          signt1=="I")
            , color = "#fa7f7f", size = .9)
  div.lyrs$ushwys <-
    geom_sf(data = filter(hwys,
                          signt1=="U")
            , color = "#882222", size = .6)

  return(div.lyrs)
}



#' flow.map.wrapper
#'
#' Big old wrapper function. Includes div layers: Place, US routes, interstates, and
#' water.
#'
#' TODO: I think main issue here is that including outside of the bounds can mess up
#' what remains in the bounds...
#'
#' @param gh A graph object. Recreates based on `sfx` if not supplied
#' @param edge.attr Edge attribute to map by. 'n' by default. If generating graph
#'   with these functions, tstr is tie strength for undirected graph and perc.to.dest
#'   or perc.from.origin is for directed graph.
#' @inheritDotParams setup.gh.wrapper
#' @inheritDotParams get.div.layers
#'
#'
#' @export flow.map.wrapper
flow.map.wrapper <- function(sfx
                             ,gh = NULL
                             ,edge.attr = 'n'
                             ,map.buffer = units::set_units(5, 'miles')
                             ,max.dst = units::set_units(5, 'miles')
                             ,...) {

  ctr <- sfx %>%
    st_transform(crs) %>%
    st_centroid()

  eligible.area <- ctr %>%
    st_buffer(map.buffer + max.dst)

  bounds.area <- ctr %>%
    st_buffer(map.buffer)

  if(is.null(gh))
    gh <- setup.gh.wrapper( sfx
                           , map.buffer = map.buffer
                           , max.dst = max.dst
                           , ...)

  lonlats <- gh %>%
    activate('nodes') %>%
    as_tibble() %>%
    st_sf() %>%
    st_coordinates()

  div.layers <- get.div.layers(frame.sf, ...)

  # get bbox for crop
  zoom.box <- st_bbox(frame.sf)

  ggbase <- flow.map.base()

  require(ggraph)

  flow.map <-
    gh %>%
    mutate(ex = get(edge.attr)) %>%
    ggraph::ggraph(layout = lonlats ) +
    geom_edge_density(
      aes(edge_fill = ex)
    ) +
    div.layers +
    geom_edge_fan(aes(edge_alpha =
                        ex
                      ,edge_width =
                        ex)) +
    ggbase +
    coord_sf( xlim =
                c( zoom.box$xmin
                   ,zoom.box$xmax)
              , ylim =
                c(zoom.box$ymin
                  ,zoom.box$ymax)
              ,expand = F
    )


  return(flow.map)
}
