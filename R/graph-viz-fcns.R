


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
#' @param bounds.sf sf area that defines bbox to get division layers from.
#' @param ... NHPN directory info passed onto `visaux::get.NHPN`
#'
#' @export get.div.layers
get.div.layers <- function( bounds.sf
                           , ...) {

  div.lyrs <- visaux::add.map.layers(bounds.sf
                                     ,add.counties = NULL
                                     ,add.places = "#583799"
                                     ,lwd = .5
                                     ,spatial.trim = st_crop)
  nhpn <- visaux::get.NHPN(st_transform(bounds.sf# need to let it take ... :,( # also need to make that fcn better so it's not crs dependent
                                        ,4326)) # ...)
  hwys <- nhpn %>% filter(signt1 %in% c("I", "U"))

  div.lyrs$ints <-
    geom_sf(data = filter(hwys,
                          signt1=="I")
            , color = "#7d3636", size = .9)
  div.lyrs$ushwys <-
    geom_sf(data = filter(hwys,
                          signt1=="U")
            , color = "#ed7272", size = .6)

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
                             ,log.transform = T
                             ,map.buffer = units::set_units(5, 'miles')
                             ,max.dst = units::set_units(5, 'miles')
                             ,...) {

  ctr <- sfx %>%
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

  # get bbox for crop
  zoom.box <- st_bbox(bounds.area)

  ggbase <- flow.map.base()

  # transform if necessary
  gh <-  gh %>%
    mutate(ex = get(edge.attr))
  if(log.transform)
    gh$ex <- log(gh$ex)

  require(ggraph)

  flow.map <-
    gh %>%
    ggraph::ggraph(layout = lonlats ) +
    geom_edge_density(
      aes(edge_fill = ex)
    ) +
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