


# making ggplot layers ----------------------------------------------------------------

flow.map.base <- function() {

  require(ggplot2)
  require(ggraph)
  list(
    scale_edge_alpha_continuous(guide = 'none'
                                ,range = c(.1,1))
    ,scale_edge_width_continuous(guide = 'none'
                                ,range = c(.1,2))
    ,scale_color_discrete(guide = 'none')
    ,theme_void()
    )
}


#' flow.map.wrapper
#'
#' Big old wrapper function. Includes div layers: Place, US routes, interstates, and
#' water.
#'
#' TODO: I think main issue here is that including outside of the bounds can mess up
#' what remains in the bounds... -> actually I think it's an issue of plot size and
#' the scale conversion to pixels.. Lines are thinner than a pixel and become
#' invisible depending on plot size. Need to think about how to handle that.
#' See: https://www.tidyverse.org/blog/2020/08/taking-control-of-plot-scaling/
#'
#' @param gh A graph object. Recreates based on `sfx` if not supplied
#' @param edge.attr Edge attribute to map by. 'n' by default. If generating graph
#'   with these functions, tstr is tie strength for undirected graph and perc.to.dest
#'   or perc.from.origin is for directed graph.
#' @inheritDotParams setup.gh.wrapper
#'
#' @export flow.map.wrapper
flow.map.wrapper <- function(sfx
                             ,gh = NULL
                             ,edge.attr = 'n'
                             ,log.transform = F
                             ,map.buffer = units::set_units(5, 'miles')
                             ,max.dst = units::set_units(5, 'miles')
                             ,...) {

  browser()

  ctr <- sfx %>%
    divM::conic.transform() %>%
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
    gh2coords()

  ggbase <- flow.map.base()

  # transform if necessary
  gh <- gh %>%
    mutate(ex = get(edge.attr))
  if(log.transform)
    gh$ex <- gh %>% gh2edges() %>% pull(ex) %>% log()

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
    sfx2coord_sf(bounds.area)

  return(flow.map)
}


#' sfx2coord_sf
#'
#'
#' moved to visaux::bbox2ggcrop
#'
sfx2coord_sf <- function(sfx, expand = F) {

  # get bbox for crop
  zoom.box <- st_bbox(sfx)

  coord_sf( xlim =
              c( zoom.box$xmin
                 ,zoom.box$xmax)
            , ylim =
              c(zoom.box$ymin
                ,zoom.box$ymax)
            ,expand = expand
            )
}
