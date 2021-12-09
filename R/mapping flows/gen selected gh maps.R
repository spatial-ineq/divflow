# ws ----------------------------------------------------------------------

library(tidyverse)
library(sf)

# option setting
sf_use_s2(T)
options(tigris_use_cache = TRUE)
Sys.setenv("VROOM_SHOW_PROGRESS"="false")

devtools::load_all()

rm(list=ls())



# workflow ---------------------------------------------------------------------------

# use trueIsolation pkg, which exports zoom boxes, to create graph maps that focus on
# particular areas of interest

trueIsolation::zoom.boxes


# ------------------------------------------------------------------------------

smpl <- trueIsolation::zoom.boxes[7,]

?setup.gh.wrapper

bbx <- smpl$bbx[[1]]
czsf <- geox::build.CZs(smpl$cz)

devtools::load_all()

trimmed.gh <-
  setup.gh.wrapper(
    cz = smpl$cz
    #sfx = bbx
    ,directed = F
    ,tracts.or.groups = 'bg'
    ,min.flows = 1
    ,tie.str.drop.deciles = 5
  )



trimmed.gh <- st_transform(trimmed.gh, 4326)
lonlats <- gh2coords(trimmed.gh)

sttm <- visaux::get.stamen.bkg(sfx = czsf # bbx
                               ,zoom = 11)
ggmap(sttm)

library(ggraph)
ggmap(sttm,
      base_layer =
        ggraph(trimmed.gh
               ,layout = lonlats)) +
  geom_edge_fan(aes(edge_alpha = tstr
                    ,edge_width = tstr)
                ,color = '#008080') +
  flow.map.base() +
  visaux::bbox2ggcrop(smpl$bbx[[1]])

# layer over a demographics chloropleth?
demos <- sfg.seg::demos.2019

# demos <- demos %>% sfg.seg::cbg.demos2tracts()

demosf <- demos %>%
  geox::geosubset(cz = smpl$cz) %>%
  geox::attach.geos(query.fcn = tigris::block_groups
                    ,year = 2019)
demosf <- demosf %>% st_transform(4326)

ggmap(sttm,
      base_layer =
        ggraph(trimmed.gh
               ,layout = lonlats)) +
  geom_sf(data = demosf
          ,aes(fill = perc.wh)
          ,color = NA
          ,alpha = .3) +
  geom_edge_fan(aes(edge_alpha = tstr
                    ,edge_width = tstr)
                ,color = '#008080') +
  flow.map.base() +
  scale_fill_viridis(option = 'C'
                     ,limits = c(0,1)) +
  visaux::bbox2ggcrop(smpl$bbx[[1]]) +
  labs(caption = smpl$an)



# wrapper function -------------------------------------------------------------

graph.map.wrapper <- function(area
                              #,save.dir =
                              ,min.flows = 10
                              ,sttm.zoom = 11
                              ,tie.str.drop.deciles = 5
                              ) {


  bbx <- area$bbx[[1]]
  czsf <- geox::build.CZs(area$cz)

  devtools::load_all()
  trimmed.gh <-
    setup.gh.wrapper(
      cz = area$cz
      #sfx = bbx
      ,directed = F
      ,tracts.or.groups = 'bg'
      ,min.flows = min.flows
      ,tie.str.drop.deciles = tie.str.drop.deciles
    )



  trimmed.gh <- st_transform(trimmed.gh, 4326)
  lonlats <- gh2coords(trimmed.gh)

  sttm <- visaux::get.stamen.bkg(sfx = czsf # bbx
                                 ,zoom = 11)

  require(ggraph)

  # layer over a demographics chloropleth?
  demos <- sfg.seg::demos.2019

  # demos <- demos %>% sfg.seg::cbg.demos2tracts()

  demosf <- demos %>%
    geox::geosubset(cz = area$cz) %>%
    geox::attach.geos(query.fcn = tigris::block_groups
                      ,year = 2019)
  demosf <- demosf %>% st_transform(4326)

  ggmap(sttm,
        base_layer =
          ggraph(trimmed.gh
                 ,layout = lonlats)) +
    geom_sf(data = demosf
            ,aes(fill = perc.wh)
            ,color = NA
            ,alpha = .3) +
    geom_edge_fan(aes(edge_alpha = tstr
                      ,edge_width = tstr)
                  ,color = '#008080') +
    flow.map.base() +
    scale_fill_viridis(option = 'C'
                       ,limits = c(0,1)) +
    visaux::bbox2ggcrop(bbx) +
    labs(caption = area$an)


    visaux::ragg.wrapper(fn = area$an)

}

# graph.map.wrapper(smpl)

# iterate through zoom boxes ---------------------------------------------------

map(1:nrow(trueIsolation::zoom.boxes)
    , ~graph.map.wrapper(trueIsolation::zoom.boxes[.,])
    )



# variations -------------------------------------------------------------------

# brooklyn had too many flows
bk <- trueIsolation::zoom.boxes %>%
  filter(an == 'brooklyn')

visaux::get.stamen.bkg(sfx = bk$bbx[[1]]
                         ,zoom = 12) %>% ggmap()


graph.map.wrapper(area = bk,
                  min.flows = 15
                  ,tie.str.drop.deciles = 6
                  ,sttm.zoom = 11)










# directed versions -------------------------------------------------------------

# (not finished, but could be interesting to play with)

bbx <- smpl$bbx[[1]]
czsf <- geox::build.CZs(smpl$cz)

devtools::load_all()
trimmed.gh <-
  setup.gh.wrapper(
    cz = smpl$cz
    #sfx = bbx
    ,directed = T
    ,tracts.or.groups = 'bg'
    ,min.flows = 1
    ,tie.str.drop.deciles = 5
  )



trimmed.gh <- st_transform(trimmed.gh, 4326)
lonlats <- gh2coords(trimmed.gh)

sttm <- visaux::get.stamen.bkg(sfx = czsf # bbx
                               ,zoom = 11)
ggmap(sttm)

library(ggraph)
ggmap(sttm,
      base_layer =
        ggraph(trimmed.gh
               ,layout = lonlats)) +
  geom_edge_fan(aes(edge_alpha = tstr
                    ,edge_width = tstr)
  ) +
  flow.map.base() +
  visaux::bbox2ggcrop(smpl$bbx[[1]])



