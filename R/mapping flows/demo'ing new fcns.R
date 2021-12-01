
# ws ---------------------------------------------------------------------------
rm(list = ls())
require(sf)
require(tidyverse)

# option setting
sf_use_s2(F)
options(tigris_use_cache = TRUE)

# sample area ------------------------------------------------------------------

sample.area <- geox::rx %>%
  filter(grepl("St. Louis", cz_name))

devtools::document()
devtools::load_all()

# build relevant distances (not necessary but helps run faster)
plc <- divM::largest.plc.in.cz %>%
  filter(cz.id %in%
           unique(sample.area$cz)
           )


map.range <- units::set_units(6 , 'miles')
max.distance <- units::set_units(10 , 'miles')

ctr <- plc %>% divM::conic.transform() %>% st_centroid()

eligible.area <- ctr %>%
  st_buffer(map.range + max.distance)

frame.area <- ctr %>%
  st_buffer(map.range)

sttm <- visaux::get.stamen.bkg(
  frame.area,
  zoom = 11
)

ggmap(sttm)


# step by step gh setup --------------------------------------------------------

devtools::load_all()
od <-
  sfx2sfg(
    sfx = eligible.area
    ,min.flows = 5
    ,tracts.or.groups = 'bg' # ct
    ,base.dir = Sys.getenv('drop_dir')
    ,sfg.dir = 'sfg-processed/orig_dest_annual/'
    ,year= 2019
    ,trim.loops = T)


gh <- sfg2gh(
  od
  ,directed = F
  )


ghsf <- spatialize.graph(
  gh
  ,frame.sf = frame.area
  ,tracts.or.groups = 'ct'
  ,directed = F
)

fgh <- apply.flow.filters(
  ghsf
  ,min.flows = 5
  , tie.str.deciles = 7
  ,frame.sf = frame.area
  ,max.dst = max.distance
)
fgh

# consider alternate flow trims ------------------------------------------------
full.gh
ghsf
fgh
ghsf %>% filter(dst > max.distance)
ghsf %>% pull(tstr) %>% quantile(seq(0,1,.1))
ghsf %>% filter(tstr > .001)

fgh %>% filter(tstr > .01)

# map stl --------------------------------------------------------------------------
to.map <- fgh
#to.map <- fgh %>% activate('nodes') %>% st_crop(frame.area)

fgh
lonlats <- to.map %>%
  activate('nodes') %>%
  as_tibble() %>%
  st_sf() %>%
  st_coordinates()

ggflbase <- flow.map.base()

devtools::load_all()

stl.flow.map <- to.map %>%
  ggraph(layout = lonlats ) +
#  geom_edge_density(
#    aes(edge_fill = n) ) +
  geom_edge_fan(aes(edge_alpha =
                      tstr
                    ,edge_width =
                      tstr
                    )) +
  scale_edge_alpha_continuous(guide = 'none'
                              ,range = c(.1,1)) +
  scale_edge_width_continuous(guide = 'none'
                             ,range = c(.1, 1.5)) +
  scale_color_discrete(guide = 'none')

stl.flow.map + ggtitle('stl flow map')

stl.flow.map + sfx2coord_sf(frame.area)



# BGs --------------------------------------------------------------------------


# step by step gh setup --------------------------------------------------------

devtools::load_all()
od <-
  sfx2sfg(
    eligible.sf = eligible.area
    ,min.flows = 10
    ,tracts.or.groups = 'bg'
    ,base.dir = Sys.getenv('drop_dir')
    ,sfg.dir = 'sfg-processed/orig_dest_annual/'
    ,year= 2019
    ,trim.loops = T)


full.gh <- sfg2gh(
  od
  ,directed = F
)

devtools::load_all()
ghsf <- spatialize.graph(
  full.gh
  ,frame.sf = frame.area
  ,tracts.or.groups = 'bg'
  ,directed = F
)

fgh <- apply.flow.filters(
  ghsf
  ,min.flows = 5
  , tie.str.deciles = 8
  ,frame.sf = frame.area
  ,max.dst = max.distance
)
ghsf
fgh

# consider alternate flow trims ------------------------------------------------
full.gh
ghsf
fgh

fgh %>% filter(tstr > .01)

# map stl --------------------------------------------------------------------------
to.map <- fgh
#to.map <- fgh %>% activate('nodes') %>% st_crop(frame.area)

fgh
lonlats <- to.map %>%
  activate('nodes') %>%
  as_tibble() %>%
  st_sf() %>%
  st_coordinates()

ggflbase <- flow.map.base()

devtools::load_all()

stl.flow.map <- to.map %>%
  ggraph(layout = lonlats ) +
    geom_edge_density(
      aes(edge_fill = tstr) ) +
  geom_edge_fan(aes(edge_alpha =
                      tstr
                    ,edge_width =
                      tstr
  )) +
  scale_edge_alpha_continuous(guide = 'none'
                              ,range = c(.2, 1)) +
  scale_edge_width_continuous(guide = 'none'
                              ,range = c(.2, 2)) +
  scale_color_discrete(guide = 'none')

#stl.flow.map + ggtitle('stl flow map')

stl.flow.map +
  sfx2coord_sf(frame.area) +
  divlyrs


# experiments ------------------------------------------------------------------


# Using ragg -------------------------------------------------------------------
library(ragg)
#sdir <- 'R/mapping flows/ragg/'
f <- knitr::fig_path('.png')
ragg::agg_png(f
              ,width = 8
              ,height = 6
              ,res = 300
              ,units = 'in'
              ,scaling = 1)
stl.flow.map +
  sfx2coord_sf(frame.area) +
  divlyrs
invisible(dev.off())
# that's cool and fast and seems like good standardisation

# stamen or google basemap -----------------------------------------------------
library(ggmap)
?register_google

lon.lat.ctr <-
  ctr %>%
  st_transform(4326) %>%
  st_coordinates()

lon.lat.ctr %>%
  as.vector() %>%
  setNames(c('lon', 'lat'))

?ggmap::get_googlemap()
ggm <- ggmap::get_googlemap(
  center = c(lon.lat.ctr %>%
               as.vector() %>%
               setNames(c('lon', 'lat')))
)

bbx <- st_bbox(st_transform(frame.area
                     ,4326))

sttm <- ggmap::get_stamenmap(
  bbox = c(left = bbx[['xmin']],
           bottom = bbx[['ymin']],
           right = bbx[['xmax']],
           top = bbx[['ymax']])
  ,zoom = 12
  ,maptype = 'toner-background'
  ,crop = T
)
?coord_sf
ggmap(ggm) +
  sfx2coord_sf(st_transform(frame.area
                            ,4326))

# reproject and re retrieve lonlats
to.map <- to.map %>% st_transform(4326)
lonlats <- to.map %>%
  activate('nodes') %>%
  as_tibble() %>%
  st_sf() %>%
  st_coordinates()

stamen.stl <-
  ggmap(sttm
      ,base_layer =
        ggraph(to.map, layout = lonlats ) ) +
  #sfx2coord_sf(st_transform(frame.area
  #                          ,4326)) +
  #geom_edge_density(
  #  aes(edge_fill = tstr) ) +
  geom_edge_fan(aes(edge_alpha =
                      log(tstr)
                    ,edge_width =
                      log(tstr)
                    #,edge_color =
                    #  log(tstr)
  ),color = '#007c91'
  ) +
  scale_edge_alpha_continuous(guide = 'none'
                              ,range = c(.1, 1)) +
  scale_edge_width_continuous(guide = 'none'
                              ,range = c(.1, 1.5)) +
  #scale_edge_color_viridis(option = 'A') +
  scale_color_discrete(guide = 'none') +
  theme_void()
#stamen.stl
f <- knitr::fig_path('.png'
                     ,number =2)
ragg::agg_png(f
              ,width = 8
              ,height = 6
              ,res = 300
              ,units = 'in'
              ,scaling = 1)
stamen.stl
invisible(dev.off())


# variants using ragg and stamen -----------------------------------------------

# different filters
fgh <- apply.flow.filters(
  ghsf
  ,min.flows = 5
  , tie.str.deciles = 8
  ,frame.sf = frame.area
  ,max.dst = max.distance
)
fgh
to.map <- fgh

# reproject and re retrieve lonlats
to.map <- to.map %>% st_transform(4326)
lonlats <- to.map %>%
  activate('nodes') %>%
  as_tibble() %>%
  st_sf() %>%
  st_coordinates()

sttm %>%
  divM::conic.transform()


stamen.stl <-
  ggmap(sttm
        ,base_layer =
          ggraph(to.map, layout = lonlats ) ) +
  #sfx2coord_sf(st_transform(frame.area
  #                          ,4326)) +
  #geom_edge_density(
  #  aes(edge_fill = tstr) ) +
  geom_edge_fan(aes(edge_alpha =
                      tstr
                    ,edge_width =
                      tstr
                    #,edge_color =
                    #  log(tstr)
  ),color = '#007c91'
  ) +
  scale_edge_alpha_continuous(guide = 'none'
                              ,range = c(.5, 1)) +
  scale_edge_width_continuous(guide = 'none'
                              ,range = c(.1, 1.2)) +
  #scale_edge_color_viridis(option = 'A') +
  scale_color_discrete(guide = 'none') +
  theme_void()
#stamen.stl
f <- knitr::fig_path('.png'
                     ,number =3)
ragg::agg_png(f
              ,width = 4
              ,height = 3
              ,res = 550
              ,units = 'in'
              ,scaling = 1)
stamen.stl
invisible(dev.off())

# varying start plcs -----------------------------------------------------------

plcs <- tigris::places(state = 29)
plcs %>% mapview::mapview()
zoom.ctr <- plcs %>% filter(NAME == 'Wellston')

# using wrapper ----------------------------------------------------------------

devtools::load_all()

# undirected cts
ugh <- setup.gh.wrapper(
  zoom.ctr
  ,map.buffer =  units::set_units(7 , 'miles')
  ,max.dst = units::set_units(8, 'miles')
  ,min.flows = 10
  ,min.str = 0 # NULL
  ,directed = F
  ,tracts.or.groups = 'bg' #'ct'#c('ct', 'bg')
  ,year = 2019
  ,base.dir = Sys.getenv('drop_dir')
  ,sfg.dir = 'sfg-processed/orig_dest_annual/'
  ,crs = '+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45'
)

# directed cts
dgh <- setup.gh.wrapper(
  zoom.ctr
  ,map.buffer =  units::set_units(7 , 'miles')
  ,max.dst = units::set_units(8, 'miles')
  ,min.flows = 10
  ,min.str =  0#NULL
  ,directed = T
  ,tracts.or.groups = 'ct'#c('ct', 'bg')
  ,year = 2019
  ,base.dir = Sys.getenv('drop_dir')
  ,sfg.dir = 'sfg-processed/orig_dest_annual/'
  ,crs = '+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45'
)



# mapping ----------------------------------------------------------------------

tmp <- fgh
tmp <- ugh

# get bbox for crop
zoom.box <- st_bbox(gh2nodes(tmp))

ggbase <- flow.map.base()

# divlyrs <- get.div.layers(gh2nodes(tmp))
# ggplot() + divlyrs # <3


# trim to x edges
ttmp <- tmp %>%
  activate('edges') %>%
  arrange(desc(tstr)) %>%
  mutate(tstr.rank = row_number()) %>%
  filter(tstr.rank <= 2000)

require(ggraph)

flow.map <-
  ttmp %>%
  mutate(ex = get('tstr')) %>%
  mutate(ex = log(ex)) %>%
  ggraph::ggraph(layout = st_coordinates(gh2nodes(.)) ) +
  #geom_edge_density(  ) +
  divlyrs +
  geom_edge_fan(aes(edge_alpha =
                      ex
                    ,edge_width =
                      ex)) +
  scale_edge_alpha_continuous(guide = 'none'
                              ,range = c(.05,1)) +
  scale_edge_width_continuous(guide = 'none'
                              ,range = c(.05,2)) +
  scale_color_discrete(guide = 'none') +
  theme_void() +
  coord_sf( xlim =
              c( zoom.box$xmin
                 ,zoom.box$xmax)
            , ylim =
              c(zoom.box$ymin
                ,zoom.box$ymax)
            ,expand = F
  )

flow.map

# vis --------------------------------------------------------------------------
devtools::load_all()

# undirected
plc <- plc %>% divM::conic.transform()
ug.flo.map <- flow.map.wrapper(sfx = plc
                 ,gh = ugh
                 ,map.buffer = map.range
                 ,max.dst = max.distance
                 ,edge.attr = 'tstr')
ug.flo.map + divlyrs

# diagnostics ------------------------------------------------------------------
dgh




# scratch ----------------------------------------------------------------------


map.buffer <- units::set_units(6 , 'miles')
max.distance <- units::set_units(3 , 'miles')
eligible.area <- ctr %>%
  st_buffer(map.buffer + max.distance)

frame.area <- ctr %>%
  st_buffer(map.buffer)
