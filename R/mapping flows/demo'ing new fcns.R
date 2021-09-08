
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

eligible.area <- ctr %>%
  st_buffer(map.range + max.distance)

frame.area <- ctr %>%
  st_buffer(map.range)

# div layers
divlyrs <- get.div.layers(bounds.sf = frame.area)
ggplot() + divlyrs # magic!


# step by step gh setup --------------------------------------------------------
devtools::load_all()
od <-
  sfx2sfg(
    eligible.sf = eligible.area
    ,min.flows = 30
    ,tracts.or.groups = 'ct'
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
  ,tracts.or.groups = 'ct'
  ,directed = F
)

fgh <- apply.flow.filters(
  ghsf
  ,frame.sf = frame.area
  ,max.dst = max.distance
)

fgh



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
