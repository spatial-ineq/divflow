#' My flow maps came a long way (as with .mapping big OD notes and hh-income-flow
#' maps). But those are both done kind of manually. This script workshopped a
#' workflow to send to Della or to use to crop from a given Place/centroid
#'
#' My strategy starting out will be:
#'  start with given sf object -->
#'   buffer around it x distance -->
#'   crop to that bbox -->
#'   take all tracts/CBGs within the bbox -->
#'    (this requires interfacing to CZs to load sfg data)
#'   then map those flows!


# ws ---------------------------------------------------------------------------
rm(list = ls())
require(sf)
require(tidyverse)

# option setting
sf_use_s2(F)
options(tigris_use_cache = TRUE)

# sample area ------------------------------------------------------------------

stl <- geox::rx %>%
  filter(grepl("St. Louis", cz_name))
# I think flexibly generating the crop area from any sf is the way, but I also
# thinking making this mappable across Places is great, so I'm goanna start with the
# StL place

plc <-
  divM::largest.plc.in.cz %>%
  filter(cz.id %in% stl$cz)

# setting parameters -----------------------------------------------------------

centroid.geom <- plc
map.buffer <- units::set_units(5
                               , 'miles')
max.dst <- units::set_units(5
                            , 'miles')
min.flows <- 10

# thinking about buffer areas --------------------------------------------------
library(mapview)

buffered.plc <- plc[1] %>%
  divM::conic.transform() %>%
  st_buffer(map.buffer)
mapview(plc) + mapview(buffered.plc)

# I think that shows me that plc centroid will actually be better.
buffered.plc <- plc[1] %>%
  divM::conic.transform() %>%
  st_centroid() %>%
  st_buffer(map.buffer)

mapview(plc) + mapview(buffered.plc)
# that looks great!

# get overlapping CZs ----------------------------------------------------------

# actually have to load within map buffer area + max distance, to include flows from
# nearby cropped neighborhoods
buffered.plc <- plc[1] %>%
  divM::conic.transform() %>%
  st_centroid() %>%
  st_buffer(map.buffer + max.dst)

cos <- map_dfr(unique(stl$statefp)
               ,~tigris::counties(.x,year = 2019))
cos <- cos %>% st_transform(st_crs(buffered.plc))
ovcos <- st_crop(cos
                 ,buffered.plc)

# check:
mapview(ovcos ) +
  mapview(buffered.plc)

czs2load <-
  geox::rx %>%
    filter(countyfp %in% ovcos$GEOID)

czs2load <- unique(czs2load$cz)

sfg <- sfg.seg::read.sfg.CZs(czs2load
                             ,sfg.dir =
                               paste0( Sys.getenv('drop_dir')
                               ,'sfg-processed/orig_dest_annual/'))

# I forget if the cz-levels are trimmed by origin or both origin AND dest.
# So I make sure not all dests are in loaded cz (they're not)
sfg

#sfg_backup <- sfg
sfg <- sfg_backup

# filter sfg to study area -----------------------------------------------------

sfg <- sfg %>%
  geox::geosubset(c('origin', 'dest') # can use this but only just fixed
                  ,cz = czs2load)
  #sfg.seg::geo.subset.cbgs(cz = czs2load)
# a check:
dests <- sfg$dest %>% unique() %>% tibble(geoid = .) %>%
  geox::attach.geos(query_fcn = tigris::block_groups)


# now filter to CROP+dst area
ellcbs <- c(sfg$origin, sfg$dest) %>% unique()
full_cbsf <- tibble(geoid = ellcbs) %>%
  geox::attach.geos(query_fcn = tigris::block_groups) %>%
  st_transform(st_crs(buffered.plc))
cbsf <- st_crop(full_cbsf,
                buffered.plc)

cbsf %>% plot() # love that :)

# get centroids
ctrs  <- cbsf %>% st_centroid()

# filter to w/in area ( although other distance filter later too)
sfg <- filter(sfg,
              origin %in% cbsf$geoid |
                dest %in% cbsf$geoid )

# transform to BGs (if you want) -----------------------------------------------
sfg.seg::cbg.flows2tracts

# drop loops, get tie strength, and filter -------------------------------------
sfg <- sfg %>% filter(origin != dest)


# get tie proportions (leaving directed for now)
tsf <- sfg %>%
  group_by(origin) %>%
  mutate(total.from = sum(n)) %>%
  group_by(dest) %>%
  mutate(total.inc = sum(n)) %>%
  ungroup()

tsf <- tsf %>%
  mutate( perc.to.dest = n / total.inc
          ,perc.from.origin = n / total.from)

# preliminary trim by just flows
tsf$n %>% quantile(seq(0,1,.1))
tsf <- tsf %>% filter(n >= min.flows)

# make gh object ---------------------------------------------------------------

library(igraph)
library(tidygraph)
gh <- tidygraph::tbl_graph(
  edges = tsf,
  directed = T
)

#' inner join to filter out of crop distance
#' (although we'll also filter by distance later)
gh <- gh %>%
  activate('nodes') %>%
  inner_join(ctrs
            ,by = c('name' = 'geoid'))


# spatialize edges
gh <- gh %>%
  activate('edges') %>%
  sfnetworks::as_sfnetwork(
    .
    ,directed = T
    ,edges_as_lines = T)

# get lengths
gh <- gh %>%
  activate('edges') %>%
  mutate(dst = st_length(geometry))

# descriptives & filtering ------------------------------------------------------------

sfe <- gh %>% divseg::get_edges()
ght <- gh
# distribution of edge attributes
tibble(sfe) %>% select(n,
               dst,
               matches('^perc')) %>%
  map( summary )

# filter by distance
ght <- ght %>%
  activate('edges') %>%
  filter(dst <= max.dst)

# filter by tie strength
ght <- ght %>%
  activate('edges') %>%
  filter_at( vars(matches('^perc'))
             ,any_vars(. > .01))

# descriptives after trim
ght
# distribution of edge attributes
ght %>% divseg::get_edges() %>%
  tibble() %>%
  select(n,
         dst,
         matches('^perc')) %>%
  map( summary )

# last spatial set ups for plot ------------------------------------------------

# get node lonlats
lonlats <- gh %>%
  divseg::get_nodes() %>%
  st_sf() %>%
  st_coordinates()

# get bbox for crop
zoom.box <- plc[1] %>%
  divM::conic.transform() %>%
  st_centroid() %>%
  st_buffer(map.buffer + max.dst
            ) %>%
  st_bbox()

ghc <- gh %>%
  activate('nodes') %>%
  st_crop(zoom.box)
lonlats <- ghc %>%
  divseg::get_nodes() %>%
  st_sf() %>%
  st_coordinates()

# plot! ------------------------------------------------------------------------

library(ggraph)
tmplot <- ght %>%
  activate('edges') %>%
  filter(n >
           100) %>%
  #100) %>%  # big filter to workshop
  mutate(ttstr =
           #log(1 + tie.strength)
         n
         ) %>%
  arrange(ttstr) %>%
  ggraph(layout = lonlats  ) +
  geom_edge_density(

    aes(edge_fill = ttstr)
                    ) +
  geom_edge_fan(aes(edge_alpha =
                      ttstr
                    #,edge_color = ttstr
                    ,edge_width = ttstr)) +
  #scale_edge_color_viridis() +
  #scale_edge_fill_viridis() +
  scale_edge_alpha_continuous(guide = 'none'
                              ,range = c(.05,1)) +
  scale_edge_width_continuous(guide = 'none'
                              ,range = c(.05,2)) +
  scale_color_discrete(guide = 'none')

tmplot +
  coord_sf( xlim =
              c( zoom.box$xmin
                 ,zoom.box$xmax)
            , ylim =
              c(zoom.box$ymin
                ,zoom.box$ymax)
            ,expand = F
  )



# other layers -----------------------------------------------------------------

# (i think on top of density but under edges is ideal..)
geom_sf(data = st_crop(wtr
                       ,zoom.box)
        , color = NA, # feature layers
        fill = '#94bdff') +
  geom_sf(data = st_crop(hwys # feature layers
                         ,zoom.box)
          , aes(color = SIGNT1)
          , size = 1.5) +
  geom_sf(data = st_crop(st_boundary(plcs),
                         zoom.box)
          , size = .5
          ,color = '#583799')



  # more basic edge plot ---------------------------------------------------------

ght %>%
  divseg::get_edges() %>%
  filter(n >
           100) %>%
  #rename(nf = n)
  ggplot() +
  geom_sf(aes( size = n
              ,alpha = n))+
  scale_color_viridis_c() +
  scale_alpha_continuous(guide = 'none'
                              ,range = c(.05,1)) +
  scale_size_continuous(guide = 'none'
                              ,range = c(.05,2.5)) +
  theme_void() %>%
  coord_sf( xlim =
              c( zoom.box$xmin
                 ,zoom.box$xmax)
            , ylim =
              c(zoom.box$ymin
                ,zoom.box$ymax)
            ,expand = F
  )




# testing functionalization ----------------------------------------------------

devtools::document()
devtools::load_all()

map.buffer <- units::set_units(5   , 'miles')
max.dst <- units::set_units(5 , 'miles')
min.flows <- 20

plc <- divM::largest.plc.in.cz %>%
  filter(cz.id %in% '24701')

ctr <- plc %>%
  divM::conic.transform() %>%
  #st_transform(crs) %>%
  st_centroid()

eligible.area <- ctr %>%
  st_buffer(map.buffer + max.dst)

bounds.area <- ctr %>%
  st_buffer(map.buffer)

tmp.sfg <- sfx2sfg(eligible.area
                   ,min.flows = 20)

bg.gh <- sfg2gh(tmp.sfg
                 ,'bg'
                 ,directed= T) #undirected ct graph
bg.gh <- spatialize.graph(bg.gh, eligible.area, 'bg', T)
bg.gh <- bg.gh %>%
  apply.flow.filters(
    bounds.sf = bounds.area
    ,directed = T
    #,min.str = 0
    )

bg.gh

bg.gh %>%
  activate('edges') %>%
  as_tibble() %>%
  map( summary )


devtools::load_all()
library(units)
gh <- setup.gh.wrapper(
   sfx = plc
  ,map.buffer = set_units(6, 'miles')
  ,max.dst = set_units(3, 'miles')
  ,directed = F
  ,tracts.or.groups = 'ct'
)

div.lyrs <- get.divlyrs

flow.map.wrapper(sfx = plc
                 ,edge.attr = 'perc.from.origin'
                 ,map.buffer = set_units(6, 'miles')
                 ,max.dst = set_units(3, 'miles')
                 ,directed = F
                 ,tracts.or.groups = 'ct')



# fixing undirected tie.str fcn ------------------------------------------------

ugh <- sfg2gh(tmp.sfg
                ,'ct'
                ,directed= F) #undirected ct graph

# as from
od <- tmp.sfg %>% filter(n >= 30)
devtools::load_all()
ugh <-
  od %>%
  sfg2gh(
    #tracts.or.groups = 'ct',
    directed = F)

dgh <-
  od %>%
  sfg2gh(
    #tracts.or.groups = 'ct',
    directed = T)

ugh
dgh


od %>%
  filter(dest == '171336001011' )
library(igraph)

tidygraph::as_tbl_graph(
  od
  ,directed = F
  ,mode = "collapse"
  ,edge.attr.comb = "sum")

dgh <- igraph::graph.data.frame(od, directed = T)
ugh <- igraph::as.undirected(dgh,
                             mode = "collapse",
                             edge.attr.comb = "sum")
ugh <-  tidygraph::as_tbl_graph(ugh)

ugh %>%
  activate('edges') %>%
  group_by(from) %>%
  mutate(total.from = sum(n)) %>%
  group_by(to) %>%
  mutate(total.inc = sum(n)) %>%
  ungroup()


ugh %>%
  activate('edges') %>%
  filter(to == 2 |
           from == 2)

#get.normaized.undirected.connectedness <- function(ugh, flow.colm = 'n') {

  # just work with edge list el
  el <- ugh %>%
    activate('edges') %>%
    as_tibble()

  # subset to ods containing either of the two tracts
  total.trips <-
    el[with(el,
            from %in% c(x, y) |
              to %in% c(x, y) ), ]


}


# thinking about how to spatial trim -------------------------------------------

# i kept a large buffer around the actual area originally, but then that makes
# trimming so difficult.

# okay this is good. cropping edges creates a new artificial node (with no name) on
# border of cropping area. Edges that are entirely outside of area are cropped.
# a lot of new nodes are added b/c aggregate effect of artifical ones added in
tmp1 <- bg.gh %>%
  #activate('nodes') %>%
  activate('edges') %>%
  st_crop(map.area) %>%
  apply.flow.filters(directed = T)

nas.after.crop <- tmp1 %>%
  activate('nodes') %>%
  filter(is.na(name))


sfe <- nas.after.crop %>%
  activate('edges') %>%
  st_as_sf()

sfe$dst %>% as.numeric() %>%  summary()


sfn <- nas.after.crop %>%
  activate('nodes') %>%
  st_as_sf()

sfe %>%
  ggplot()+
  geom_sf(aes( size = log(n)
              ,alpha = log(n))) +
  scale_alpha_continuous(guide = 'none'
                        ,range = c(.5, 1)) +
  scale_size_continuous(guide = 'none'
                        ,range = c(.5, 2)) +
  theme_void()


sfn %>%
  mutate(id = row_number()) %>%
  mapview()

lonlats <- nas.after.crop %>%
  activate('nodes') %>%
  st_coordinates()

nas.after.crop %>%
  ggraph::ggraph(layout = lonlats ) +
  geom_edge_fan(aes(edge_alpha =
                      n
                    ,edge_width =
                      n)) +
  scale_edge_alpha_continuous(guide = 'none'
                              ,range = c(.05, 1)) +
  scale_edge_width_continuous(guide = 'none'
                              ,range = c(.05, 2)) +
#  scale_color_discrete(guide = 'none') +
  theme_void()

flow.map <- flow.map.wrapper(map.area,
                 gh = bg.gh,
                 edge.attr =
                   'perc.to.dest'
                 )
flow.map


zoom.box <- map.area %>% st_bbox()



bg.ghT <-  bg.gh %>%
  activate('nodes') %>%
  st_crop(zoom.box)

lonlats <- bg.gh %>%
  activate('nodes') %>%
  st_as_sf() %>%
  st_coordinates()

bg.gh %>%
  activate('edges') %>%
  filter(n > 30) %>%
  mutate(ex =
           log(perc.to.dest)) %>%
  ggraph::ggraph(layout = lonlats ) +
  #geom_edge_density(
  #  aes(edge_fill = ex)
  #) +
  #div.layers +
  geom_edge_fan(aes(edge_alpha =
                      ex
                    ,edge_width =
                      ex)) +
  scale_edge_alpha_continuous(guide = 'none'
                              ,range = c(.1, 1)) +
  scale_edge_width_continuous(guide = 'none'
                             ,range = c(.1, 2)) +
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



# reminder to fix this tiestr fcn... (hardcoded colnames, clunky)
tmp.sfg %>%
  head(1000) %>%
  mutate(tstr =
           map2_dbl(origin, dest,
                    ~get.normalized.undirected.connectedness(.x, .y,
                                                             el =
                                                               head(tmp.sfg
                                                                    ,1000)))
  )
divseg::

# how to deal with MULTIPOLYGON bgs? -------------------------------------------

# ah, they were created from cropping
sfn <- gh %>%
  divseg::get_nodes()
# i think dropping them is appropriate, b/c they're on edge of cropped area
# or i just let st_centroid handle them

tmp <- sfn %>%
  filter(st_geometry_type(geometry) %in%
           'MULTIPOLYGON') %>%
  st_sf()

# illustrates why i feel comfy dropping
full_cbsf %>%
  filter(geoid %in%
           tmp$name) %>%
  ggplot() +
  geom_sf() +
  geom_sf(data = st_boundary(buffered.plc)
          ,color = 'red')

mapview(tmp) +
  mapview(st_buffer(st_centroid(tmp), 100), color = 'red'
          ,cex = 5)

lonlats <- gh %>%
  divseg::get_nodes() %>%
  st_sf() %>%
  st_coordinates()
tmp


# thinking about reorganizing tie str fcn --------------------------------------
?divseg::get.normalized.undirected.connectedness


# this does something, leaving directionality, which I kind of like.
tsf <- sfg %>%
  group_by(origin) %>%
  mutate(total.from = sum(n)) %>%
  group_by(dest) %>%
  mutate(total.inc = sum(n)) %>%
  ungroup()

tsf <- tsf %>%
  mutate( perc.to.dest = n / total.inc
         ,perc.from.origin = n / total.from)

tsf %>%
  group_by(dest) %>%
  summarise(sum(perc.to.dest))
tsf %>%
  mutate(tstr =

           )

# scratch ----------------------------------------------------------------------


ctsf <- stl$countyfp %>%
  map_dfr(
    ~tigris::tracts(state = substr(.x, 1,2)
                    ,county = substr(.x, 3,5)
                    ,year= 2019)
  )

ctsf
