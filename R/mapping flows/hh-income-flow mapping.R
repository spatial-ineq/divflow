# ws ----------------------------------------------------------------------

library(tidyverse)

# option setting
library(sf)
sf_use_s2(F)
options(tigris_use_cache = TRUE)
Sys.setenv("VROOM_SHOW_PROGRESS"="false")

devtools::load_all()
#library(sfg.seg)

rm(list=ls())


# build demos table -------------------------------------------------------

# preppd locally (with divseg)
res <- vroom::vroom('.quick exports and viz/res-vars.csv') %>% select(-`...1`)

# denver cbsa

# get a graph ------------------------------------------------------------------

library(tidygraph)

gh <-
  divseg::setup.cts.as.network(cbsa.id = '19740'
                               ,directed =T
                               ,flow.floor = 10
                               ,add.centroid.geos = F
                               ,turn2tracts = T
                               ,round.flows = F)

gh

# add other demos
gh <- gh %>%
  activate('nodes') %>%
  left_join(res
            ,by = c('name' = 'geoid19'))

gh <- gh %>%
  activate('nodes') %>%
  mutate(prop.fam.75p =
           cut(n.fam75p / hh,
               seq(0,1,.25)))


# add geos
gh <- gh %>%
  spatialize.nodes()
gh <- gh %>%
  activate('edges') %>%
  sfnetworks::as_sfnetwork(., directed = T,
                           edges_as_lines = T)
gh <- gh %>%
  activate('edges') %>%
  mutate(len = # in km
           as.numeric(st_length(geometry)) / 1000)




# filters ----------------------------------------------------------------------


# trim..
full.gh <- gh
# gh <- full.gh

# ..by distance
gh %>% activate('edges') %>% pull(len) %>% summary()
# 40 mile (64 km) max
gh <- gh %>%
  activate('edges') %>%
  filter(len <= 64)

# ..by trips
gh <- gh %>%
  activate('edges') %>%
  filter(n >= 15)

el  <- divseg::get_edges(gh)
"gh <- gh %>%
  activate('edges') %>%
  mutate(tie.strength =
           map2_dbl(el$from, el$to,
                    ~get.normalized.undirected.connectedness(.x, .y,
                                                             el = el))
  )
# ..by tie
gh <-
  gh %>%
  activate('edges') %>%
  filter(tie.strength >
           quantile(tie.strength,
                    seq(0,1,.25))[2]) # drop n quartiles
"

sfn <- divseg::get_nodes(gh)
sfe <- divseg::get_edges(gh)



# map --------------------------------------------------------------------------

sfn %>% mapview::mapview(zcol = 'prop.fam.75p')
sfe$n %>% summary()
# get coords
lonlats <- sfn %>% st_coordinates()

gh %>%
  activate('nodes') %>%
  filter(is.na(prop.fam.75p))

sfe$tie.strength %>% quantile(seq(0,1,.1))

library(ggraph)

sfn <- ght %>% divseg::get_nodes()
sfe <- ght %>% divseg::get_edges()

# get coords
lonlats <- sfn %>% st_coordinates()


gh %>%
  activate('edges') %>%
  filter(tie.strength > .005) %>%
  mutate(ttstr =
           log(1 + tie.strength)) %>%
  mutate(to_prop.fam.75p =
           .N()$prop.fam.75p[to]) %>%
  ggraph(layout = lonlats  ) +
  #geom_sf(data = wtr, color = NA, # feature layers
  #        fill = '#94bdff') +
  #geom_sf(data = hwys, aes(color = SIGNT1)) + # feature layers
  #geom_edge_density(aes(edge_fill = ttstr)) +
  geom_edge_fan(aes(edge_alpha =
                      ttstr#n
                    ,edge_color =
                      to_prop.fam.75p
                    ,edge_width =
                      ttstr#n
                    )) +
  #scale_edge_color_viridis() +
  scale_edge_color_viridis(discrete = T) +
  scale_edge_alpha_continuous(guide = 'none'
                              ,range = c(.1,1)) +
  scale_edge_width_continuous(guide = 'none'
                              ,range = c(.05,1)) +
  scale_color_discrete(guide = 'none') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "null")
        #,legend.position = 'none'
        )



# zoom -------------------------------------------------------------------------

ctr <- sfn[sfn$name == '08031002701',]
zoom.lvls <- list(
  m = visaux::cntr2bbx(ctr, 9, st_crs(sfn)),
  l = visaux::cntr2bbx(ctr, 15, st_crs(sfn))
)


ght <- gh %>%
  activate('nodes') %>%
  st_crop(zoom.lvls$l)

ght %>% divseg::get_nodes() %>% mapview::mapview()

ght %>% divseg::get_edges() %>% pull(tie.strength) %>% quantile(seq(0,1,.25))

ght
sfn <- ght %>% divseg::get_nodes()
sfe <- ght %>% divseg::get_edges()


# get some of the extra layers
options(tigris_use_cache = TRUE)
sfn <- st_transform(sfn
                    ,4326)

wtr <- visaux::water.wrapper(x = sfn)

hwypath <- paste0(Sys.getenv('drop_dir'),
                  'shapefiles/nhpn/') %>% list.files(full.names = T,
                                                     pattern = "shp$")


hwys <- st_read(hwypath
                ,wkt_filter=
                  spocc::bbox2wkt(bbox = st_bbox(sfn))
)

hwys <-
  hwys %>%
  filter_at(vars(matches("^SIGNT[1-3]"))
            ,any_vars(. %in% c('I', "U")))

map(list(sfn, wtr, hwys)
    ,st_crs)


# get coords
lonlats <- sfn %>% st_coordinates()

ghplot <- ght %>%
  activate('edges') %>%
  filter(tie.strength > .006) %>%
  mutate(ttstr =
           log(1 + tie.strength)) %>%
  mutate(to_prop.fam.75p =
           .N()$prop.fam.75p[to]) %>%
  ggraph(layout = lonlats  ) +
  geom_edge_density(aes(edge_fill = to_prop.fam.75p)) +

  geom_sf(data = st_crop(wtr
                         ,zoom.lvls$l)
          , color = NA, # feature layers
          fill = '#94bdff') +
  geom_sf(data = st_crop(hwys # feature layers
                         ,zoom.lvls$l)
          , aes(color = SIGNT1)
          , size = 1.5) +

  geom_edge_fan(aes(edge_alpha =
                      ttstr#n
                    ,edge_color =
                      to_prop.fam.75p
                    ,edge_width =
                      ttstr#n
  )) +
  scale_edge_fill_viridis(discrete = T)  +
  scale_edge_color_viridis(discrete = T) +
  scale_edge_alpha_continuous(guide = 'none'
                              ,range = c(.1,1)) +
  scale_edge_width_continuous(guide = 'none'
                              ,range = c(.05,1)) +
  scale_color_discrete(guide = 'none') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "null")
        ,legend.position = 'bottom'
  )
ghplot


ggsave(plot=ghplot,
       "denver-hhincome-flows-rstr.tiff"
       , dpi = 400
       , device = "tiff")



# zooming to >.5 ---------------------------------------------------------------

sfn %>%
  mutate(prop.fam.75p =
           n.fam75p/hh) %>%
  mapview::mapview(zcol = 'prop.fam.75p')

sfn$n.fam75p
# people visiting only tracts >.5% of fams are in top 25% of distrubtion
gh %>%
  activate('edges') %>%
  filter({.N()$n.fam75p[to] /
           .N()$hh[to]} > .5
           ) %>%
  divseg::getno('nodes')



sfn %>%
  mutate(prop.fam.75p =
           n.fam75p / hh)




# ------------------------------------------------------------------------------

sfn <- gh %>%
  divseg::get_nodes()

sfn
# get coords
sfn <- sfn %>%
  cbind(
  st_coordinates(sfn))


options(tigris_use_cache = TRUE)

geox::rx
brm <-
  res %>%
  rename(geoid = geoid19) %>%
  sfg.seg::geo.subset.cbgs('geoid',
                           cbsa_id = '13820') %>%
  geox::attach.geos()

brml <-
  brm %>%
  tibble() %>%
  select(geoid, pop, hh, matches('wh|hsp|bl|asian')) %>%
    pivot_longer(matches("^n\\.")) %>%
    mutate(name = gsub("n.", "", name))


brml <- brml %>%
  left_join(brm['geoid'])

brml <-
  cbind(brml
      ,st_coordinates(
        st_centroid(
          divM::conic.transform(
          st_sf(brml$geometry))))
      )

brml %>%
  ggplot(  ) +
  geom_density_2d_filled(
    aes( x=X
        ,y=Y
        #,color = name
        ,fill = name
        ,group = name
        ,alpha = value
    )
    #,stat='identity'
  ) #+

?after_stat()
ggplot(  ) +
  geom_den
  stat_density2d(aes(x=X ,y=Y, z=value, color=name, alpha=..level..
                     ,fill = name),
               data=brml, size=2, contour=TRUE)


geom_sf(data = wtr
          ,fill = "blue"
          ,color = NA)
  ?geom_bin_2d()
