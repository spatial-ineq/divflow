
# prep data fcns ---------------------------------------------------------------




#' bundle.divseg.cts
#'
#' Bundles data from a variety of sources for a given cz/cbsa. Data will include:
#'
#' `ctd` - prepped in divseg; all control variables and incoming/outgoing aggregated
#' flows by demographic characteristics
#'
#' `divc` - preppedin divM; measures of "dividedness"; whether tracts are across a
#' division from one another
#'
#' & tract centroid geometries retrieved with `tigris` and prepped with `sf`
#'
#' @inheritParams get.btwn.tract.dividedness
#' @inheritDotParams get.btwn.tract.dividedness
#' @inheritDotParams get.divseg.prepped
#'
bundle.divseg.cts <- function(region.ids,
                              add.centroid.geos = F,
                              ddir = Sys.getenv("drop_dir"),
                              ...) {

  # get tract-level sfg/seg data
  ctd <- get.divseg.prepped("ct", ...)

  # add perc.acres.sfz
  ctd <- ctd %>% rename(perc.acres.sfz = acres_perc_sfr)

  # population per sq km
  ctd$pop.dens <- with(ctd, pop / aland)

  # get tract/btwn-tract division measures
  divc <- get.btwn.tract.dividedness(region.ids,
                                     ddir, ...)

  # join all attributes
  out <- left_join(divc,
                     ctd)

  # port in geographies if appropriate
  if(add.centroid.geos) {
    require(sf)

    ctsf <- divM::tracts.from.region(region.ids) %>% select(geoid)
    ctsf <- ctsf %>% divM::conic.transform()
    ctrs <- st_centroid(ctsf)

    out <- left_join(out, ctrs)
  }
  return(out)
}


#' setup.od.as.network
#'
#' Pulls origin-destination flows, estimated from safegraph data and sets up as a
#' graph/network. Most parameters are for intermediate filters and wrapping other
#' steps
#'
#' @param cz,cbsa region to setup CT graph for (cbsa not implemented yet)
#' @param directed Whether to return directed graph (undirected if false)
#' @param flow.floor Because we do daily averages with annual data to get estimated
#'   flows, there are many o-d flows < 1 (for example, about 70% of the 13 million
#'   directed cross-tract flows in Philadelphia cz are less than 1 and 95% are less
#'   than 3.5). Option here to filter based on a floor immediately upon loading.
#' @param pop.floor population floor for each tract
#' @param turn2tracts Safegraph origin-destination data is saved at block groups. If
#'   this is true, flow counts will be aggregated to tracts by dropping the last
#'   digit of the block group identifiers.
#' @param add.centroid.geos Adds tract centroid geometries to node attributes. Will
#'   only work for CZs just to avoid too much added complexity here
#' @param round.flows the ergm models give a warning (and estimation seems to be
#'   slowed down?) when valued edges are doubles instead of integers. Will return
#'   integer flows unless this is false.
#' @param drop.loops Whether or not to drop loops / intratract flows.
#' @param dropbox.dir dropbox directory, where all the data will be pulled from.
#'
#' @return a `tidygraph` graph with nodes representing block groups or tracts and
#'   valued edges representing flows
#'
setup.cts.as.network <- function(cz.id = NULL, cbsa.id = NULL,
                                 directed = F,
                                 flow.floor = 5,
                                 pop.floor = 0,
                                 add.centroid.geos = T,
                                 turn2tracts = F,
                                 round.flows = T,
                                 drop.loops = T,
                                 dropbox.dir = Sys.getenv("drop_dir"),
                                 ...) {

  # get identifiers for region
  r <- divM::get.region.identifiers(cz = cz.id,
                                    cbsa = cbsa.id)

  # get sfg od data from dropbox

  # they are saved by CZ, so load differently based on region type
  # to implement CBSAs can adapt from: https://github.com/spatial-ineq/sfg.seg/blob/main/R/Della-wrapper-fcns.R

  sfgdir <- paste0(dropbox.dir,
                   "sfg-processed/orig_dest_annual/")

  if( !is.null(cbsa.id) )
    czs2load <-
      xwalks::ctx %>%
      select(contains("cz"),
             contains("cbsa")) %>%
      filter(cbsa == cbsa.id) %>%
      distinct() %>%
      pull(cz)
  else
    czs2load <-
      cz.id

  od <- sfg.seg::read.sfg.CZs(czs2load,
                              sfgdir,
                              year = "2019")
  # filters and fixes
  od <- od %>% select(origin = 1, dest = 2, n = 3, 4)

  od <- od[od$n > flow.floor, ]

  od <- od %>%
    mutate(across(1:2,
                  ~stringr::str_pad(as.character(.x),
                                    12, "left", pad = "0")))

  # trim to trips w/in cz (start and end)
  od <- sfg.seg::geo.subset.cbgs(od,
                                 subset.cols = c("origin", "dest"),
                                 cz = cz.id,
                                 cbsa = cbsa.id)
  # turn into tracts
  if(turn2tracts)
    od <- sfg.seg::cbg.flows2tracts(od)

  # transform to graph
  require(igraph)
  require(tidygraph)

  gh <- tidygraph::as_tbl_graph(od,
                                directed = directed)

  # get attrs
  attrs <- bundle.divseg.cts(r,
                             add.centroid.geos = add.centroid.geos,
                             dropbox.dir,
                             ...) %>%
    pivot_wider(names_from = div.type,
                names_prefix = "poly.id_",
                values_from = poly.id)

  # add in "city center" attribute (city center defined as largest Plc in CZ)
  if(r$region.type == "cz")
    attrs <- is.in.city.center(attrs, ...)

  gh <- gh %>%
    activate("nodes") %>%
    left_join(attrs,
              by = c("name" = "geoid"))

  # final cleans
  gh <- gh %>%
    activate("nodes") %>%
    filter(pop >= pop.floor)

  if(round.flows)
    gh <- gh %>%
    activate("edges") %>%
    mutate(n =
             as.integer(round(n)))

  if(drop.loops)
    gh <- gh %>%
    activate("edges") %>%
    filter(from != to)

  return(gh)
}



# setup for running ergm -------------------------------------------------------


#' get.normalized.undirected.connectedness
#'
#' Gets a measure of "tie strength" of total visits between tracts, scaled by total
#' flows into and out of each tract. Takes individual node/tract identifiers and an
#' edgelist containing all flows in between them.
#'
#' ..I'd like to make this fcn faster/more efficient
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

  # drop geo; factors -> ints
  gh <- gh %>%
    activate("nodes") %>%
    select(-geometry) %>%
    mutate(across(where(is.factor),
                  as.integer))

  # apply population floor
  if(is.numeric(population.floor))
    gh <- gh %>%
      activate("nodes") %>%
      filter(pop > population.floor)

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

  pop.center <-
    divM::largest.plc.in.cz %>%
    filter(cz.id %in% cz) %>%
    st_transform(st_crs(sfx))

  sbgp <- st_intersects(
    st_point_on_surface(sfx),
    pop.center)

  sfx <-  sfx %>%
    mutate(in.cc =
             lengths(sbgp) > 0)

  if(filter2cc)
    sfx <- sfx %>%
    filter(in.cc)

  return(sfx)
}


# tidygraph helper fcns --------------------------------------------------------


#' get_nodes
#'
#' Convenience fcn to get nodes as tibble from `tidygraph`
#'
#' @export get_nodes
get_nodes <- function(tbl_grph) {

  tbl_grph %>%
    activate("nodes") %>%
    as_tibble()
}


#' get_edges
#'
#' Convenience fcn to get edges as tibble from `tidygraph`
#'
#' @export get_edges
get_edges <- function(tbl_grph) {

  tbl_grph %>%
    activate("edges") %>%
    as_tibble()
}







# Wrapper -- prep for ERGM -------------------------------------------------


#' Wrapper_prep4ergm
#'
#' Wraps all steps to load and prep data for ERGM. Warning that this function
#' can be a little wonky across systems due to the data coming from so many
#' different places. Parameters `dropbox.dir` should point to directory
#' containing relevant datasets in our dropbox, including `dividedness-measures`
#' and `sfg-processed` with O-D flows. Parameter `ddir` should point to a
#' directory containing the prepped datasets that this package can be used to
#' generate, as with the `glomerating-datasets` scripts.
#'
#' Other arguments can alter how data is prepped.
#'
#' @inheritDotParams setup.cts.as.network
#' @inheritDotParams build.pairwise.dst.matrix
#' @inheritDotParams is.in.city.center
#' @param build.dst.mat Whether or not to build a distance matrix using
#'   `build.pairwise.dst.matrix`. Can be false to skip this step if a model
#'   doesn't want distance or wants to do a function of distance not
#'   accommodated by that function yet (exponential decay; power law decay).
#'   (Daraganova et al 2012 has some possibilities for thinking about different
#'   "spatial interaction fcns" for spatial networks.)
#'
Wrapper_prep4ergm <- function(...,
                               build.dst.mat = T) {

  # browser()

  # pull together data (cz/cbsa arguments passed on)
  gh <-
    suppressMessages(
      setup.cts.as.network(...))

  # get point geometrries for nodes & use to build distance/spw matrix
  ctrs <- gh %>% get_nodes() %>% st_sf() %>% select(geoid = name)

  if(build.dst.mat)
    dst.mat <<- build.pairwise.dst.matrix(ctrs, ...)

  # get edges and tie strength measures
  el <- gh %>% get_edges()

  gh <- gh %>%
    activate("edges") %>%
    mutate(tie.strength =
             map2_dbl(el$from, el$to,
                      ~get.normalized.undirected.connectedness(.x, .y,
                                                               el = el))
           )

 return(gh)
}


# running model ----------------------------------------------------------------



#' tbl_graph2ERGM
#'
#' From a tidygraph form, coerce to statnet Network form and run an ERGM.
#' Consolidates final steps to run ERGM after data is prepped. To be called within
#' `Wrapper_prep.and.run.ERGM` or separately on a similarly prepped graph object.
#'
#' Glossary of ERGM terms here:
#' https://cran.r-project.org/web/packages/ergm/vignettes/ergm-term-crossRef.html
#'
#' @param formula a formula specification for the ergm model. Names of attributes and
#'   graph dataset in the model should match those by earlier calls in this function.
#' @param seed random seed, for reproducability
#' @param ... additional arguments to pass onto `ergm` call
#'
#' @export tbl_graph2ERGM
tbl_graph2ERGM <- function(prepped_graph,
                           formula,
                           seed = 1234,
                           ...) {

  require(statnet)

  # define in global env for formula access (i think another way too but this is
  # ez)
  ergh <<- ergm.clean.n.convert(prepped_graph)

  # set seed and run
  set.seed(seed)

  mo <-
    ergm(
      formula,
      response = "tie.strength",
      reference = ~Poisson,
      ...
    )

  return(mo)
}



# convenience functions --------------------------------------------------------



#' spatialize.nodes
#'
#' I save prepped graphs without geometries for portability (they would sometimes get
#' corrupted moving across systems otherwise). This re-attaches tract geometries
#' based on name colm
#'
#' @param gh tidygraph object
#' @param centroids whether to attach geometries as centroids or polygons (default
#'   true)
#' @param region.ids possible region identifiers as retrieved from
#'   `divM::get.region.identifiers`. Defaults to inferring CZ from graph if left NULL
#' @param geoid.colm string for column name that contains tract geoids in the graph
#'   object.
#'
spatialize.nodes <- function(gh, centroids = T,
                             region.ids = NULL, geoid.colm = "name"
                             , crs = 4326) {

  if(is.null(region.ids)) {
    cz.id <- gh %>% activate("nodes") %>% pull(cz) %>% unique()
    region.ids <- divM::get.region.identifiers(cz = cz.id)
  }
  ctsf <- divM::tracts.from.region(region.ids)
  ctsf <- ctsf %>% divM::conic.transform()
  if(centroids)
    ctsf <- st_centroid(ctsf)
  ctsf <- ctsf %>% st_transform(crs)

  gh <- gh %>%
    activate("nodes") %>%
    left_join(ctsf["geoid"],
              by = setNames("geoid", geoid.colm))
  return(gh)
}


#' get_nodes
#'
#' Convenience fcn to get nodes as tibble from `tidygraph`
#'
#' @export get_nodes
get_nodes <- function(tbl_grph) {

  tbl_grph %>%
    activate("nodes") %>%
    as_tibble()
}


#' get_edges
#'
#' Convenience fcn to get edges as tibble from `tidygraph`
#'
#' @export get_edges
get_edges <- function(tbl_grph) {

  tbl_grph %>%
    activate("edges") %>%
    as_tibble()
}



# in progress -------------------------------------------------------------


# resources and refernces -----------------------------------------------------------


#' http://statnet.org/Workshops/valued.html#32_modeling_ordinal_relational_data_using_ergmrank
#'
#' https://cran.r-project.org/web/packages/ergm/vignettes/ergm-term-crossRef.html
