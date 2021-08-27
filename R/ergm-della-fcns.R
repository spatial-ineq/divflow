#' Wrapper fcns to send to della w/ rslurm



# by cz -------------------------------------------------------------------


#' Della.wrapper_cz.ergm.variations
#'
#' A metawrapper. Preps graphs, runs a series of ERGM models, saves
#' visualizations, and saves all three sets of objects. A lot is hardcoded for
#' now. In particular, NHPN hwy data should be sent along with the fcn using the
#' `rslurm::slurm_apply` `add_objects` argument (called `nhpn`), as well as a
#' list of ergm specifications titled `model.variations.`
#'
#' Will maybe change so it requires a proximity matrix or list of prx.matrices
#' at some point..
#'
#' @inheritParams Wrapper_prep4ergm
#'
Della.wrapper_cz.ergm.variations <- function(
  cz,
  directed = F,
  flow.floor = 1,
  drop.loops = T,
  pop.floor = 500,
  dropbox.dir = "/scratch/gpfs/km31/"
  ,ddir =
    "/home/km31/all/divseg/.input-datasets/"
  ,seed = 1234
) {

  require(tidyverse)
  require(sf)
  require(statnet)
  require(igraph)
  require(tidygraph)
  require(divseg)
  requireNamespace("divM")

  # refresh hwy crs
  st_crs(nhpn) <- "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45"

  # get region identifiers
  r <- divM::get.region.identifiers(cz = cz)
  # use to build a base filename for saving objects
  base.fn <- paste(r, collapse = "-")

  print(base.fn)

  gh <-
    Wrapper_prep4ergm(
      cz =
        cz
      ,directed = directed
      ,flow.floor = flow.floor
      ,drop.loops = drop.loops
      ,pop.floor = pop.floor
      ,dropbox.dir = dropbox.dir
      ,ddir =
        ddir
      ,build.dst.mat = FALSE
    )

  # get point geometries for nodes & use to build distance/spw matrix i do it
  # outside of wrapper to provide options
  ctrs <- gh %>% get_nodes() %>% st_sf() %>% select(geoid = name)
  dst.mat <<- build.pairwise.dst.matrix(ctrs)
  decay.mat <<- build.pairwise.dst.matrix(ctrs, decay = T, h0 = 32.2)

  # add ln population dens
  gh <- gh %>%
    activate("nodes") %>%
    mutate(ln.pop.dens = log(pop.dens))

  # save graph
  gh.base.dir <- "/scratch/gpfs/km31/ergms/cz-graphs/"
  gh.sv.path <- paste0(gh.base.dir,
                       base.fn, ".rds")
  saveRDS(gh, file = gh.sv.path)

  # switch to statnet libraries
  ergh <<- ergm.clean.n.convert(gh)

  # set seed
  set.seed(seed)

  # run models
  model.base.dir <<- "/scratch/gpfs/km31/ergms/models/"

  for(fo.str in names(model.variations)) {

    print(paste(base.fn, fo.str))

    mo.save.path <- paste0(model.base.dir,
                           fo.str, "-", base.fn, ".rds")
    mo <-
      # Wraps in try/catch
      tryCatch(
        ergm(
          model.variations[[fo.str]],
          response = "tie.strength",
          reference = ~Poisson
        )
        ,error = function(cond) {
          message(cond)
          # Choose a return value in case of error
          return(cond)
        })
    # saves regardless--- so error message will be saved as rds if arises
    saveRDS(mo, file = mo.save.path)
  }

  # also save a vis
  print("prepping vis..")
  vis.base.dir <- "/scratch/gpfs/km31/ergms/viz/"
  tryCatch(
    Wrapper_make.graph.map(gh,
                           cz.id = cz,
                           trim2plc = T,
                           save.dir = vis.base.dir)
    ,error = function(cond) {
      message(cond)
      # Choose a return value in case of error
      return(cond)
    })

  print(paste("finished", base.fn))
  return(1)
}

#' Della.wrapper_sif.variations
#'
#' Not a graceful fcn; include `model.variations` & `prx.mats` in
#' `include_objects` rslurm call.
#'
Della.wrapper_sif.variations <- function(cz, seed = 1234) {

  require(tidyverse)
  require(sf)
  require(statnet)
  require(igraph)
  require(tidygraph)
  require(divseg)
  requireNamespace("divM")

  gh <-
    divseg::Wrapper_prep4ergm(
      cz = cz,
      directed = F,
      flow.floor = 1,
      pop.floor = 500,
      dropbox.dir = "/scratch/gpfs/km31/"
      ,ddir =
        "/home/km31/all/divseg/.input-datasets/"
    )

  # add ln population dens
  gh <- gh %>%
    activate("nodes") %>%
    mutate(ln.pop.dens = log(pop.dens))

  ergh <<- ergm.clean.n.convert(gh)

  # get region identifiers
  r <- divM::get.region.identifiers(cz = cz)
  # use to build a base filename for saving objects
  base.fn <- paste(r, collapse = "-")

  ergm.dir <- "/scratch/gpfs/km31/ergms/"
  sif.svpath <-
    paste0(ergm.dir,
           "sif-models/")

  set.seed(seed)

  for(fo.str in names(model.variations)) {

    print(paste(base.fn, fo.str))

    mo.save.path <- paste0(sif.svpath,
                           fo.str, "-", base.fn, ".rds")
    mo <-
      # Wraps in try/catch
      tryCatch(
        ergm(
          model.variations[[fo.str]],
          response = "tie.strength",
          reference = ~Poisson
        )
        ,error = function(cond) {
          message(cond)
          # Choose a return value in case of error
          return(cond)
        })
    # saves regardless--- so error message will be saved as rds if arises
    saveRDS(mo, file = mo.save.path)
  }

  print(paste("finished sif variations", base.fn))
  return(1)
}


# reading output fcns -----------------------------------------------------


#' fn2tags
#'
#' Gets region and/or formula tags from a filename saved accoring to conventions
#' used by Della wrapper fcns here.
#'
#' @inheritParams model2tidy
#'
fn2tags <- function(fn, tag = c("region", "formula", "all")) {

  require(tidyverse)
  tags <-  fn %>%
    gsub(".rds$", "", .) %>%
    strsplit("-cz-") %>% unlist()

  if(tag[1] == "region")
    return(tags[2])
  if(tag[1] == "formula")
    return(tags[1])

  tags <- paste(tags, collapse = "-")
  return(tags)
}


#' model2tidy
#'
#' Loads a model .rds saved with the della script, extracts some diagnostics and
#' coefficient estimates, and saves to a running .csv with all model outputs.
#' Designed to map across filenames
#'
#' @param fn,mdir a filename of saved model .rds, and it's directory.
#' @param svpath Path to a .csv to which to write mass model summary
#'
#'
model2tidy <- function(fn,
                       mdir = "/scratch/gpfs/km31/ergms/models/"
                       ,svpath =
                         "/scratch/gpfs/km31/ergms/model summaries/cz-ergms.csv"
) {

  #browser()

  require(tidyverse)
  require(statnet)

  # read
  mo <- readRDS(paste0(mdir, fn))

  # get labels for region & formula
  fos <- fn2tags(fn, "formula")
  regions <- fn2tags(fn, "region")

  print(paste("parsing coef estimates from", regions, fos))

  # initialize tibble w/ labels
  mss <- tibble(
    region = regions,
    fo = fos
  )

  # add model
  mss$model <- list(mo)

  # add coefficients and unnest
  mss$coefs <-
    list(broom::tidy(mo))

  mss <- mss %>% unnest(cols = coefs)

  # add AIC
  mss$aic <- summary(mo)$aic

  # drop model object for summary table
  mss <- mss %>%
    select(-model)

  # use earlier fcn for writing a running table
  sfg.seg::write.running.table(
    mss, svpath
  )

  return(1)
}

