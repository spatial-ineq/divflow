# Wrapper -- prep for ERGM -------------------------------------------------


#' Wrapper_prep4ergm
#'
#' This function is outdated and I will have to revise to reflect other changes
#'
#' old: Wraps all steps to load and prep data for ERGM. Warning that this function
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
      # setup.cts.as.network(...))

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
