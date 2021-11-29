



# building distance and proximity matrices --------------------------------


#' build.pairwise.dst.matrix
#'
#' Builds a pairwise distance matrix from an sf object. Option to get bisq distance
#' wights. Mostly just wraps `sf::st_distance`. See
#' http://mjh4.blogspot.com/2012/09/ergm-edgecov-and-dyadcov-specifications.html or
#' https://cran.r-project.org/web/packages/ergm/vignettes/ergm-term-crossRef.html for
#' `edgecov`/`dyadcov` documentation.
#'
#' @param sfx sf object to generate distances from
#' @param decay Whether to get decayed distances with bisq weights
#' @param h0 if `decay` is `TRUE`, the distance in KM at which spatial weights are 0.
#'
#' @return distance or spatial weights matrix. If distance, matrix is of numerics
#'   representing KMs
#'
#' @export build.pairwise.dst.matrix
build.pairwise.dst.matrix <- function(sfx, decay = F, h0 = 32.2, ...) {

  dst.mat <- st_distance(sfx)
  dst.dim <- dim(dst.mat)

  # drop explicit units tag and convert to KM
  dst.mat <- as.numeric(dst.mat) / 1000

  if(decay)
    dst.mat <-
    spgwr::gwr.bisquare(dst.mat^2, h0)

  dst.mat <-
    matrix(dst.mat,
           ncol = dst.dim[1])

  return(dst.mat)

}

# ties across distance, relative to all centroid distances ----------------

#' relative.tie_dst.frequency
#'
#' Gets distribution of ties relative to distribution of centroids. Because centroids
#' are not an even continuous distribution, this might show a more meaningful spatial
#' distribution of ties.
#'
#' Implements technique in Daraganova 2012, referred to as  "binned relative
#' frequency of ties against distance interval"
#'
#' @param sfn,sfe,gh Spatial nodes, spatial edges or graph object. Save load time by
#'   providing spatial nodes / edges; otherwise just supply graph gh.
#' @param n.bins number of distance bins
#'
#' @export relative.tie_dst.frequency
relative.tie_dst.frequency <- function(sfn = NULL, sfe = NULL,
                                       gh = NULL, n.bins = 50, dst.mat = NULL) {

  if(is.null(sfn) | is.null(sfe)) {
    # spatial nodes
    sfn <- gh %>% activate('nodes') %>%
      as_tibble() %>% st_sf()

    # spatial edges
    sfe <-  gh %>%
      activate("edges") %>%
      sfnetworks::as_sfnetwork(., directed = T,
                               edges_as_lines = T) %>%
      get_edges() %>%
      st_sf()

    # edge distances
    sfe$dst <- st_length(sfe$geometry)
    sfe$dst <- as.numeric(sfe$dst) / 1000 # to km
  }

  # browser()

  # build raw distance matrix
  if(is.null(dst.mat))
    dst.mat <- build.pairwise.dst.matrix(sfn)

  # all xwise dsts
  xwise <-
    tibble(dst = dst.mat[lower.tri(dst.mat)] )

  # get dst breaks
  padded.max <- ceiling(max(xwise$dst)+1e-5)
  dst.brks <-
    seq(0, padded.max,
        by = padded.max /
          n.bins)

  # counts of pairwise dsts for each distance bin (the denominator in finding
  # average flows per distance, based on underlying nodewise geography.)
  xwise$dst.bin <-
    cut(xwise$dst,
        breaks = dst.brks,
        include.lowest = T,
        dig.lab = 5)
  binned.xwise <- xwise %>% count(dst.bin, name = "xwise.count")

  # put actual flows in those same bins
  sfe <-
    tibble(sfe) %>%
    mutate(dst.bin =
             cut(dst,
                 breaks = dst.brks,
                 include.lowest = T,
                 dig.lab = 5))
  dst.freq <-
    sfe %>%
    group_by(dst.bin) %>%
    summarise(flow.count = sum(n))

  # merge in xwise counts
  dst.freq <- dst.freq %>%
    left_join(binned.xwise)

  # average flows for each distance bin
  dst.freq$avg.flow <-
    with(dst.freq, flow.count/xwise.count)

  # add dst midpoints (can be nicer for plotting etc.)
  dst.freq$dst.lvl <-
    appHelpers::get_mean_from_interval(dst.freq$dst.bin)

  # rearrange cols
  dst.freq <- dst.freq %>%
    select(dst.bin, dst.lvl, xwise.count, flow.count, avg.flow)

  return(dst.freq)
}



# maximum likelihood fitting ----------------------------------------------

#' for the parameters for the spatial interaction fucntion

#' power.law.loss.fcn
#'
#' Loss function that wraps the `power.law` fcn. Uses the `dst.freq` object, as
#' returned from `relative.tie_dst.frequency`. This function is basically copied
#' from https://rpubs.com/YaRrr/MLTutorial
#'
#' @param par vector of parameters passed onto `power.law`
#'
power.law.loss.fcn <- function(par) {

  pb <- par[1]
  gamma <- par[2]
  alpha <- par[3]
  err.sigma <- par[4]

  if(err.sigma < 0)
    deviance <- 10000000
  # If the error standard deviation is valid (i.e.; > 0), then calculate the deviance...
  if(err.sigma > 0) {

    # Calculate the likelihood of each data point.
    # Here is where you specify the model and how you calculate likelihoods.
    likelihoods <-
      dnorm(dst.freq$avg.flow,
            mean =
              #attenuated.power.law
              power.law(dst.freq$dst.lvl,
                        pb, gamma, alpha),
            sd = err.sigma)

    # Because we don't like likelihoods of 0 (which should rarely happen), convert 0s
    # to a very small number
    likelihoods[likelihoods == 0] <- .001

    # Now let's convert the vector of likelihoods to a summary deviance score (-2 times sub log-lik)

    # Calculate log-likelihoods of each data point
    log.likelihoods <- log(likelihoods)

    # Calculate the deviance (-2 times sum of log-likelihoods)
    deviance <- -2 * sum(log.likelihoods)

  }
  return(deviance)
}



# decay and SIF functions ------------------------------------------------------



#' power law
#'
#' @param d distance
#'
power.law <- function(d, pb = 1, gamma = 1, alpha = 1
                      ,...) {
  pb * (1 + alpha*d) ^ -gamma
}


#' attenuated power law
#'
#' "Unlike the simple power law model....the attenuated variant exhibits
#' negative local curvature. This model thus tends to produce both a high proportion
#' of long-range edges and a high degree of local cohesion (the scale of these
#' effects being dependent upon α). Where (αd)γ » 1, both variants exhibit similar
#' behavior."
#'
#' @inheritParams power.law
#'
attenuated.power.law <- function(d, pb = 1, gamma = 1, alpha = 1
                                 ,...) {
  pb * (1 + (alpha*d) ^ -gamma)
}


#' inverse exponential decay
#'
#' Fd(x = pb exp(–αx), α >=0 scaling param
#'
#' @inheritParams power.law
#'
exponential.decay <- function(d, pb = 1, alpha = 0.5
                              ,...) {
  {pb * exp(-alpha * d)}
}


#' bisq.dist2weights
#'
#'
#' Uses a Bisquare transformation to turn distances to spatial weights. Transfroms by
#' weight = (1-(ndist/H)^2)^2 for distances less than or equal to H, 0 otherwise.
#' Mimics Bi-square option in `lctools::moransI`.
#'
#' @inheritParams power.law
#' @param cutoff `H` in the formula, where weight is 0.
#'
#'
bisq.dist2weights <- function(d
                              , cutoff = dst.ceiling
                              , ... ) {

  if_else( d >= cutoff
           ,0
           ,(1-as.numeric(d/cutoff)^2)^2)

}

#' neg.exp
#'
#' e^-distance. Used by Massey and Denton for the Clustering dimension of
#' segregation. Notably they include tract where i==j and use an estimated
#' self-distance
#'
#' @inheritParams power.law
#'
neg.exp <- function(d, ...) {
  exp(-d)
}
