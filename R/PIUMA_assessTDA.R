#' @title Assessment of Scale-Free model fitting
#'
#' @description This function assesses the fitting to a scale-free net model.
#'
#' @param x A TDAobj object, processed by the  \code{\link{jaccardMatrix}}
#'
#' @param showPlot Whether the  plot has to be generated. Default: FALSE
#'
#' @return A list containing:
#' \itemize{
#'    \item connectivity of the resulting graph
#'    \item the estimated gamma value
#'    \item the correlation between degree \(k\) and its distribution \(p(k)\).
#'    \item The p-value of the correlation between the k and the degree
#'    distribution p(k).
#'    \item The correlation between the logarithm (base 10) of k and the
#'    logarithm (base 10) of the degree distribution p(k).
#'    \item The p-value of the correlation between the logarithm (base 10) of k
#'    and the logarithm (base 10) of the degree distribution p(k).
#'    \item A composite score reflecting how strongly power‐law behavior
#'    coexists with graph cohesion, computed as the absolute product between
#'    cor(P(k)*k) and connectivity}
#'
#' @details The scale-free networks show a high negative correlation beztween k
#' and p(k).
#'
#' @author  Mattia Chiesa, Laura Ballarini, Luca Piacentini, Carlo Leonardi
#'
#' @examples
#' ## use example data:
#' data(tda_test_data)
#' netModel <- checkScaleFreeModel(tda_test_data)
#' print(netModel)
#'
#' @seealso
#' \code{\link{makeTDAobj}},
#' \code{\link{dfToDistance}},
#' \code{\link{dfToProjection}},
#' \code{\link{mapperCore}},
#' \code{\link{jaccardMatrix}}
#'
#' @export
checkScaleFreeModel <- function(x, showPlot = FALSE) {
  if (!is(x, "TDAobj")) {
    stop("'x' argument must be a TDAobj object")
  }
  jaccIndexes <- getJacc(x)
  if (missing(jaccIndexes)) {
    stop("'jaccIndexes' argument must be provided")
  }
  if (!is.matrix(jaccIndexes)) {
    stop("'jaccIndexes' argument must be a matrix")
  }
  if (!is.logical(showPlot)) {
    stop("'showPlot' argument must be TRUE or FALSE")
  }
  if (!all((jaccIndexes >= 0 & jaccIndexes <= 1) | is.na(jaccIndexes))) {
    stop("'jaccIndexes' must be in [0; 1] or NA")
  }
  if (length(showPlot) > 1) {
    stop("length(showPlot) must be equal to 1")
  }
  if (any(is.infinite(jaccIndexes))) {
    stop("Inf values are not allowed in the 'jaccIndexes' matrix")
  }

  infoModelRes <- list()
  adjDataRes <- jaccIndexes
  adjDataRes[is.na(adjDataRes)] <- 0
  adjDataRes[adjDataRes > 0] <- 1
  graphFromAdjMatr <- igraph::graph_from_adjacency_matrix(adjDataRes,
                                                          mode = "undirected")
  clusters <- igraph::clusters(graphFromAdjMatr)
  largest_comp_size <- max(clusters$csize)
  total_nodes <- igraph::vcount(graphFromAdjMatr)
  infoModelRes[["Connectivity"]] <- largest_comp_size / total_nodes
  dataPl <- igraph::degree(graphFromAdjMatr, mode = "in")
  dataPlDist <- data.frame(k = 0:max(dataPl),
                           pk = igraph::degree_distribution(graphFromAdjMatr))
  dataPlDist <- dataPlDist[dataPlDist$pk > 0 & dataPlDist$k > 0, ]

  if (max(dataPl) >= 5) {
    powLawFit <- igraph::fit_power_law(dataPl + 1, round(max(dataPl) / 5))
    infoModelRes[["gamma"]] <- powLawFit$alpha
    logk <- log10(dataPlDist$k)
    logpk <- log10(dataPlDist$pk)
    k <- (dataPlDist$k)
    pk <- (dataPlDist$pk)

    if (length(dataPlDist$k) >= 5) {
      rCorkpk <- Hmisc::rcorr(k, pk)
      rCorlogklogpk <- Hmisc::rcorr(logk, logpk)
      infoModelRes[["corkpk"]] <- rCorkpk[["r"]][2]
      infoModelRes[["pValkpk"]] <- rCorkpk[["P"]][2]
      infoModelRes[["corlogklogpk"]] <- rCorlogklogpk[["r"]][2]
      infoModelRes[["pVallogklogpk"]] <- rCorlogklogpk[["P"]][2]
    } else {
      infoModelRes[["corkpk"]] <- NA
      infoModelRes[["pValkpk"]] <- NA
      infoModelRes[["corlogklogpk"]] <- NA
      infoModelRes[["pVallogklogpk"]] <- NA
      warning("Insufficient data points (k < 5). Scale-free metrics set to NA.")
    }

    infoModelRes[["ProductScore"]] <- abs(infoModelRes$corlogklogpk) * infoModelRes$Connectivity

    if (showPlot) {
      plot_ScaleFreeLaw(dataPlDist, rCorkpk, rCorlogklogpk)
    }
  } else {
    infoModelRes[["gamma"]] <- NA
    infoModelRes[["corkpk"]] <- NA
    infoModelRes[["pValkpk"]] <- NA
    infoModelRes[["corlogklogpk"]] <- NA
    infoModelRes[["pVallogklogpk"]] <- NA
    infoModelRes[["ProductScore"]] <- NA
    warning("Network too sparse (max degree < 5). All metrics set to NA.")
  }

  return(infoModelRes)
}



#' @title Compute the Network Entropy
#'
#' @description This function computes the average of the entropies for
#' each node of a network.
#'
#' @param outcome_vect A vector containing the average outcome values for each
#' node
#' of a network.
#'
#' @return The network entropy using each node of a network.
#'
#' @details The average of the entropies is related to the amount of information
#' stored in the network.
#'
#' @author Mattia Chiesa, Laura Ballarini, Luca Piacentini
#'
#' @examples
#' # use example data:
#' set.seed(1)
#' entropy <- checkNetEntropy(round(runif(10), 0))
#'
#' @seealso
#' \code{\link{makeTDAobj}},
#' \code{\link{dfToDistance}},
#' \code{\link{dfToProjection}},
#' \code{\link{mapperCore}},
#' \code{\link{jaccardMatrix}},
#' \code{\link{tdaDfEnrichment}}
#'
#' @export
checkNetEntropy <- function(outcome_vect) {
  # checks----------------------------------------------------------------------
  # check missing arguments
  if (missing(outcome_vect)) {
    stop("'outcome_vect' argument must be provided")
  }

  # check the type of argument
  if (!is.numeric(outcome_vect)) {
    stop("'outcome_vect' argument must be numeric")
  }

  # check the presence of NA or Inf
  if (any(is.na(outcome_vect))) {
    stop("NA values are not allowed in the 'vectFreq'")
  }

  if (any(is.infinite((outcome_vect)))) {
    stop("Inf values are not allowed in the 'vectFreq'")
  }

  # body -----
  outcome_vect_r <- round(outcome_vect, 0)
  uniq_classes <- unique(sort(outcome_vect_r))

  p_logp <- c()
  for (i in uniq_classes) {
    prob_class <- length(which(outcome_vect_r %in% i)) / length(outcome_vect_r)
    p_logp[i] <- prob_class * log(prob_class, base = 2)
  }
  outcome_entropy <- round(-sum(p_logp), 3)
  return(outcome_entropy)
}


#' Predict Mapper Graph Geometry (lightweight)
#'
#' @description
#' Infer a geometry label for \code{x@graph$igraph} using fast heuristics.
#' Writes only \code{x@graph$predicted$class} (one of
#' \code{c("SF","RGG","WS","ER","SBM","CM")}).
#'
#' #' @details
#' Heuristics (hierarchical decision):
#' \itemize{
#'   \item \strong{SF (relaxed)}: rely on \code{checkScaleFreeModel(x)}. Declare
#'         scale-free if at least one of the following holds:
#'         \code{|corlogklogpk| >= 0.55}, \code{|corkpk| >= 0.70},
#'         \code{1.6 <= gamma <= 3.6}, or \code{Connectivity >= 0.40};
#'         alternatively accept SF if the product score
#'         \code{|corlogklogpk| * Connectivity >= 0.2}.
#'   \item \strong{WS}: small-world index \code{sigma > 1.2} with
#'         \code{C/C_ER >= 3} and \code{L/L_ER <= 1.2}.
#'   \item \strong{RGG}: very high clustering vs ER (\code{C/C_ER >= 5}),
#'         longer paths (\code{L/L_ER >= 1.3}), and positive degree
#'         assortativity (\code{r >= 0.10}).
#'   \item \strong{ER}: Poisson-like degree dispersion
#'         \code{VMR = var(k)/mean(k) ~ 1} (within 30\%),
#'         \code{|C - p| <= 0.05}, \code{|r| <= 0.05}, and
#'         \code{0.8 <= sigma <= 1.2}.
#'   \item \strong{SBM}: strong modular structure, \code{Q >= 0.40} with
#'         \code{>= 3} communities.
#'   \item \strong{CM}: heterogeneous degrees (\code{var(k)/mean(k) > 2}) with
#'         clustering close to ER (\code{|C - C_ER| <= 0.05}); otherwise use a
#'         sigma-based fallback (WS if \code{sigma > 1.2}, else ER).
#' }
#' The function sets only \code{x@graph$predicted}. It is intentionally
#' lightweight for fast computation.
#'
#' @param x A \code{TDAobj} with \code{x@graph$igraph} set.
#' @param verbose Logical; print the chosen label. Default \code{FALSE}.
#'
#' @return The input \code{TDAobj} with \code{x@graph$predicted} set.
#'
#' @author Carlo Leonardi, Mattia Chiesa
#'
#' @examples
#' data(tda_test_data)
#' #tda_test_data <- predict_mapper_class(tda_test_data)
#'
#' @export
predict_mapper_class <- function(x, verbose = FALSE) {
  if (!methods::is(x, "TDAobj")) {
    stop("'x' must be a TDAobj")
  }
  g <- x@graph$igraph
  if (is.null(g) || !igraph::is_igraph(g)) {
    stop("`x@graph$igraph` is not set;
         build it via jaccardMatrix() and setGraph().")
  }
  comp <- igraph::components(g)
  gc_id <- which.max(comp$csize)
  g_gc  <- igraph::induced_subgraph(g, vids = which(comp$membership == gc_id))
  n <- igraph::vcount(g_gc)
  m <- igraph::ecount(g_gc)
  if (n < 3L || m == 0L) {
    x@graph$predicted <- "ER"
    if (isTRUE(verbose)) message("predict_mapper_class: ER (degenerate)")
    return(invisible(x))
  }
  k     <- igraph::degree(g_gc)
  kbar  <- mean(k)
  p     <- if (n > 1L) (2 * m) / (n * (n - 1L)) else NA_real_
  C_obs <- igraph::transitivity(g_gc, type = "global")
  L_obs <- igraph::average.path.length(g_gc, directed = FALSE)
  C_er  <- p
  L_er  <- if (is.finite(kbar) && kbar > 1) log(n) / log(kbar) else NA_real_
  sigma <- if (is.finite(C_obs) && is.finite(C_er) && is.finite(L_obs) &&
               is.finite(L_er)  && C_er > 0 && L_er > 0) {
    (C_obs / C_er) / (L_obs / L_er)
  } else NA_real_
  r_deg <- igraph::assortativity_degree(g_gc, directed = FALSE)
  vmr   <- if (is.finite(kbar) && kbar > 0) stats::var(k) / kbar else NA_real_

  com_fg <- try(igraph::cluster_fast_greedy(g_gc), silent = TRUE)
  Q      <- if (!inherits(com_fg, "try-error")) igraph::modularity(com_fg) else
    NA_real_
  ncom   <- if (!inherits(com_fg, "try-error")) {
    length(unique(igraph::membership(com_fg)))
  } else NA_integer_

  th_sf_gamma    <- c(1.6, 3.6)
  th_sf_cor_ll   <- 0.55
  th_sf_cor_lin  <- 0.70
  th_sf_conn_min <- 0.40
  th_sf_prod_min <- 0.2
  th_sigma_ws      <- 1.2
  th_Cratio_ws_min <- 3.0
  th_Lratio_ws_max <- 1.2
  th_Cratio_rgg    <- 5.0
  th_Lratio_rgg    <- 1.3
  th_rdeg_rgg_min  <- 0.10
  th_vmr_er_tol    <- 0.30
  th_Cdiff_er      <- 0.05
  th_sigma_er_low  <- 0.8
  th_sigma_er_high <- 1.2
  th_Q_sbm         <- 0.40
  th_ncom_sbm      <- 3L
  th_vmr_cm        <- 2.0
  th_Cdiff_cm      <- 0.05

  sf_ok <- FALSE
  sf_res <- try(checkScaleFreeModel(x), silent = TRUE)
  if (!inherits(sf_res, "try-error")) {
    cor_ll <- sf_res[["corlogklogpk"]]
    cor_ln <- sf_res[["corkpk"]]
    gamma  <- sf_res[["gamma"]]
    conn   <- sf_res[["Connectivity"]]
    prod   <- abs(cor_ll) * conn

    votes <- sum(c(
      is.finite(cor_ll) && abs(cor_ll) >= th_sf_cor_ll,
      is.finite(cor_ln) && abs(cor_ln) >= th_sf_cor_lin,
      is.finite(gamma)  && gamma >= th_sf_gamma[1L] && gamma <= th_sf_gamma[2L]
    ))
    sf_ok <- is.finite(conn) && conn >= th_sf_conn_min &&
      (votes >= 1L || (is.finite(prod) && prod >= th_sf_prod_min))
  }

  label <- NA_character_

  if (isTRUE(sf_ok)) {
    label <- "SF"
  } else if (is.finite(sigma) &&
             sigma > th_sigma_ws &&
             is.finite(C_obs) && is.finite(C_er) && (C_obs / C_er) >=
             th_Cratio_ws_min &&
             is.finite(L_obs) && is.finite(L_er) && (L_obs / L_er) <=
             th_Lratio_ws_max) {
    label <- "WS"
  } else if (is.finite(C_obs) && is.finite(C_er) && (C_obs / C_er) >=
             th_Cratio_rgg &&
             is.finite(L_obs) && is.finite(L_er) && (L_obs / L_er) >=
             th_Lratio_rgg &&
             is.finite(r_deg) && r_deg >= th_rdeg_rgg_min) {
    label <- "RGG"
  } else if (is.finite(vmr) && abs(vmr - 1) <= th_vmr_er_tol &&
             is.finite(C_obs) && is.finite(p)   && abs(C_obs - p) <=
             th_Cdiff_er &&
             is.finite(sigma) && sigma >= th_sigma_er_low && sigma <=
             th_sigma_er_high &&
             is.finite(r_deg) && abs(r_deg) <= 0.05) {
    label <- "ER"
  } else if (is.finite(Q) && Q >= th_Q_sbm &&
             is.finite(ncom) && ncom >= th_ncom_sbm) {
    label <- "SBM"
  } else if (is.finite(vmr) && vmr > th_vmr_cm &&
             is.finite(C_obs) && is.finite(C_er) && abs(C_obs - C_er) <=
             th_Cdiff_cm) {
    label <- "CM"
  } else {
    label <- if (is.finite(sigma) && sigma > th_sigma_ws) "WS" else "ER"
  }
  x@graph$predicted <- label
  if (isTRUE(verbose)) message("predict_mapper_class: ", label)
  invisible(x)
}
