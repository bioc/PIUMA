#' @title Automatic Clustering of a Mapper Graph by Predicted Geometry
#' (with kNN tie-break)
#' @description
#' Cluster Mapper nodes (x@graph$igraph) with a chosen community algorithm or
#' an automatic selection based on predicted graph geometry
#' (x@graph$predicted).
#' Assign observations either by k-NN tie-breaking (default) or by pure
#' topological label concatenation.
#'
#' @param x A \code{TDAobj} with \code{x@graph$igraph} and
#'                \code{x@graph$predicted$class} set.
#' @param method  One of \code{"automatic","fast_greedy","walktrap",
#'                "edge_betweenness","optimal","label_propagation"}.
#'                Default \code{"automatic"}.
#' @param k Integer >=1 or \code{FALSE}. Default \code{5L}. If numeric, use
#'                k-NN mean distance to break ties; if \code{FALSE}, concatenate
#'                topological labels (e.g. "2_8", "1_2_3").
#'
#' @return The input \code{TDAobj} \emph{invisibly}, with \code{x@clustering}
#' updated:
#' \describe{
#'   \item{\code{nodes_cluster}}{Data frame with columns \code{node},
#'   \code{obs}, \code{cluster}.}
#'   \item{\code{obs_cluster}}{Data frame with columns \code{obs},
#'   \code{cluster}.}
#' }
#'
#' @details
#' In \code{method = "automatic"}, the algorithm is chosen from the
#' predicted geometry:
#' \describe{
#'   \item{\code{SF} / \code{CM}}{Use \emph{fast greedy} modularity
#'   optimization.}
#'   \item{\code{WS}}{Use \emph{Walktrap} (short random walks).}
#'   \item{\code{RGG}}{Use \emph{edge betweenness} (bridge detection).}
#'   \item{\code{SBM}}{Prefer \emph{optimal} (exact modularity) for small
#'   graphs; falls back for larger ones.}
#'   \item{\code{ER}}{Use \emph{label propagation} (fast, parameter-free).}
#' }
#' Isolated nodes (degree = 0) become singletons with unique labels.
#'
#' @author Carlo Leonardi, Mattia Chiesa
#'
#' @examples
#' data(vascEC_norm)
#' data(vascEC_meta)
#' #df_TDA <- cbind(vascEC_meta, vascEC_norm)
#' #df_TDA <- makeTDAobj(df_TDA,outcomes = c("stage","zone"))
#' #df_TDA <- dfToDistance(df_TDA,'euclidean')
#' #df_TDA <- dfToProjection(df_TDA, "UMAP", nComp = 2)
#' #df_TDA <- mapperCore(df_TDA,
#'  #         nBins = 20, overlap = 0.3,
#'  #         mClustNode = 2, clustMeth = "kmeans")
#' #df_TDA <- jaccardMatrix(df_TDA)
#' #df_TDA <- setGraph(df_TDA)
#' #df_TDA <- predict_mapper_class(df_TDA)
#' #df_TDA <- autoClusterMapper(df_TDA,method = 'walktrap')
#'
#' @seealso
#' \code{\link{mapperCore}}
#'
#' @export
autoClusterMapper <- function(x,
                              method = c("automatic", "fast_greedy", "walktrap",
                                         "edge_betweenness", "optimal",
                                         "label_propagation"),
                              k = 5L) {
  method <- match.arg(method)
  if (!inherits(x, "TDAobj")) stop("`x` must be a TDAobj")

  g <- x@graph$igraph
  if (is.null(g) || !igraph::is_igraph(g)) {
    stop("`x@graph$igraph` is not set; run jaccardMatrix(x) and
         setGraph(x, ...) first.")
  }
  if (is.null(x@graph$predicted) ||
      !is.character(x@graph$predicted) ||
      length(x@graph$predicted) != 1) {
    stop("`x@graph$predicted$class` is not set;
         run predict_mapper_class(x) first.")
  }
  geom <- x@graph$predicted
  if (!geom %in% c("ER","SF","WS","RGG","SBM","CM")) {
    stop("`x@graph$predicted` must be one of: ER, SF, WS, RGG, SBM, CM")
  }

  # Ensure names exist and are unique
  vnames <- igraph::V(g)$name
  if (is.null(vnames)) {
    vnames <- paste0("v", seq_len(igraph::vcount(g)))
    igraph::V(g)$name <- vnames
  }

  # Robust: undirected + simplify on the subgraph used for clustering
  if (igraph::is_directed(g)) g <- igraph::as.undirected(g, mode = "collapse")
  degs     <- igraph::degree(g)
  isolates <- which(degs == 0)
  noniso   <- which(degs > 0)

  # Pick community function
  comFun <- if (method == "automatic") {
    switch(geom,
           SF  = igraph::cluster_fast_greedy,
           CM  = igraph::cluster_fast_greedy,
           WS  = igraph::cluster_walktrap,
           RGG = igraph::cluster_edge_betweenness,
           SBM = igraph::cluster_optimal,
           ER  = igraph::cluster_label_prop
    )
  } else {
    switch(method,
           fast_greedy       = igraph::cluster_fast_greedy,
           walktrap          = igraph::cluster_walktrap,
           edge_betweenness  = igraph::cluster_edge_betweenness,
           optimal           = igraph::cluster_optimal,
           label_propagation = igraph::cluster_label_prop,
           stop("Unknown method: ", method)
    )
  }

  # Cluster only non-isolates
  g_sub <- igraph::induced_subgraph(g, vids = noniso)
  g_sub <- igraph::simplify(g_sub, remove.multiple = TRUE, remove.loops = TRUE)

  # Guard: cluster_optimal is exact & expensive
  if (identical(comFun, igraph::cluster_optimal) &&
      igraph::vcount(g_sub) > 200) {
    warning("SBM to cluster_optimal skipped for >200 vertices;
            using Walktrap fallback.")
    comFun <- igraph::cluster_walktrap
  }

  com_sub  <- if (igraph::vcount(g_sub) > 0) comFun(g_sub) else NULL
  memb_sub <- if (!is.null(com_sub)) igraph::membership(com_sub) else
    integer(0L)

  # Map membership by vertex names (not positions)
  full_memb <- setNames(rep(NA_character_, length(vnames)), vnames)
  if (length(memb_sub)) {
    sub_names <- igraph::V(g_sub)$name
    if (is.null(names(memb_sub))) names(memb_sub) <- sub_names
    full_memb[names(memb_sub)] <- as.character(memb_sub)
  }

  # Label isolates
  if (length(isolates)) {
    iso_names <- vnames[isolates]
    iso_names <- sort(iso_names, method = "radix")
    full_memb[iso_names] <- paste0("Singleton_", seq_along(iso_names))
  }
  # Node to obs expansion
  dfm <- getDfMapper(x)
  nodes_cluster <- do.call(rbind, lapply(seq_len(nrow(dfm)), function(i) {
    node <- rownames(dfm)[i]
    raw  <- dfm[i, 1]
    # support character or list-of-character
    obs  <- if (is.list(raw)) unlist(raw[[1]], use.names = FALSE) else
      strsplit(raw, "\\s+")[[1]]
    obs  <- trimws(obs)
    data.frame(node = node, obs = obs, stringsAsFactors = FALSE)
  }))
  nodes_cluster$cluster <- full_memb[nodes_cluster$node]

  # If k = FALSE then concatenate clusters per obs
  if (identical(k, FALSE)) {
    oc <- split(nodes_cluster$cluster, nodes_cluster$obs)
    lab <- vapply(oc, function(cl) paste(sort(unique(cl)), collapse = "_"),
                  character(1))
    obs_cluster <- data.frame(obs = names(lab), cluster = unname(lab),
                              stringsAsFactors = FALSE)
  } else {
    # k-NN tie-break
    if (!is.numeric(k) || length(k) != 1 || k < 1)
      stop("`k` must be integer >= 1 or FALSE")
    k <- as.integer(k)
    dist_mat <- as.matrix(getDistMat(x))
    if (is.null(rownames(dist_mat)) || is.null(colnames(dist_mat))) {
      stop("Distance matrix must have row/column names.")
    }
    all_obs <- unique(nodes_cluster$obs)
    if (!all(all_obs %in% rownames(dist_mat))) {
      missing <- setdiff(all_obs, rownames(dist_mat))
      stop("Some observations have no distances: ",
           paste(head(missing, 10), collapse = ", "),
           if (length(missing) > 10) " ...")
    }

    chooseCluster <- function(obs, candidates) {
      d <- dist_mat[obs, ]
      means <- vapply(candidates, function(c) {
        others <- setdiff(nodes_cluster$obs[nodes_cluster$cluster == c], obs)
        if (!length(others)) return(Inf)
        mean(utils::head(sort(d[others]), k))
      }, numeric(1))
      candidates[which.min(means)]
    }

    oc_list <- split(nodes_cluster$cluster, nodes_cluster$obs)
    obs_cluster <- do.call(rbind, lapply(names(oc_list), function(obs) {
      cl <- unique(oc_list[[obs]])
      final <- if (length(cl) == 1L) cl else chooseCluster(obs, cl)
      data.frame(obs = obs, cluster = final, stringsAsFactors = FALSE)
    }))
  }

  x@clustering <- list(
    nodes_cluster = nodes_cluster,
    obs_cluster   = obs_cluster
  )
  invisible(x)
}
