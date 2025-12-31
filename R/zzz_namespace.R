#' @keywords internal
"_PACKAGE"

#' @importFrom grDevices rainbow
#' @importFrom methods is new show
#' @importFrom stats approx as.dist cmdscale complete.cases cor cutree dist hclust kmeans prcomp setNames
#' @importFrom utils head
#' @importFrom SummarizedExperiment SummarizedExperiment colData assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom igraph is_igraph is_directed as.undirected simplify degree induced_subgraph membership V
#' @importFrom igraph average.path.length degree_distribution gorder graph_from_adjacency_matrix transitivity fit_power_law
#' @importFrom igraph cluster_fast_greedy cluster_walktrap cluster_edge_betweenness cluster_optimal cluster_label_prop
#' @importFrom cluster daisy
#' @importFrom dbscan dbscan extractDBSCAN kNNdist optics
#' @importFrom umap umap umap.defaults
#' @importFrom tsne tsne
#' @importFrom kernlab kpca
#' @importFrom vegan isomap
#' @importFrom ggplot2 aes annotation_logticks element_text geom_line geom_point ggplot ggtitle labs
#' @importFrom ggplot2 scale_color_manual scale_fill_gradientn scale_fill_manual scale_x_continuous
#' @importFrom ggplot2 scale_x_log10 scale_y_log10 stat_smooth theme theme_bw
#' @importFrom igraph is_igraph is_directed as.undirected induced_subgraph
#' @importFrom igraph components vcount ecount degree transitivity
#' @importFrom igraph average.path.length assortativity_degree
#' @importFrom igraph cluster_fast_greedy modularity membership
#' @importFrom stats var setNames
NULL
