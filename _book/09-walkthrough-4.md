---
title: 'Education Dataset Analysis Pipeline: Walkthrough #4'
output:
  pdf_document: default
  html_document: default
---



## Background

Relationships are important to us. In the case of many research techniques, relationships are&mdash;deservedly&mdash;the focus of analyses. It is not very difficult to imagine *qualitative* techniques to study relationships: one could ask other individuals about who their friends are, why they are their friends, and what they like to do when with them. 

Increasingly, it is also not hard to imagine *quantitative* techniques to study relationships, too. In a way, the same questions that could be used qualitatively can serve as the basis for the quantitative study of relationships. Indeed, **social network analysis** uses these relations in a range of visualizations as well as statistical models.

## Creating a network graph

Now, we'll use the **igraph** package to create a *graph* of our simulated data's network. Once we have this graph, we can calculate many useful descriptive statistics associated with the network's features:  

- *Diameter*: The length of the longest geodesic; that is, the max distance between two nodes in the network.  
- *Density*: The ratio of the number of edges and the number of possible edges for a network of that size.  
- *Transitivity*: The balance of connections. Also called the clustering coefficient. The probability that the adjacent vertices of a vertex are connected. When the clustering coefficient is large it implies that a graph is highly clustered around a few nodes; when it is low it implies that the links in the graph are relatively evenly spread among all the nodes (Hogan, 2017).  
- *Reciprocity*: The proportion of mutual connections (in a directed network). The probability that the opposite counterpart of a directed edge is also included in the graph.  
- *Degree*: The number of connections someone has with others nominating or being a nominee (Kadushin, 2012).

Keep in mind that measures such as diameter and density “can be misleading when comparing graphs of substantially different sizes” (Hogan, 2017, p. 255). Therefore, these measures should only be used to compare networks of similar size, or the same network at different points in time.

First, let's create the network graph, called `sim_graph`:


```r
library(igraph)
```

```
## Error in library(igraph): there is no package called 'igraph'
```

```r
sim_graph <- data %>% select(nominator, nominee) %>%  # this creates an edgelist
    as.matrix %>% 
    graph_from_edgelist(directed=TRUE) %>%
    set_vertex_attr(name='degree', value=degree(., mode='total', loops=FALSE)) %>% 
    set_vertex_attr(name='in_degree', value=degree(., mode='in', loops=FALSE)) %>% 
    set_vertex_attr(name='out_degree', value=degree(., mode='out', loops=FALSE))
```

```
## Error in data %>% select(nominator, nominee) %>% as.matrix %>% graph_from_edgelist(directed = TRUE) %>% : could not find function "%>%"
```

```r
network_summary <- sim_graph %>% V %>% length %>%  # number of vertices/nodes
    rbind(sim_graph %>% gsize) %>%  # number of edges
    rbind(sim_graph %>% diameter) %>%  # max distance between two vertices
    rbind({sim_graph %>% edge_density * 100} %>% round(2)) %>%     
    rbind({sim_graph %>% transitivity("global") * 100} %>% round(2)) %>% 
    rbind({sim_graph %>% reciprocity * 100} %>% round(2))  %>%
    rbind(sim_graph %>% vertex_attr('degree') %>% mean %>% round(2)) %>% 
    rbind(sim_graph %>% vertex_attr('degree') %>% sd %>% round(2)) %>% 
    rbind(sim_graph %>% vertex_attr('degree') %>% median) %>% 
    rbind(sim_graph %>% vertex_attr('degree') %>% min) %>%
    rbind(sim_graph %>% vertex_attr('degree') %>% max)
```

```
## Error in sim_graph %>% V %>% length %>% rbind(sim_graph %>% gsize) %>% : could not find function "%>%"
```

```r
colnames(network_summary) <- c("")
```

```
## Error in colnames(network_summary) <- c(""): object 'network_summary' not found
```

```r
rownames(network_summary) <- c("Number of nodes: ", "Number of edges: ", "Diameter: ",
                               "Density: ", "Transitivity: ", "Reciprocity: ",
                               "Mean degree: ", "SD degree: ", "Median degree: ",
                               "Min degree: ", "Max degree: ")
```

```
## Error in rownames(network_summary) <- c("Number of nodes: ", "Number of edges: ", : object 'network_summary' not found
```

```r
network_summary
```

```
## Error in eval(expr, envir, enclos): object 'network_summary' not found
```

## Clustering

With a graph of our simulated network, we can see if *clustering* occurs in this network and describe some characteristics of these clusters.  

First, a definition: a *cluster*&mdash;also called a community or group&mdash;is a set of nodes with many edges inside the community and few edges between outside it (i.e. between the community itself and the rest of the graph). There are numerous methods for determining network clusters, but here we use the *spinglass clustering algorithm*, which maps community detection onto finding the ground state of an infinite range spin glass. Csardi, Nepusz, and Airoldi (2016, pp. 132-133) explained:

>
The clustering method of [Reichardt and Bornholdt](https://arxiv.org/abs/cond-mat/0603718) (2006) is motivated by spin glass models from statistical physics. Such models are used to describe and explain magnetism at the microscopic scale at finite temperatures. Reichardt and Bornholdt (2006) drew an analogy between spin glass models and the problem of community detection on graphs and proposed an algorithm based on the simulated annealing of the spin glass model to obtain well-defined communities in a graph. A spin glass model consists of a set of particles called spins that are coupled by ferromagnetic or antiferromagnetic bonds. Each spin can be in one of k possible states. The well-known Potts model then defines the total energy of the spin glass with a given spin configuration... Spins and interactions in the Potts model are very similar to graphs: each spin in the model corresponds to a vertex, and each interaction corresponds to an edge... Reichardt and Bornholdt (2006) gave efficient update rules for the above energy function, making it possible to apply a simulated annealing procedure to find the ground state of the model that corresponds to a low energy configuration. Their algorithm starts from a random configuration of spins and tries to flip all the spins once in each time step. After each individual spin flip, the energy of the new configuration is evaluated.
>

In other words, the spinglass clustering algorithm partitions nodes into communities by optimizing an energy function. The energy is optimized using the following function (Reichardt and Bornholdt, 2008): 
$$H({\sigma}) = -\sum(a_{ij} \textrm{internal links}) + \sum(b_{ij}\textrm{internal non-links}) + \sum(c_{ij}\textrm{external links}) - \sum(d_{ij}\textrm{external non-links})$$. 

This function penalizes missing edges or non-links of people/nodes in the same cluster and present links or edges between people/nodes in different clusters. It also rewards present links or edges between people/nodes in the same cluster and missing links or edges between people/nodes in different clusters. Thus, a lower score (i.e., lower energy level) is better as it means that the internal links and external non-links have more weightage in that model. In other words, in a strong model, members within clusters are strongly linked and members in separate clusters are weakly linked . Here, $a_{ij}, b_{ij}, c_{ij}, d_{ij}$ represent the individual weights of the four components. 

The initial R code to produce spinglass clusters is straightforward. First, we identify one "giant" cluster&mdash;basically, all nodes that have even a loose connection to each other and create a new network graph by removing any nodes that are not part of the giant cluster. Let's call the graph of the giant cluster `giant_cluster`:


```r
library(igraph)
```

```
## Error in library(igraph): there is no package called 'igraph'
```

```r
giant_cluster <- sim_graph %>% 
    set_vertex_attr(name='membership', 
                    value = clusters(sim_graph) %>% as.data.frame %>% pull(membership)
    ) %>%
    delete_vertices({vertex_attr(., 'membership') != 1} %>% which)
```

```
## Error in sim_graph %>% set_vertex_attr(name = "membership", value = clusters(sim_graph) %>% : could not find function "%>%"
```

Next, we separate the giant cluster into more meaningful clusters, as determined by the spinglass clustering algorithm. We store the cluster information in `csg0`:


```r
<<<<<<< HEAD
library(tidyverse)

create_yvar2 = function(yv1) {
  # Creates yvar2 as a linear outcome of yvar1 
  # Args: 
  #   yv1: yvar1
  yv1 * .25 + rnorm(n = 1, mean = 10, sd = 100)
}

#------------------------------------------------------------------------------

# Make the dataset 
high_y_values <- tibble(
  yvar1 = sample(1000:10000, 100, replace = TRUE)) %>% 
  mutate(yvar2 = map_dbl(yvar1, create_yvar2))
```

```
## Error: <text>:1:1: unexpected input
## 1: <<
##     ^
```


```r
library(simstudy)
=======
csg_0 <- giant_cluster %>% cluster_spinglass  # creates the clusters; 'csg' = cluster spinglass
csg_0$membership %>% unique %>% length  # number of clusters/communities/groups
```

```
## Error: <text>:2:1: unexpected '=='
## 1: library(simstudy)
## 2: ==
##    ^
```
>>>>>>> 58d1f3d32df76c709efd19069c532b125f465a3d

One of the important outcomes of this method is the _modularity_ value $M$. Modularity measures how good the division is, or how separated are the different vertex types from each other. The spinglass algorithm looks for the modularity of the optimal partition. For a given network, the partition with maximum modularity corresponds to the optimal community structure (i.e., a higher $M$ is better).

The maximum modularity score is +1 and according to Hogan (2017), networks with a modularity above 0.3 as "very modular," meanig that most edges are within communities. Note also that if $M$ = 0, all nodes belong to one group.

<<<<<<< HEAD
data1 <- genData(500, def)
data1
=======

```r
csg_0$modularity
>>>>>>> 58d1f3d32df76c709efd19069c532b125f465a3d
```

```
## Error: <text>:2:1: unexpected '>'
## 1: csg_0$modularity
## 2: >
##    ^
```

### Identifying the "typical" number of clusters returned with the spinglass algorithm

It is important to note that a different result is returned each time the spinglass clustering algorithm is run. For this reason, we needed to run a number of simulations to see what the "typical" number of clusters are. We ran the algorithm 100 times and looked at the mean and median number of clusters obtained. We made a note of a _seed_ that produced the median number of clusters, confirmed that this was reproducible, and then set this seed so that all future work will be run with this same clustering configuration.


```r
csg_matrix <- matrix(NA, nrow=1, ncol=100)
for (i in 1:100) {
    print(i)
    set.seed(i)
    csg = giant_cluster %>% cluster_spinglass
    csg_matrix[1,i] <- max(csg$membership)
}
```

```
## [1] 1
```

```
## Error in giant_cluster %>% cluster_spinglass: could not find function "%>%"
```

```r
csg_matrix_summary <- csg_matrix %>% length %>%
    rbind(csg_matrix %>% mean %>% round(2)) %>% 
    rbind(csg_matrix %>% sd %>% round(2)) %>% 
    rbind(csg_matrix %>% median) %>% 
    rbind(csg_matrix %>% min) %>%
    rbind(csg_matrix %>% max)
```

```
## Error in csg_matrix %>% length %>% rbind(csg_matrix %>% mean %>% round(2)) %>% : could not find function "%>%"
```

```r
colnames(csg_matrix_summary) <- c("")
```

```
## Error in colnames(csg_matrix_summary) <- c(""): object 'csg_matrix_summary' not found
```

```r
rownames(csg_matrix_summary) <- c("number of tests: ", "mean: ", "sd: ",  "median: ", "min: ", "max: ")
```

```
## Error in rownames(csg_matrix_summary) <- c("number of tests: ", "mean: ", : object 'csg_matrix_summary' not found
```

```r
csg_matrix_summary
```

```
## Error in eval(expr, envir, enclos): object 'csg_matrix_summary' not found
```


```r
## select a seed from this list which reproduces the median number of clusters
seeds <-{as.vector(csg_matrix) == median(csg_matrix)} %>% which
```

```
## Error in {: could not find function "%>%"
```

```r
our_seed <- seeds[1]
```

```
## Error in eval(expr, envir, enclos): object 'seeds' not found
```

```r
set.seed(our_seed)  # set the seed
```

```
## Error in set.seed(our_seed): object 'our_seed' not found
```

```r
csg <- giant_cluster %>% cluster_spinglass
```

```
## Error in giant_cluster %>% cluster_spinglass: could not find function "%>%"
```

```r
csg_summary <- csg$vcount %>% 
    rbind(giant_cluster %>% gsize) %>% 
    rbind(csg$csize %>% length) %>% 
    rbind(csg$modularity %>% round(4))
```

```
## Error in csg$vcount %>% rbind(giant_cluster %>% gsize) %>% rbind(csg$csize %>% : could not find function "%>%"
```

```r
colnames(csg_summary) <- c("")
```

```
## Error in colnames(csg_summary) <- c(""): object 'csg_summary' not found
```

```r
rownames(csg_summary) <- c("Number of nodes: ", "Number of edges: ", "Number of clusters: ", "Modularity: ")
```

```
## Error in rownames(csg_summary) <- c("Number of nodes: ", "Number of edges: ", : object 'csg_summary' not found
```

```r
csg_summary
```

```
## Error in eval(expr, envir, enclos): object 'csg_summary' not found
```

```r
print("Size of each cluster: ", quote=FALSE); print(csg$csize)
```

```
## [1] Size of each cluster:
```

```
## Error in print(csg$csize): object 'csg' not found
```

### Test of statistical significance for spinglass clusters

The test for statistical significance for spinglass clustering is a bit different than the familiar tests that return $p$-values (Csardi, Nepusz, & Airoldi (2016, pp. 132-138).

The idea behind this test of significance is that a random network of equal size and degree distribution as our observed network should have a lower modularity score--that is, if the observed network does in fact have statistically significant clustering.

The following R procedure generates 100 randomized instances of our network (with the same size and degree distribution) using the `sample_ degseq()` function. The `method = 'vl'` ensures that there are no loop edges in the randomly generated networks. We then applied the spinglass clustering algorithm to each of the 100 randomized instances of the network.

A '0' result from this procudure indicates that no randomized networks have community structure with a modularity score that is higher than the one obtained from the original, observed network. Hence a '0' result means that our network has significant community structure; any non-zero results means that the detected spinglass clusters are not statistically significant.


```r
degrees <- giant_cluster %>% as.undirected %>% degree(mode='all', loops=FALSE)
```

```
## Error in giant_cluster %>% as.undirected %>% degree(mode = "all", loops = FALSE): could not find function "%>%"
```

```r
qr_vl <- replicate(100, sample_degseq(degrees, method="vl"), 
                   simplify=FALSE) %>%
    lapply(cluster_spinglass) %>%
    sapply(modularity) 
```

```
## Error in replicate(100, sample_degseq(degrees, method = "vl"), simplify = FALSE) %>% : could not find function "%>%"
```

```r
sum(qr_vl > csg$modularity) / 100
```

```
## Error in eval(expr, envir, enclos): object 'qr_vl' not found
```

## Network Visualization with Clusters

Visualizations of social networks are interesting and powerful--and increasingly common.

Here, we create a visualization of our network structure, using the color palette generated by our spinglass clustering.


```r
## color-blind palette
## source: https://jacksonlab.agronomy.wisc.edu/2016/05/23/15-level-colorblind-friendly-palette/
palette <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
             "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
             "#920000","#924900","#db6d00","#24ff24","#ffff6d")
csg_palette <- palette[csg$membership]
```

```
## Error in eval(expr, envir, enclos): object 'csg' not found
```

Here, we used the _Fruchterman-Reingold layout algorithm_ (`layout = 'fr'`), which is appropriate for large (but still with less than 1,000 nodes), potentially disconnected networks.


```r
library(ggraph)
```

```
## Error in library(ggraph): there is no package called 'ggraph'
```

```r
library(gridExtra)
```

```
## Error in library(gridExtra): there is no package called 'gridExtra'
```

```r
layout_randomly <- giant_cluster %>% create_layout(layout='randomly')
```

```
## Error in giant_cluster %>% create_layout(layout = "randomly"): could not find function "%>%"
```

```r
layout_mds <- giant_cluster %>% create_layout(layout='mds')
```

```
## Error in giant_cluster %>% create_layout(layout = "mds"): could not find function "%>%"
```

```r
layout_fr <- giant_cluster %>% create_layout(layout='fr')
```

```
## Error in giant_cluster %>% create_layout(layout = "fr"): could not find function "%>%"
```

```r
layout_drl <- giant_cluster %>% create_layout(layout='drl')
```

```
## Error in giant_cluster %>% create_layout(layout = "drl"): could not find function "%>%"
```

```r
layout_kk <- giant_cluster %>% create_layout(layout='kk')
```

```
## Error in giant_cluster %>% create_layout(layout = "kk"): could not find function "%>%"
```

```r
layout_sugiyama <- giant_cluster %>% create_layout(layout='sugiyama')
```

```
## Error in giant_cluster %>% create_layout(layout = "sugiyama"): could not find function "%>%"
```

```r
## Additional layout algorithms to try:
#layout_auto <- giant_cluster %>% create_layout(layout='nicely')
#layout_lgl <- giant_cluster %>% create_layout(layout='lgl')
#layout_dh <- giant_cluster %>% create_layout(layout='dh')
#layout_graphopt <- giant_cluster %>% create_layout(layout='graphopt')

viz_random <- ggraph(layout_randomly) +
    geom_edge_link(width=.1, arrow = arrow(length=unit(1, 'mm'))) +
    geom_node_point(alpha=.75, size=4, color=csg_palette) +
    theme_bw() +
    theme(plot.background = element_blank(), panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), legend.position="none"
    )
```

```
## Error in ggraph(layout_randomly): could not find function "ggraph"
```

```r
viz_mds <- ggraph(layout_mds) +
    geom_edge_link(width=.1, arrow = arrow(length=unit(1, 'mm'))) +
    geom_node_point(alpha=.75, size=4, color=csg_palette) +
    theme_bw() +
    theme(plot.background = element_blank(), panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), legend.position="none"
    )
```

```
## Error in ggraph(layout_mds): could not find function "ggraph"
```

```r
viz_fr <- ggraph(layout_fr) +
    geom_edge_link(width=.1, arrow = arrow(length=unit(1, 'mm'))) +
    geom_node_point(alpha=.75, size=4, color=csg_palette) +
    theme_bw() +
    theme(plot.background = element_blank(), panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), legend.position="none"
    )
```

```
## Error in ggraph(layout_fr): could not find function "ggraph"
```

```r
viz_drl <- ggraph(layout_drl) +
    geom_edge_link(width=.1, arrow = arrow(length=unit(1, 'mm'))) +
    geom_node_point(alpha=.75, size=4, color=csg_palette) +
    theme_bw() +
    theme(plot.background = element_blank(), panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), legend.position="none"
    )
```

```
## Error in ggraph(layout_drl): could not find function "ggraph"
```

```r
viz_kk <- ggraph(layout_kk) +
    geom_edge_link(width=.1, arrow = arrow(length=unit(1, 'mm'))) +
    geom_node_point(alpha=.75, size=4, color=csg_palette) +
    theme_bw() +
    theme(plot.background = element_blank(), panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), legend.position="none"
    )
```

```
## Error in ggraph(layout_kk): could not find function "ggraph"
```

```r
viz_sugiyama <- ggraph(layout_sugiyama) +
    geom_edge_link(width=.1, arrow = arrow(length=unit(1, 'mm'))) +
    geom_node_point(alpha=.75, size=4, color=csg_palette) +
    theme_bw() +
    theme(plot.background = element_blank(), panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), legend.position="none"
    )
```

```
## Error in ggraph(layout_sugiyama): could not find function "ggraph"
```

```r
grid.arrange(viz_random, viz_mds, viz_fr,
             viz_drl, viz_kk, viz_sugiyama, 
             nrow = 3)
```

```
## Error in grid.arrange(viz_random, viz_mds, viz_fr, viz_drl, viz_kk, viz_sugiyama, : could not find function "grid.arrange"
```

<<<<<<< HEAD
=======
The layouts above are as follows: the first row has 1.) a *random* plotting of nodes and 2.) the *MDS* algorithm. The second row has 3.) the *FR* algorithm and 4.) the *DRL* algorithm. The third row as 5.) the *KK* algorithm, and 6.) the *Sugiyama* algorithm. Each of these is appropriate in different situations; if unsure of which to use, the *nicely* option automatically selects an algorithm based on the network size.

## Selection and influence

Behind these visualizations, though, there are also statistical models and methods that can help to understand what is going on with respect to particular relationships in a network in additional ways.

One way to consider these models and methods is in terms of two *processes* at play in our relationships (cite). These two processes are commonly (though not exclusively) the focus of statistical analyses of networks. In addition to not being exclusive, they do not interact independently: they affect each other reciprocally (Xu, Frank, & Penuel, 2018). They are:

- *Selection*: the processes regarding who chooses to have a relationship with whom
- *Infuence*: the processes regarding how who we have relationships with affects our behavior

While these are complex, they can be studied with the type of data collected from asking people about their relationships (and possibly asking them about or studying their behavior--or measuring some outcome). Happily, the use of these methods has expanded along with **R**: many of the best tools for studying social networks are in the form of long-standing R packages. Additionally, while there are many potential naunces to studying selection and influence, these are models that can fundamentally be carried out with regression, or the linear model (or extensions of it).

In this walkthrough, the influence model is the focus. Nevertheless, we provide some direction for how to carry out selection modeling, too, at the end. 

## Creating example data in the form of an edgelist 

First, let's create three different data frames. Here is what they should contain:

- A data frame indicating who the *nominator* and *nominee* for the relation (i.e., if Stefanie says that José is her friend, then Stefanie is the nominator and José the nominee) - as well as an optional variable indicating the weight, or strength, of their relation.
- This data frame and its type can be considered the basis for many types of social network analysis and is a common structure for network data: it is an *edgelist*.
- Data frames indicating the values of some behavior - an outcome - at two different time points.

## An example of influence

In this example, we create some example data that can be used to explore questions about how influence works. Note that Joshua Rosenberg and Sarah Galey initially wrote the following code for a walkthrough shared on Ken Frank's website [here](https://msu.edu/~kenfrank/resources.htm).

>>>>>>> 58d1f3d32df76c709efd19069c532b125f465a3d
Let's take a look at the merged data. What this data now contains is the first data frame, `data1`, with each nominees' outcome at time 1 (`yvar1`). Note that we will find each nominators' outcome at time 2 later on.


```r
data <- as_tibble(data)
```

```
## Error in as_tibble(data): could not find function "as_tibble"
```

```r
data
```

```
## function (..., list = character(), package = NULL, lib.loc = NULL, 
##     verbose = getOption("verbose"), envir = .GlobalEnv, overwrite = TRUE) 
## {
##     fileExt <- function(x) {
##         db <- grepl("\\.[^.]+\\.(gz|bz2|xz)$", x)
##         ans <- sub(".*\\.", "", x)
##         ans[db] <- sub(".*\\.([^.]+\\.)(gz|bz2|xz)$", "\\1\\2", 
##             x[db])
##         ans
##     }
##     names <- c(as.character(substitute(list(...))[-1L]), list)
##     if (!is.null(package)) {
##         if (!is.character(package)) 
##             stop("'package' must be a character string or NULL")
##         if (any(package %in% "base")) 
##             warning("datasets have been moved from package 'base' to package 'datasets'")
##         if (any(package %in% "stats")) 
##             warning("datasets have been moved from package 'stats' to package 'datasets'")
##         package[package %in% c("base", "stats")] <- "datasets"
##     }
##     paths <- find.package(package, lib.loc, verbose = verbose)
##     if (is.null(lib.loc)) 
##         paths <- c(path.package(package, TRUE), if (!length(package)) getwd(), 
##             paths)
##     paths <- unique(normalizePath(paths[file.exists(paths)]))
##     paths <- paths[dir.exists(file.path(paths, "data"))]
##     dataExts <- tools:::.make_file_exts("data")
##     if (length(names) == 0L) {
##         db <- matrix(character(), nrow = 0L, ncol = 4L)
##         for (path in paths) {
##             entries <- NULL
##             packageName <- if (file_test("-f", file.path(path, 
##                 "DESCRIPTION"))) 
##                 basename(path)
##             else "."
##             if (file_test("-f", INDEX <- file.path(path, "Meta", 
##                 "data.rds"))) {
##                 entries <- readRDS(INDEX)
##             }
##             else {
##                 dataDir <- file.path(path, "data")
##                 entries <- tools::list_files_with_type(dataDir, 
##                   "data")
##                 if (length(entries)) {
##                   entries <- unique(tools::file_path_sans_ext(basename(entries)))
##                   entries <- cbind(entries, "")
##                 }
##             }
##             if (NROW(entries)) {
##                 if (is.matrix(entries) && ncol(entries) == 2L) 
##                   db <- rbind(db, cbind(packageName, dirname(path), 
##                     entries))
##                 else warning(gettextf("data index for package %s is invalid and will be ignored", 
##                   sQuote(packageName)), domain = NA, call. = FALSE)
##             }
##         }
##         colnames(db) <- c("Package", "LibPath", "Item", "Title")
##         footer <- if (missing(package)) 
##             paste0("Use ", sQuote(paste("data(package =", ".packages(all.available = TRUE))")), 
##                 "\n", "to list the data sets in all *available* packages.")
##         else NULL
##         y <- list(title = "Data sets", header = NULL, results = db, 
##             footer = footer)
##         class(y) <- "packageIQR"
##         return(y)
##     }
##     paths <- file.path(paths, "data")
##     for (name in names) {
##         found <- FALSE
##         for (p in paths) {
##             tmp_env <- if (overwrite) 
##                 envir
##             else new.env()
##             if (file_test("-f", file.path(p, "Rdata.rds"))) {
##                 rds <- readRDS(file.path(p, "Rdata.rds"))
##                 if (name %in% names(rds)) {
##                   found <- TRUE
##                   if (verbose) 
##                     message(sprintf("name=%s:\t found in Rdata.rds", 
##                       name), domain = NA)
##                   thispkg <- sub(".*/([^/]*)/data$", "\\1", p)
##                   thispkg <- sub("_.*$", "", thispkg)
##                   thispkg <- paste0("package:", thispkg)
##                   objs <- rds[[name]]
##                   lazyLoad(file.path(p, "Rdata"), envir = tmp_env, 
##                     filter = function(x) x %in% objs)
##                   break
##                 }
##                 else if (verbose) 
##                   message(sprintf("name=%s:\t NOT found in names() of Rdata.rds, i.e.,\n\t%s\n", 
##                     name, paste(names(rds), collapse = ",")), 
##                     domain = NA)
##             }
##             if (file_test("-f", file.path(p, "Rdata.zip"))) {
##                 warning("zipped data found for package ", sQuote(basename(dirname(p))), 
##                   ".\nThat is defunct, so please re-install the package.", 
##                   domain = NA)
##                 if (file_test("-f", fp <- file.path(p, "filelist"))) 
##                   files <- file.path(p, scan(fp, what = "", quiet = TRUE))
##                 else {
##                   warning(gettextf("file 'filelist' is missing for directory %s", 
##                     sQuote(p)), domain = NA)
##                   next
##                 }
##             }
##             else {
##                 files <- list.files(p, full.names = TRUE)
##             }
##             files <- files[grep(name, files, fixed = TRUE)]
##             if (length(files) > 1L) {
##                 o <- match(fileExt(files), dataExts, nomatch = 100L)
##                 paths0 <- dirname(files)
##                 paths0 <- factor(paths0, levels = unique(paths0))
##                 files <- files[order(paths0, o)]
##             }
##             if (length(files)) {
##                 for (file in files) {
##                   if (verbose) 
##                     message("name=", name, ":\t file= ...", .Platform$file.sep, 
##                       basename(file), "::\t", appendLF = FALSE, 
##                       domain = NA)
##                   ext <- fileExt(file)
##                   if (basename(file) != paste0(name, ".", ext)) 
##                     found <- FALSE
##                   else {
##                     found <- TRUE
##                     zfile <- file
##                     zipname <- file.path(dirname(file), "Rdata.zip")
##                     if (file.exists(zipname)) {
##                       Rdatadir <- tempfile("Rdata")
##                       dir.create(Rdatadir, showWarnings = FALSE)
##                       topic <- basename(file)
##                       rc <- .External(C_unzip, zipname, topic, 
##                         Rdatadir, FALSE, TRUE, FALSE, FALSE)
##                       if (rc == 0L) 
##                         zfile <- file.path(Rdatadir, topic)
##                     }
##                     if (zfile != file) 
##                       on.exit(unlink(zfile))
##                     switch(ext, R = , r = {
##                       library("utils")
##                       sys.source(zfile, chdir = TRUE, envir = tmp_env)
##                     }, RData = , rdata = , rda = load(zfile, 
##                       envir = tmp_env), TXT = , txt = , tab = , 
##                       tab.gz = , tab.bz2 = , tab.xz = , txt.gz = , 
##                       txt.bz2 = , txt.xz = assign(name, read.table(zfile, 
##                         header = TRUE, as.is = FALSE), envir = tmp_env), 
##                       CSV = , csv = , csv.gz = , csv.bz2 = , 
##                       csv.xz = assign(name, read.table(zfile, 
##                         header = TRUE, sep = ";", as.is = FALSE), 
##                         envir = tmp_env), found <- FALSE)
##                   }
##                   if (found) 
##                     break
##                 }
##                 if (verbose) 
##                   message(if (!found) 
##                     "*NOT* ", "found", domain = NA)
##             }
##             if (found) 
##                 break
##         }
##         if (!found) {
##             warning(gettextf("data set %s not found", sQuote(name)), 
##                 domain = NA)
##         }
##         else if (!overwrite) {
##             for (o in ls(envir = tmp_env, all.names = TRUE)) {
##                 if (exists(o, envir = envir, inherits = FALSE)) 
##                   warning(gettextf("an object named %s already exists and will not be overwritten", 
##                     sQuote(o)))
##                 else assign(o, get(o, envir = tmp_env, inherits = FALSE), 
##                   envir = envir)
##             }
##             rm(tmp_env)
##         }
##     }
##     invisible(names)
## }
## <bytecode: 0x7fc572bbfac8>
## <environment: namespace:utils>
```

### Calculating an exposure term

This is the key step that makes this model - a regression, or linear, model - one that is special. It is creating an exposure term. The idea is that the exposure term "captures" how your interactions with someone, over some period of time (between the first and second time points) impact some outcome. This model accounts for an individual's initial report of the outcome, i.e., their time 1 prior value, so it is a model for *change* in some outcome.


```r
# Calculating exposure
data$exposure <- data$relate * data$yvar1

# Calculating mean exposure
mean_exposure <- data %>%
    group_by(nominator) %>%
    summarize(exposure_mean = mean(exposure))
<<<<<<< HEAD

=======
>>>>>>> 58d1f3d32df76c709efd19069c532b125f465a3d
```

```
## Error: <text>:8:1: unexpected input
## 7:     summarize(exposure_mean = mean(exposure))
## 8: <<
##    ^
```

What this data frame - `mean_exposure` - contains is the mean of the outcome (in this case, `yvar1`) for all of the individuals the nominator had a relation with.

As we need a final data set with `mean_exposure`, `mean_exposure_plus`, `degree`, `yvar1`, and `yvar2` added, we'll process the data a bit more.


```r
mean_exposure_terms <- left_join(mean_exposure, mean_exposure_plus, by = "nominator")
```

```
## Error in left_join(mean_exposure, mean_exposure_plus, by = "nominator"): could not find function "left_join"
```

```r
names(data2) <- c("nominator", "yvar1") # rename nominee as nominator to merge these
```

```
## Error in names(data2) <- c("nominator", "yvar1"): object 'data2' not found
```

```r
final_data <- left_join(mean_exposure_terms, data2, by = "nominator")
```

```
## Error in left_join(mean_exposure_terms, data2, by = "nominator"): could not find function "left_join"
```

```r
final_data <- left_join(final_data, data3, by = "nominator") # data3 already has nominator, so no need to change
```

```
## Error in left_join(final_data, data3, by = "nominator"): could not find function "left_join"
```

### Regression (linear models)

Calculating the exposure term is the most distinctive and important step in carrying out influence models. Now, we can simply use a linear model to find out how much relations - as captured by the influence term - affect some outcome.


```r
model1 <- lm(yvar2 ~ yvar1 + exposure_mean, data = final_data)
```

```
## Error in is.data.frame(data): object 'final_data' not found
```

```r
summary(model1)
```

```
## Error in summary(model1): object 'model1' not found
```

Note that these models show ... [add]

So, the influence model is used to study a key process for social network analysis, but it is one that is useful, because you can quantify, given what you measure and how you measure it, *the network effect*, something that is sometimes not considered, especially in education (Sweet, 2017). It's also fundamentally a regression. That's really it, as the majority of the work goes into calculating the exposure term.

## Selection models

While this tutorial focused on influence models, selection models are also commonly used - and are commonly of interest not only to researchers but also to administrators and teachers (and even to youth and students). 

Here, we briefly describe a few possible approaches for using a selection model.

At its core, the selection model is a regression - albeit, one that is a generalization of one, namely, a logistic regression (sometimes termed a generalized linear model, because it is *basically* a regression but is one with an outcome that consists just of 0's and 1's). Thus, the most straight-away way to use a selection model is to use a logistic regression where all of the relations (note the `relate` variable in `data1` above) are indicated with a 1. But, here is the important and challenging step: all of the *possible relations* (i.e., all of the relations that are possible between all of the individuals in a network) are indicated with a 0 in an edgelist. Note that, again, an edgelist is the preferred data structure for carrying out this analysis. This step involves some data wrangling, especially the idea of widening or lengthening a data frame.

<!-- May want to add a short bit of code on this using `gather()` and `spread()` -->

Once all of the relations are indicated with a 1 or a 0, then a simple linear regression can be used. Imagine that we are interested in whether individuals from the *same* group are more or less likely to interact than those from different groups; same could be created in the data frame based upon knowing which group both nominator and nominee are from:


```r
m_selection <- glm(relate ~ 1 + same, data = edgelist1)
```

While this is a straightforward way to carry out a selection model, there are some limitations to it. Namely, it does not account for individuals who send more (or less) nominations overall--and not considering this may mean other effects, like the one associated with being from the *same* group, are not accurate. A few extensions of the linear model - including those that can use data for which relationships are indicated with weights, not just 1's and 0's, have been developed. 

One type of model extends the logistic regression. It can be used for data that is not only 1's and 0's but also data that is normally distributed or has fixed-ranks. It is the **amen** package available [here](https://cran.r-project.org/web/packages/amen/index.html).

A particularly common one is an Exponential Random Graph Model, or an ERGM. An R package that makes estimating these easy is available [here](https://cran.r-project.org/web/packages/ergm/index.html). That R package, **ergm**, is part of a powerful and often-used collection of packages, including those for working with network data (data that can begin with an edgelist, but may need additional processing that is challenging to do with edgelist data), **statnet**. A link to the statnet packages is [here](https://statnet.org/).

## References

Bates, D., Maechler, M., Bolker, B., & Walker, S. (2018). lme4: Linear mixed-effects models using 'Eigen' and S4 (Version 1.1-19) [R package]. Retrieved from https://cran.r-project.org/package=lme4

Csardi, G. (2018). igraph: Network analysis and visualization (Version 1.2.2) [R package]. Retrieved from https://CRAN.R-project.org/package=igraph

Csardi, G., Nepusz, T., & Airoldi, E. M. (2016). *Statistical network analysis with igraph*. New York, NY: Springer. Retrieved from https://sites.fas.harvard.edu/~airoldi/pub/books/BookDraft-CsardiNepuszAiroldi2016.pdf

Hogan, B. (2017). Online social networks: Concepts for data collection and analysis. In N. G. Fielding, R. M. Lee, & G. Blank (Eds.), *The SAGE handbook of online research methods* (2nd ed., pp. ). London, UK: SAGE.

Kadushin, C. (2012). *Understanding social networks: Theories, concepts, and findings.* New York, NY: Oxford University Press.

Pedersen, T. L. (2018). ggraph: An implementation of grammar of graphics for graphs and networks (Version 1.0.2) [R package]. Retrieved from https://CRAN.R-project.org/package=ggraph

R Core Team. (2018). R: A language and environment for statistical computing (Version 3.5.0) [Computer software]. Vienna, Austria: R Foundation for Statistical Computing. Retrieved from https://www.R-project.org/

Reichardt, J., & Bornholdt, S. (2006). Statistical mechanics of community detection. *Physical Review E, 74*(1), 016110. Retrieved from https://arxiv.org/abs/cond-mat/0603718

Wickham, H., Chang, W., & RStudio. (2016). ggplot2: Create elegant data visualisations using the grammar of graphics (Version 2.2.1) [R package]. Retrieved from https://CRAN.R-project.org/package=ggplot2

Wickham, H., Francois, R., Henry, L., Muller, K., & RStudio. (2018). dplyr: A grammar of data manipulation (Version 0.7.6) [R package]. Retrieved from https://CRAN.R-project.org/package=dplyr
