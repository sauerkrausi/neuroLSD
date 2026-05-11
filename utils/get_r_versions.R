packages <- c(
  "Cairo",
  "CalNetExploreR",
  "ComplexHeatmap",
  "ComplexUpset",
  "NatParksPalettes",
  "RColorBrewer",
  "UpSetR",
  "bigstatsr",
  "broom",
  "circlize",
  "cowplot",
  "data.table",
  "devtools",
  "dplyr",
  "factoextra",
  "forcats",
  "furrr",
  "ggVennDiagram",
  "ggbiplot",
  "ggdendro",
  "ggdist",
  "gghalves",
  "gghighlight",
  "ggplot2",
  "ggpmisc",
  "ggpubr",
  "ggraph",
  "ggrepel",
  "ggsignif",
  "gptstudio",
  "grid",
  "gridExtra",
  "igraph",
  "irlba",
  "limma",
  "lintr",
  "magick",
  "org.Hs.eg.db",
  "patchwork",
  "pheatmap",
  "plotly",
  "plyr",
  "png",
  "qs",
  "readr",
  "reshape2",
  "rlang",
  "scales",
  "signal",
  "stringr",
  "superheat",
  "tidyplots",
  "tidyr",
  "tidyverse",
  "topGO",
  "umap",
  "uwot",
  "viridis",
  "zoo"
)

versions <- sapply(packages, function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    as.character(packageVersion(pkg))
  } else {
    NA
  }
})

versions_df <- data.frame(
  Package = names(versions),
  Version = unname(versions),
  row.names = NULL
)
write.csv(
  versions_df,
  "/Users/felix/Desktop/r_package_versions.csv",
  row.names = FALSE
)
