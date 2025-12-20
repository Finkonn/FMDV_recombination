library('ggtree')
library('dplyr')
library('ggnewscale')
library('phytools')
library('ape')
library('ggplot2')

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("modified_gradients.R")

trees = list.files(path = "all_sequences_treefiles/", full.names = TRUE, pattern = ".treefile$")
meta_path = "metadata_all_sequences_pools.csv"
info = read.csv(meta_path)

unique_serotypes <- unique(info$serotype)
unique_topotypes <- unique(info$Topotype) 
unique_lineages <- unique(info$lineage)
unique_pools <- unique(info$pool)

serotype_colors <- setNames(c("#b0b006","#d84a1c" ,"#43998c", "#857fa9","#9ec744", "#d3973f", "#6194b5"), 
                            unique_serotypes)

lineage_colors <- c(
  "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#321faf", 
  "#F0027F", "#BF5B17", "#842f16", "#1B9E77", "#A1aa00",
  "#7570B3", "#5168EE", "#66A6fE"
)

topotype_colors <- c(
  "#FF7F00", "#A6761D", "#A6CEE3", "#1F78B4", "#c90303",
  "#33A02C", "#FB9A99", "#1111ff", "#FDBF6F", "#E6AB02",
  "#CAB2D6", "#6A3D9A", "#B15928", "#8000FF", "#B3CDE3",
  "#bbE535", "#FFFF00", "#A9A2F9", "#2FaFCC", "#A1D2AD",
  "#FDDAEC", "#F2a2F2", "#a1f8aa", "#FDCDAC", "#1EA5E8",
  "#00FFFF", "#E6F5C9", "#FF0080", "#32CD32", "#fa31f1",
  "#f14e2f"
)
unique_pools

pool_colors <- c("#ffffff", "#fed580ff","#ff8080ff","#e280ffff", "#cdcdcdff", "#7598beff", "#80b8ffff", "#fefe80ff",
                 "#bf5b80ff", "#a5f280ff", "#ca7842ff")

names(lineage_colors) <- unique_lineages
names(topotype_colors) <- unique_topotypes
names(pool_colors) <- unique_pools

plot_tree_with_gradient_and_heatmap = function(tree_file, meta, serotype_colors, topotype_colors, lineage_colors, pool_colors) {
  tree = read.tree(tree_file)
  tree_rooted = midpoint.root(tree)
  info = read.csv(meta)
  
  info$Topotype <- case_when(
    info$Topotype == "EUROPE-SOUTH-AMERICA-EURO-SA-SOUTH" ~ "EURO-SA",
    info$Topotype == "EAST-AFRICA-1-EA-1" ~ "EA-1",
    info$Topotype == "EAST-AFRICA-2-EA-2" ~ "EA-2",
    info$Topotype == "EAST-AFRICA-3-EA-3" ~ "EA-3",
    info$Topotype == "EAST-AFRICA-4-EA-4" ~ "EA-4",
    info$Topotype == "I-NORTHWEST-ZIMBABWE-NWZ" ~ "I-NWZ",
    info$Topotype == "I-SOUTHEAST-ZIMBABWE-SEZ" ~ "I-SEZ",
    info$Topotype == "III-NORTHWEST-ZIMBABWE-NWZ" ~ "III-NWZ",
    info$Topotype == "II-WESTERN-ZIMBABWE-WZ" ~ "II-WZ",
    info$Topotype == "III-WESTERN-ZIMBABWE-WZ" ~ "III-WZ",
    info$Topotype == "INDONESIA-1-ISA-1-1" ~ "ISA-1",
    info$Topotype == "IV-EAST-AFRICA1-EA-1" ~ "IV-EA-1",
    info$Topotype == "MIDDLE-EAST-SOUTH-ASIA-ME-SA" ~ "ME-SA",
    info$Topotype == "V-East-Africa-EA" ~ "V-EA",
    info$Topotype == "VII-EAST-AFRICA-2-EA-2" ~ "VII-EA-2",
    info$Topotype == "WEST-AFRICA-WA" ~ "WA",
    info$Topotype == "SOUTHEAST-ASIA-SEA" ~ "SEA",
    info$Topotype == "EUROPE-SOUTH-AMERICA-EURO-SA" ~ "EURO-SA",
    TRUE ~ info$Topotype
  )

  unique_topotypes_abbreviated <- unique(info$Topotype)
  topotype_colors_abbreviated <- topotype_colors
  names(topotype_colors_abbreviated) <- sapply(names(topotype_colors), function(x) {
    case_when(
      x == "EUROPE-SOUTH-AMERICA-EURO-SA-SOUTH" ~ "EURO-SA",
      x == "EAST-AFRICA-1-EA-1" ~ "EA-1",
      x == "EAST-AFRICA-2-EA-2" ~ "EA-2",
      x == "EAST-AFRICA-3-EA-3" ~ "EA-3",
      x == "EAST-AFRICA-4-EA-4" ~ "EA-4",
      x == "I-NORTHWEST-ZIMBABWE-NWZ" ~ "I-NWZ",
      x == "I-SOUTHEAST-ZIMBABWE-SEZ" ~ "I-SEZ",
      x == "III-NORTHWEST-ZIMBABWE-NWZ" ~ "III-NWZ",
      x == "II-WESTERN-ZIMBABWE-WZ" ~ "II-WZ",
      x == "III-WESTERN-ZIMBABWE-WZ" ~ "III-WZ",
      x == "INDONESIA-1-ISA-1-1" ~ "ISA-1",
      x == "IV-EAST-AFRICA1-EA-1" ~ "IV-EA-1",
      x == "MIDDLE-EAST-SOUTH-ASIA-ME-SA" ~ "ME-SA",
      x == "V-East-Africa-EA" ~ "V-EA",
      x == "VII-EAST-AFRICA-2-EA-2" ~ "VII-EA-2",
      x == "WEST-AFRICA-WA" ~ "WA",
      x == "SOUTHEAST-ASIA-SEA" ~ "SEA",
      x == "EUROPE-SOUTH-AMERICA-EURO-SA" ~ "EURO-SA",
      TRUE ~ x
    )
  })
  
  heat_data = data.frame(
    "serotype" = info$serotype,
    "Topotype" = info$Topotype,
    "Lineage" = info$lineage,
    "Pool" = info$pool
  )
  rownames(heat_data) <- info$GBAC
  
  t = ggtree(tree_rooted, size = 0.3) %<+% info +
    geom_point2(aes(label=label, 
                    subset = !is.na(as.numeric(label)) & as.numeric(label) < 95), size=0.1, color="red",alpha=0.5) +
    geom_tiplab(size = 0.4, aes(color=label)) +
    scale_color_manual(values=info$color, guide='none') + 
    theme(legend.position = "none") 
  
  t_with_serotype <- gheatmap(t, heat_data["serotype"], 
                              offset = 0.5,
                              width = 0.1,
                              color = NULL,
                              colnames = FALSE) +
    scale_fill_manual(values = serotype_colors, name = "Serotype")

  t_with_serotype <- t_with_serotype + ggnewscale::new_scale_fill()

  t_with_topotype <- gheatmap(t_with_serotype, heat_data["Topotype"],
                              offset = 0.7,
                              width = 0.1,
                              color = NULL,
                              colnames = FALSE) +
    scale_fill_manual(values = topotype_colors_abbreviated, name = "Topotype", na.value = "white")

  t_with_topotype <- t_with_topotype + ggnewscale::new_scale_fill()

  t_with_lineage <- gheatmap(t_with_topotype, heat_data["Lineage"],
                             offset = 0.9,
                             width = 0.1,
                             color = NULL,
                             colnames = FALSE) +
    scale_fill_manual(values = lineage_colors, name = "Lineage", na.value = "white")

  t_with_lineage <- t_with_lineage + ggnewscale::new_scale_fill()
  
  t_with_pool <- gheatmap(t_with_lineage, heat_data["Pool"],
                             offset = 1.1,
                             width = 0.1,
                             color = NULL,
                             colnames = FALSE) +
    scale_fill_manual(values = pool_colors, name = "Pool", na.value = "white")

  
  return(t_with_pool)
}


for (file in trees) {
  g = plot_tree_with_gradient_and_heatmap(file, meta_path, serotype_colors, topotype_colors, lineage_colors, pool_colors) +
    ggtitle(strsplit(basename(file), '.treefile')[[1]][1]) +
    guides(
      fill = guide_legend(
        ncol = 1,
        keywidth = unit(0.1, "cm"),
        keyheight = unit(0.1, "cm"),
        override.aes = list(size = 2)
      )
    ) +
    theme(
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 18)
    )
  
  legend <- cowplot::get_legend(g)
  g_nolegend <- g + theme(legend.position = "none")
  
  svg_file <- paste0("Plots/", basename(file), ".svg")
  
  ggsave(svg_file, g_nolegend, height = 16, width = 7, dpi = 600)
  
}

ggsave(paste0("Plots/", basename(file), "_legend.svg"), 
       plot = legend, height = 16, width = 7, dpi = 600)
