library('randomcoloR')
library('colorspace')
library('RColorBrewer')
library('dplyr')

#' Adds colors for taxa labels to metadata.
#' 
#' @param order_files  Path to text file with paths to files with taxa order (=order files).
#' If there are several text paths, they should be separated by \n.
#' The taxa in order files should be separated by \n.
#' 
#' @param metadata Path to text file with metadata for taxa (csv)
#' @return metadata dataframe with additional 'color' column 
add_colors2meta = function(order_files, metadata){
  
  meta = read.csv(metadata)
  
  # Read file with paths to order files
  order_files = read.table(order_files)["V1"]
  list_taxa_df = list()
  for (i in 1:nrow(order_files)){
    taxa = read.table(order_files[i,1])
    colnames(taxa) = "GBAC"
    list_taxa_df = append(list_taxa_df,taxa)
  }
  
  # Generate base colors using RColorBrewer
  base_colors = brewer.pal(max(3, length(list_taxa_df)), "Set3")
  print(base_colors)
  
  for (i in 1:length(list_taxa_df)){
    num_colors = length(list_taxa_df[[i]])
    # Lightened and darkened color range
    color1 = lighten(base_colors[i %% length(base_colors) + 1], 0.3)
    color2 = darken(base_colors[i %% length(base_colors) + 1], 0.3)
    
    print(paste("Base color:", base_colors[i %% length(base_colors) + 1]))
    print(paste("Lightened color:", color1))
    print(paste("Darkened color:", color2))
    
    # Create gradient palette
    color_range <- colorRampPalette(c(color1, color2))
    colors_clade <- color_range(num_colors)
    
    # Create dataframe of colors and corresponding names
    df = as.data.frame(cbind(list_taxa_df[[i]], colors_clade))
    colnames(df) = c("GBAC","color")
    
    # Update list
    list_taxa_df[[i]] = df
  }
  
  # Combine list of dataframes into one dataframe by row
  df_colors = bind_rows(list_taxa_df)
  meta_upd = full_join(meta, df_colors, by="GBAC")
  return(meta_upd)
}
