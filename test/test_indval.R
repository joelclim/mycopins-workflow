library(readr)
library(labdsv)
library(vegan)

working_folder_path <- "C:/Users/Joel/work/kean-stme-2903-11/mycopins/shiny-apps/mycopins/"
data_file <- c(
  paste0(working_folder_path, "data/FA23_complete_counts.csv"),
  paste0(working_folder_path, "data/SP24_complete_counts.csv"),
  paste0(working_folder_path, "data/transectC_complete_counts.csv"),
  paste0(working_folder_path, "data/FA23_complete_counts_fungi.csv"),
  paste0(working_folder_path, "data/SP24_complete_counts_fungi.csv"),
  paste0(working_folder_path, "data/transectC_counts_fungi.csv"),
  paste0(working_folder_path, "data/FA23_complete_counts_wood_saprotroph.csv"),
  paste0(working_folder_path, "data/SP24_complete_counts_wood_saprotroph.csv"),
  paste0(working_folder_path, "data/transectC_counts_wood_saprotroph.csv")
)

get_counts_data <- function(team, species) {
  index <- as.numeric(team) + (as.numeric(species) - 1) * 3
  
  return(read_csv(data_file[index]))
}

###############################################################################

team <- 3 # FA23 + SP24
species <- 1 # All species
df <- get_counts_data(team, species)

discriminator_index <- which(names(df) == "_Splitter_")

wood <- df[, (discriminator_index+1):ncol(df)]
wood.env <- df[, 1:(discriminator_index-1)]

wood.groups <- wood.env$Wood.Texture

indval_results <- indval(wood, wood.groups, permutations = 1000)

# relative frequency of species in classes
relfrq_df <- data.frame(indval_results$relfrq)

# relative abundance of species in classes
relabu_df <- data.frame(indval_results$relabu)

# the indicator value for each species
indval_df <- data.frame(indval_results$indval)

# Indicator Species Analysis Results
# RelAbu, pValue, Class
# the class each species has maximum indicator value for
maxcls_df <- data.frame(indval_results$maxcls)

# the indicator value for each species to its maximum class
indcls_df <- data.frame(indval_results$indcls)

# the probability of obtaining as high an indicator values as observed over the specified iterations
pval_df <- data.frame(indval_results$pval)

indval_df <- data.frame(indval_results$indval)
groups <- names(indval_df)[indval_results$maxcls]

values <- c()
for (row in 1:length(indval_results$maxcls)) {
  col <- indval_results$maxcls[row]
  values <- c(values, indval_results$relfrq[row, col])
}
frqs <- data.frame(Frequency=values)

values <- c()
for (row in 1:length(indval_results$maxcls)) {
  col <- indval_results$maxcls[row]
  values <- c(values, indval_results$relabu[row, col])
}
abus <- data.frame(Abundance=values)

results_df <- data.frame(Group=groups, # group
                         frqs,
                         abus,
                         Indicator.Value = indval_results$indcls, # indicator value 
                         P_value = indval_results$pval    # p-Value
                         )
names(results_df) <- c("Group",
                       
                       "Indicator Value", 
                       "P-Value")

method <- "bray"
wood.dist <- vegdist(wood, method = method)

# Perform hierarchical clustering
k <- 3
hc <- hclust(wood.dist)
# Cutting the dendrogram to form k clusters
clusters <- cutree(hc, k)
# The clusters vector now contains the cluster membership for each sample
print(clusters)

indval_results <- indval(wood, clusters, permutations = 1000)
