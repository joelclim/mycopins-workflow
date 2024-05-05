library(readr)
library(vegan)
library(viridis)

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

method <- "bray"
group <- 2 # FA23 + SP24
species <- 1 # All species

df <- get_counts_data(group, species)
discriminator_index <- which(names(df) == "_Splitter_")
wood <- df[, (discriminator_index+1):ncol(df)]
wood.env <- df[, 1:(discriminator_index-1)]

wood.dist <- vegdist(wood, method = method)
wood.mds <- metaMDS(wood, try=1000, distance = method, k = 2)

###############################################################################
# Ordination Plot
type <- "none" # none, points, text
groups <- wood.env$Wood.Texture
colors <- viridis(length(unique(groups)))
ordiplot(wood.mds, type=type)
ordiellipse(wood.mds, groups=groups, draw="polygon", col=colors, label=TRUE)
#orditorp(wood.mds, display="sites", cex=0.78, air=0.01)
#orditorp(wood.mds, display="species", cex=0.50, col="red", air=0.01)

env <- data.frame(wood.env$Days.Elapsed)
names(env) <- "Days.Elapsed"
env_vectors <- envfit(wood.mds, env, perm=1000, arrows = TRUE)
plot(env_vectors, cex=0.80, add=TRUE)

###############################################################################
# Analysis of Similarities
wood.ano <- with(wood.env, anosim(wood.dist, Wood.Texture))
wood.ano

###############################################################################
# PERMANOVA

#model <- "wood.dist ~ Days.Elapsed + Weeks.Elapsed + Months.Elapsed + Bimonthly.Elapsed + Season + Wood.Type + Wood.Texture + Days.Elapsed:Weeks.Elapsed + Days.Elapsed:Months.Elapsed + Days.Elapsed:Bimonthly.Elapsed + Days.Elapsed:Season + Days.Elapsed:Wood.Type + Days.Elapsed:Wood.Texture + Weeks.Elapsed:Months.Elapsed + Weeks.Elapsed:Bimonthly.Elapsed + Weeks.Elapsed:Season + Weeks.Elapsed:Wood.Type + Weeks.Elapsed:Wood.Texture + Months.Elapsed:Bimonthly.Elapsed + Months.Elapsed:Season + Months.Elapsed:Wood.Type + Months.Elapsed:Wood.Texture + Bimonthly.Elapsed:Season + Bimonthly.Elapsed:Wood.Type + Bimonthly.Elapsed:Wood.Texture + Season:Wood.Type + Season:Wood.Texture + Wood.Type:Wood.Texture + Days.Elapsed:Weeks.Elapsed:Months.Elapsed + Days.Elapsed:Weeks.Elapsed:Bimonthly.Elapsed + Days.Elapsed:Weeks.Elapsed:Season + Days.Elapsed:Weeks.Elapsed:Wood.Type + Days.Elapsed:Weeks.Elapsed:Wood.Texture + Days.Elapsed:Months.Elapsed:Bimonthly.Elapsed + Days.Elapsed:Months.Elapsed:Season + Days.Elapsed:Months.Elapsed:Wood.Type + Days.Elapsed:Months.Elapsed:Wood.Texture + Days.Elapsed:Bimonthly.Elapsed:Season + Days.Elapsed:Bimonthly.Elapsed:Wood.Type + Days.Elapsed:Bimonthly.Elapsed:Wood.Texture + Days.Elapsed:Season:Wood.Type + Days.Elapsed:Season:Wood.Texture + Days.Elapsed:Wood.Type:Wood.Texture + Weeks.Elapsed:Months.Elapsed:Bimonthly.Elapsed + Weeks.Elapsed:Months.Elapsed:Season + Weeks.Elapsed:Months.Elapsed:Wood.Type + Weeks.Elapsed:Months.Elapsed:Wood.Texture + Weeks.Elapsed:Bimonthly.Elapsed:Season + Weeks.Elapsed:Bimonthly.Elapsed:Wood.Type + Weeks.Elapsed:Bimonthly.Elapsed:Wood.Texture + Weeks.Elapsed:Season:Wood.Type + Weeks.Elapsed:Season:Wood.Texture + Weeks.Elapsed:Wood.Type:Wood.Texture + Months.Elapsed:Bimonthly.Elapsed:Season + Months.Elapsed:Bimonthly.Elapsed:Wood.Type + Months.Elapsed:Bimonthly.Elapsed:Wood.Texture + Months.Elapsed:Season:Wood.Type + Months.Elapsed:Season:Wood.Texture + Months.Elapsed:Wood.Type:Wood.Texture + Bimonthly.Elapsed:Season:Wood.Type + Bimonthly.Elapsed:Season:Wood.Texture + Bimonthly.Elapsed:Wood.Type:Wood.Texture + Season:Wood.Type:Wood.Texture + Days.Elapsed:Weeks.Elapsed:Months.Elapsed:Bimonthly.Elapsed + Days.Elapsed:Weeks.Elapsed:Months.Elapsed:Season + Days.Elapsed:Weeks.Elapsed:Months.Elapsed:Wood.Type + Days.Elapsed:Weeks.Elapsed:Months.Elapsed:Wood.Texture + Days.Elapsed:Weeks.Elapsed:Bimonthly.Elapsed:Season + Days.Elapsed:Weeks.Elapsed:Bimonthly.Elapsed:Wood.Type + Days.Elapsed:Weeks.Elapsed:Bimonthly.Elapsed:Wood.Texture + Days.Elapsed:Weeks.Elapsed:Season:Wood.Type + Days.Elapsed:Weeks.Elapsed:Season:Wood.Texture + Days.Elapsed:Weeks.Elapsed:Wood.Type:Wood.Texture + Days.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Season + Days.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Wood.Type + Days.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Wood.Texture + Days.Elapsed:Months.Elapsed:Season:Wood.Type + Days.Elapsed:Months.Elapsed:Season:Wood.Texture + Days.Elapsed:Months.Elapsed:Wood.Type:Wood.Texture + Days.Elapsed:Bimonthly.Elapsed:Season:Wood.Type + Days.Elapsed:Bimonthly.Elapsed:Season:Wood.Texture + Days.Elapsed:Bimonthly.Elapsed:Wood.Type:Wood.Texture + Days.Elapsed:Season:Wood.Type:Wood.Texture + Weeks.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Season + Weeks.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Wood.Type + Weeks.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Wood.Texture + Weeks.Elapsed:Months.Elapsed:Season:Wood.Type + Weeks.Elapsed:Months.Elapsed:Season:Wood.Texture + Weeks.Elapsed:Months.Elapsed:Wood.Type:Wood.Texture + Weeks.Elapsed:Bimonthly.Elapsed:Season:Wood.Type + Weeks.Elapsed:Bimonthly.Elapsed:Season:Wood.Texture + Weeks.Elapsed:Bimonthly.Elapsed:Wood.Type:Wood.Texture + Weeks.Elapsed:Season:Wood.Type:Wood.Texture + Months.Elapsed:Bimonthly.Elapsed:Season:Wood.Type + Months.Elapsed:Bimonthly.Elapsed:Season:Wood.Texture + Months.Elapsed:Bimonthly.Elapsed:Wood.Type:Wood.Texture + Months.Elapsed:Season:Wood.Type:Wood.Texture + Bimonthly.Elapsed:Season:Wood.Type:Wood.Texture + Days.Elapsed:Weeks.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Season + Days.Elapsed:Weeks.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Wood.Type + Days.Elapsed:Weeks.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Wood.Texture + Days.Elapsed:Weeks.Elapsed:Months.Elapsed:Season:Wood.Type + Days.Elapsed:Weeks.Elapsed:Months.Elapsed:Season:Wood.Texture + Days.Elapsed:Weeks.Elapsed:Months.Elapsed:Wood.Type:Wood.Texture + Days.Elapsed:Weeks.Elapsed:Bimonthly.Elapsed:Season:Wood.Type + Days.Elapsed:Weeks.Elapsed:Bimonthly.Elapsed:Season:Wood.Texture + Days.Elapsed:Weeks.Elapsed:Bimonthly.Elapsed:Wood.Type:Wood.Texture + Days.Elapsed:Weeks.Elapsed:Season:Wood.Type:Wood.Texture + Days.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Season:Wood.Type + Days.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Season:Wood.Texture + Days.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Wood.Type:Wood.Texture + Days.Elapsed:Months.Elapsed:Season:Wood.Type:Wood.Texture + Days.Elapsed:Bimonthly.Elapsed:Season:Wood.Type:Wood.Texture + Weeks.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Season:Wood.Type + Weeks.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Season:Wood.Texture + Weeks.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Wood.Type:Wood.Texture + Weeks.Elapsed:Months.Elapsed:Season:Wood.Type:Wood.Texture + Weeks.Elapsed:Bimonthly.Elapsed:Season:Wood.Type:Wood.Texture" # + Months.Elapsed:Bimonthly.Elapsed:Season:Wood.Type:Wood.Texture + Days.Elapsed:Weeks.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Season:Wood.Type + Days.Elapsed:Weeks.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Season:Wood.Texture + Days.Elapsed:Weeks.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Wood.Type:Wood.Texture + Days.Elapsed:Weeks.Elapsed:Months.Elapsed:Season:Wood.Type:Wood.Texture + Days.Elapsed:Weeks.Elapsed:Bimonthly.Elapsed:Season:Wood.Type:Wood.Texture + Days.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Season:Wood.Type:Wood.Texture + Weeks.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Season:Wood.Type:Wood.Texture + Days.Elapsed:Weeks.Elapsed:Months.Elapsed:Bimonthly.Elapsed:Season:Wood.Type:Wood.Texture"
model <- "wood.dist ~ Days.Elapsed + Season + Wood.Texture + Days.Elapsed:Season + Days.Elapsed:Wood.Texture + Season:Wood.Texture + Days.Elapsed:Season:Wood.Texture"
wood.ado <- adonis2(as.formula(model), data = wood.env)
wood.ado

###############################################################################
# Alpha Diversity Index

index = "shannon"
features_df <- data.frame(Weeks.Elapsed = wood.env$Weeks.Elapsed, 
                          Wood.Texture = wood.env$Wood.Texture)
alpha_indices <- diversity(wood, index = index, 
                    groups = as.formula(features_df$Weeks.Elapsed),
                    MARGIN=1)

alpha_index_df <- data.frame(Feature = as.numeric(names(alpha_indices)), 
                             Alpha_Index = alpha_indices)
alpha_index_df$Feature <- factor(alpha_index_df$Feature, 
                       levels = alpha_index_df$Feature[
                         order(as.numeric(alpha_index_df$Feature))])

ggplot(alpha_index_df, aes(x = Feature , 
                           y = Alpha_Index, 
                           fill = Feature)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = round(Alpha_Index, 2)), 
            position = position_stack(vjust = 1.05),
            color = "black", size = 3.5) +
  theme_minimal() +
  labs(x = "Time", 
       y = index,
       fill = "Weeks.Elapsed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###############################################################################

# List of predictors
predictors <- c("Days.Elapsed", "Weeks.Elapsed", 
                "Months.Elapsed", "Bimonthly.Elapsed",
                "Season", "Wood.Type", "Wood.Texture")

# Generate all non-empty subsets of the predictors list
subsets <- unlist(lapply(1:length(predictors), function(n) {
  combn(predictors, n, simplify = FALSE)
}), recursive = FALSE)

# Flatten the list
subsets <- lapply(subsets, function(x) unlist(x))

# Create interaction terms
interaction_terms <- lapply(subsets, function(set) {
  paste(set, collapse = ":")
})

interaction_terms <- interaction_terms[1:50]

paste(interaction_terms, collapse = " + ")
# Print all formula strings
print(formula_strings)



