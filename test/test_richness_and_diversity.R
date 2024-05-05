library(readr)
library(vegan)
library(viridis)
library(dplyr)

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

group <- 2 # FA23 + SP24
species <- 1 # All species

df <- get_counts_data(group, species)
discriminator_index <- which(names(df) == "_Splitter_")
wood <- df[, (discriminator_index+1):ncol(df)]
wood.env <- df[, 1:(discriminator_index-1)]

###############################################################################
# Species Richness
feature_label <- "Days.Elapsed"

richness <- specnumber(wood)
ggplot(df, aes(x=as.factor(Days.Elapsed), 
               y=richness, 
               fill=Days.Elapsed)) +
  geom_boxplot(fill=viridis(1, alpha=0.70)) +
  xlab(feature_label) +
  ylab("Species Richness") +
  theme_classic()

summary_stats_by_factor <- df %>%
  group_by(as.factor(Days.Elapsed)) %>%
  summarise(
    Mean = mean(richness),
    SD = sd(richness),
    Min = min(richness),
    Max = max(richness),
    .groups = 'drop' 
  )

###############################################################################
# Alpha Diversity


###############################################################################
# Beta Diversity

###############################################################################
# Gamma Diversity

df <- read_csv(data_file[4])

discriminator_index <- which(names(df) == "_Splitter_")
features_df <- df[, 1:(discriminator_index-1)]
transectC_df <- df[, (discriminator_index+1):ncol(df)]
any(rowSums(transectC_df, na.rm = TRUE) == 0)

x <- transectC_df[apply(transectC_df!=0, 1, all), ]

df <- cbind(features_df, transectC_df)

dissimilarity_index <- "bray"
wood <- transectC_df
wood.env <- features_df
wood.dist <- vegdist(transectC_df, method = dissimilarity_index)


wood.mds <- metaMDS(transectC_df, try=1000, distance = dissimilarity_index, k = 2)
wood.ano <- with(wood.env, anosim(wood.dist, Wood.Texture))
wood.ado <- adonis2(wood ~ Season + Wood.Texture, data = wood.env)
#ado <- wood.ado[,c('var', cols)]


#paste0("ANOSIM statistic R: ", wood.ano["statistic"], 
#       "\n     Significance:", wood.ano["signif"])

#wood.mds$stress

library(ggplot2)
data <- t(wood)
species <- rownames(data)
data <- data.frame(species=species, data)
richness <- specnumber(data)
richness <- data.frame(Group=species, Richness=richness)
richness <- richness[order(-richness$Richness), ]
ggplot(richness[1:10, ], aes(x = Group, y = Richness)) +
  geom_barplot() +
  labs(
    title = "Boxplot of Species Richness",
    x = "Group",
    y = "Species Richness"
  ) +
  theme_minimal()


library(ggplot2)

# Example species richness data (replace with your actual data)
species_richness <- data.frame(
  Group = c("Group A", "Group A", "Group B", "Group B", "Group C"),
  Richness = c(10, 15, 12, 18, 14)
)

# Create boxplot
ggplot(species_richness, aes(x = Group, y = Richness)) +
  geom_boxplot() +
  labs(
    title = "Boxplot of Species Richness",
    x = "Group",
    y = "Species Richness"
  ) +
  theme_minimal()