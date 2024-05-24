library(readr)
library(vegan)
library(viridis)
library(dplyr)

working_folder_path <- "C:/Users/Joel/work/kean-stme-2903-11/github.com/joelclim"
setwd(working_folder_path)


get_community_data <- function(transect, species) {
  if (1 == species) {
    community_file <- paste0("./data/", transect, "/mycopins_community.csv")
  } else if (2 == species) { # genus/species
    community_file <- paste0("./data/", transect, "/mycopins_community_fungi.csv")
  } else if (3 == species) { # wood saprotrophs fungi
    community_file <- paste0("./data/", transect, "/mycopins_community_wood_saprotroph.csv")
  }
  
  community_df <- read_csv(community_file, show_col_types = FALSE)
  
  return(community_df)
}


get_environment_data <- function(transect, species) {
  if (1 == species) {
    environment_file <- paste0("./data/", transect, "/mycopins_environment.csv")
  } else if (2 == species) { # genus/species
    environment_file <- paste0("./data/", transect, "/mycopins_environment_fungi.csv")
  } else if (3 == species) { # wood saprotrophs fungi
    environment_file <- paste0("./data/", transect, "/mycopins_environment_wood_saprotroph.csv")
  }
  
  environment_df <- read_csv(environment_file, show_col_types = FALSE)
  
  return(environment_df)
}


method <- "bray"
batch <- "All"
species <- 1 # All species

mycopins <- get_community_data(batch, species)
mycopins.env <- get_environment_data(batch, species)

mycopins.dist <- vegdist(mycopins, method = method)
mycopins.mds <- metaMDS(mycopins, try=1000, distance = method, k = 2)

###############################################################################
# Species Frequency
# mycopins.data <- cbind(mycopins.env, mycopins)
# x_label <- "Days.Elapsed"
# species_frequency <- t(mycopins.data) %>%
#   group_by_at(x_label) %>%
#   summarize(
#     Count = n()
#   )
# species_frequency



###############################################################################
# Species Richness
x_label <- "Days.Elapsed"
y_label <- "SpeciesRichness"

richness <- specnumber(mycopins)

data <- mycopins
data.env <- mycopins.env
data.env$IndependentVariable <- data.env[[x_label]]
data.env$SpeciesRichness <- log2(specnumber(data))
ggplot(data.env, aes(x = as.factor(IndependentVariable), y = SpeciesRichness)) +
  geom_boxplot() +
  labs(x = x_label, y = y_label) +
  theme_classic()

box_plot_data <- data.env[, c(x_label, y_label)]

summary_data <- box_plot_data %>%
  group_by_at(x_label) %>%
  summarize(
    Count = n(),
    Mean = mean(SpeciesRichness),
    Median = median(SpeciesRichness),
    Min = min(SpeciesRichness),
    Max = max(SpeciesRichness),
    SD = sd(SpeciesRichness)
  )
summary_data


###############################################################################
# Alpha Diversity


###############################################################################
# Beta Diversity
beta_diversity_df <- as.data.frame(as.matrix(mycopins.dist))

###############################################################################
# Gamma Diversity

df <- read_csv(data_file[4])

discriminator_index <- which(names(df) == "_Splitter_")
env_df <- df[, 1:(discriminator_index-1)]
transectC_df <- df[, (discriminator_index+1):ncol(df)]
any(rowSums(transectC_df, na.rm = TRUE) == 0)

x <- transectC_df[apply(transectC_df!=0, 1, all), ]

df <- cbind(env_df, transectC_df)

dissimilarity_index <- "bray"
mycopins.dist <- vegdist(mycopins, method = dissimilarity_index)
mycopins.mds <- metaMDS(mycopins, try=1000, distance = dissimilarity_index, k = 2)
mycopins.ano <- with(mycopins.env, anosim(wood.dist, Wood.Texture))
mycopins.ado <- adonis2(mycopins ~ Season + Wood.Texture, data = mycopins.env)
#ado <- wood.ado[,c('var', cols)]


#paste0("ANOSIM statistic R: ", wood.ano["statistic"],
#       "\n     Significance:", wood.ano["signif"])

#wood.mds$stress

library(ggplot2)
data <- t(mycopins)
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