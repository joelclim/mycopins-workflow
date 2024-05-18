library(readr)
library(vegan)
library(viridis)
library(dplyr)

working_folder_path <- "C:/Users/Joel/work/kean-stme-2903-11/github.com/joelclim"
setwd(working_folder_path)

get_counts_data <- function(batch, species) {
  if (1 == species) {
    counts_file <- paste0("./data/", batch, "/complete_counts.csv")
  } else if (2 == species) {        # genus/species
    counts_file <- paste0("./data/", batch, "/complete_counts_fungi_only.csv")
  } else if (3 == species) { # wood saprotrophs fungi
    counts_file <- paste0("./data/", batch, "/complete_counts_wood_saprotrophs_only.csv")
  }
  
  counts_df <- read_csv(counts_file)
  
  return(counts_df)
}

method <- "bray"
batch <- "SP24-transectC-1"
species <- 1 # All species

df <- get_counts_data(batch, species)
discriminator_index <- which(names(df) == "_Splitter_")
wood <- df[, (discriminator_index+1):ncol(df)]
wood.env <- df[, 1:(discriminator_index-1)]

wood.dist <- vegdist(wood, method = method)
wood.mds <- metaMDS(wood, try=1000, distance = method, k = 2)

###############################################################################
# Species Richness
x_label <- "Days.Elapsed"
y_label <- "SpeciesRichness"

richness <- specnumber(wood)

data <- wood
data.env <- wood.env
data.env$IndependentVariable <- data.env[[x_label]]
data.env$SpeciesRichness <- specnumber(data)

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
beta_diversity_df <- as.data.frame(as.matrix(wood.dist))

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
wood <- transectC_df
wood.env <- env_df
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