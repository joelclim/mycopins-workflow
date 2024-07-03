working_directory <- "C:/Users/Joel/work/kean-stme-2903-11/github.com/joelclim"
data_directory <- paste0(working_directory, "/data/All")
community_file <- paste0(data_directory, "/mycopins_community_fungi.csv")
environment_file <- paste0(data_directory, "/mycopins_environment_fungi.csv")
gbif_otu_table_file <- paste0(data_directory, "/OTU_table.tsv")

library(dplyr)
library(jsonlite)
library(readr)

mycopins <- read_csv(community_file, show_col_types = FALSE,
                     locale = locale(encoding = "UTF-8"))
mycopins.env <- read_csv(environment_file, show_col_types = FALSE,
                         locale = locale(encoding = "UTF-8"))

otu_table <- cbind(mycopins.env$Sample.Number, mycopins)
otu_table <- t(otu_table)
otu_table <- as.data.frame(otu_table)
colnames(otu_table) <- otu_table[1,]
otu_table <- otu_table[-1,]
otu_table <- otu_table %>%
  mutate_all(as.numeric)
otu_table <- cbind(id=rownames(otu_table), otu_table)

write.table(otu_table, gbif_otu_table_file, sep = "\t", 
            row.names = FALSE, na = "", fileEncoding = "UTF-8")


# Issues:
# We don't use OTU ids. OTUs in scata can be combined if they are determined to represent the same species