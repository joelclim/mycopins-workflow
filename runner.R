working_directory <- "C:/Users/Joel/work/kean-stme-2903-11/github.com/joelclim/mycopins"
setwd(working_directory)

source("./workflow/mycopins-workflow.R")

###############################################################################
# Batch to run
#
#batch_name <- "FA23-transectC"
#scata_dataset_name <- "FA23-C-Mycopins"
#scata_job_id <- "scata6398"

#batch_name <- "SP24-transectC-1"
#scata_dataset_name <- "SP24-RFI-TranC \\(Full\\)"
#scata_job_id <- "scata6397"

batch_name <- "SP24-transectC-2"
scata_dataset_name <- "SP24-transectC-2"
scata_job_id <- "scata6466"

###############################################################################
# Run workflow
#
configuration <- mycopins_config(batch_name, scata_dataset_name, scata_job_id)

mycopins_preprocess(configuration)
mycopins_search(configuration)
mycopins_identify(configuration)
mycopins_annotate(configuration)


###############################################################################
# Auxiliary functions
#
mycopins_gbif_taxonomy(configuration)
mycopins_fungal_traits(configuration)
