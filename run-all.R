# Set the working directory
working_directory <- "C:/Users/Joel/work/kean-stme-2903-11/github.com/joelclim"
setwd(working_directory)

# Set the workflow directory
workflow_directory <- "./mycopins-workflow/"


################################################################################
# workflow
#
source("./mycopins-workflow/workflow.R")

batch_name <- "FA23-transectC"
scata_dataset_name <- "FA23-C-Mycopins"
scata_job_id <- "scata6398"
run_workflow(batch_name, scata_dataset_name, scata_job_id, prompt=FALSE)

batch_name <- "SP24-transectC"
scata_dataset_name <- "SP24-transectC"
scata_job_id <- "scata6471"
run_workflow(batch_name, scata_dataset_name, scata_job_id, prompt=FALSE)

batch_name <- "SP24-transectA-1"
scata_dataset_name <- "SP24-transectA-1"
scata_job_id <- "scata6474"
run_workflow(batch_name, scata_dataset_name, scata_job_id, prompt=FALSE)

batch_name <- "SP24-transectA-2"
scata_dataset_name <- "SP24-transectA-2"
scata_job_id <- "scata6475"
run_workflow(batch_name, scata_dataset_name, scata_job_id, prompt=FALSE)

batch_name <- "SU24-transectA-1"
scata_dataset_name <- "SU24-transectA-1"
scata_job_id <- "scata6540"
run_workflow(batch_name, scata_dataset_name, scata_job_id, prompt=FALSE)

batch_name <- "SU24-transectA-2"
scata_dataset_name <- "SU24-transectA-2"
scata_job_id <- "scata6541"
run_workflow(batch_name, scata_dataset_name, scata_job_id, prompt=FALSE)

################################################################################
# merge
#
source("./mycopins-workflow/merge.R")

configuration <- mycopins_merge_config("TransectC",
                                       "FA23-transectC",
                                       "SP24-transectC")
mycopins_merge(configuration)


configuration <- mycopins_merge_config("SP24-TransectA",
                                       "SP24-transectA-1",
                                       "SP24-transectA-2")
mycopins_merge(configuration)

configuration <- mycopins_merge_config("SU24-TransectA",
                                       "SU24-transectA-1",
                                       "SU24-transectA-2")
mycopins_merge(configuration)

configuration <- mycopins_merge_config("TransectA",
                                       "SP24-transectA",
                                       "SU24-transectA")
mycopins_merge(configuration)

configuration <- mycopins_merge_config("All",
                                       "TransectA",
                                       "TransectC")
mycopins_merge(configuration)


################################################################################
# IPT-GBIF
#
source("./mycopins-workflow/ipt-gbif.R")

configuration <- mycopins_ipt_gbif_config()
mycopins_ipt_gbif_generate(configuration)
