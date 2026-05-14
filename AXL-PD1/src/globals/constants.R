##PROJECT CONSTANTS
get_working_dir_path <- function() {
  nodename <- Sys.info()["nodename"]
  if(nodename == "MacBook-Pro.local") {
    return("~/Documents/BGU/AXL-PD1")
  } else if (nodename == "bhn1095") {
    return("~/AXL-PD1")
  } else {
    stop("BASE_DIR_PATH is not defined for this environment")
  }
} 

BASE_DIR_PATH <- get_working_dir_path()
PROJECT_NAME <- "AXL-PD1"
CODE_DIR_PATH <- file.path(BASE_DIR_PATH, "src")
DATA_DIR_PATH <- file.path(BASE_DIR_PATH, "data")
OUTPUT_DIR_PATH <- file.path(BASE_DIR_PATH, "output")
DROPBOX_DIR_PATH <- file.path("/local/ofir/OCohenLab/Dropbox")
CANCER_DB_DIR_PATH <- file.path(DROPBOX_DIR_PATH, "Cancer_DB")