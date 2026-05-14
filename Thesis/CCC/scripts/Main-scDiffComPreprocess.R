# --- Master Runner Config ---
script_path <- "/gpfs0/bgu-ofircohen/users/gigiadir/Thesis/CCC/scripts/scDiffComPreprocess.R"
output_path <- "/gpfs0/bgu-ofircohen/users/gigiadir/CCC-PreProcess"
genes_path <- file.path("/gpfs0/bgu-ofircohen/users/gigiadir/Thesis/CCC/outputs/data/Complexes.Oncogenes.OncoKB.Cosmic.NCG.rds")

results_path <- file.path(output_path, "results")
logs_path <- file.path(output_path, "logs")

script_dir <- dirname(script_path)

#gene_rds_path <- file.path(script_dir, "final_gene_selection.rds")

#genes <- c("AXL", "ERBB2", "EGFR", "HLA-A", "HLA-B", "HLA-C", "ESR1", "MKI67", "GATA3", "CTLA4", "CDH1", "CDH2", "CDH11", "CTNNB1", "SDC4", "THBS1", "VWF", "COL1A1", "COL1A2", "VEGFA", "IGF1", "IGF2", "CSF1", "CSF1R", "CSF2")
# genes <- c()
 # genes <- readRDS(genes_path)
 # genes <- genes[50:100]
genes <- c(
  "ABI1", "ACTG1", "APH1A", "ARID1B", "ARID2", "BARD1", "BAZ1A", "BLM", 
  "BRCA1", "BRIP1", "BUB1B", "CCNC", "CDK8", "CDKN2A", "CUL3", "CUL4A", 
  "CYFIP1", "FANCE", "FANCL", "GNAQ", "GNB1", "LDB1", "MLH1", "NBN", 
  "NDC80", "NDUFB9", "NPM1", "PARP1", "PBRM1", "PSMB2", "RAD21", "RAD50", 
  "RAD51B", "RAD51C", "SMARCA2", "SMC1A", "STAG1", "TCEB1", "VHL", 
  "XRCC1", "XRCC2"
)

bash_file_path <- file.path(script_dir, "submit_jobs.sh")

all_commands <- c("#!/bin/bash", "") 

# Updated destination path
new_data_dir <- "/gpfs0/bgu-ofircohen/users/gigiadir/scObjects"

# Define datasets with updated paths
datasets_config <- list(
  # Breast
  # list(
  #   path = file.path(new_data_dir, "Bassez2021_Breast_step4.RDS"),
  #   column_id = "ident",
  #   cell_type_col = "Consensus_Cell_Type"
  # ),
  # list(
  #   path = file.path(new_data_dir, "Qian2020_Breast_step4.RDS"),
  #   column_id = "ident",
  #   cell_type_col = "Consensus_Cell_Type"
  # ),
  # list(
  #   path = file.path(new_data_dir, "Wu2021_Breast_step4.RDS"),
  #   column_id = "ident",
  #   cell_type_col = "Consensus_Cell_Type"
  # )
  
  # Lung
  # list(
  #   path = file.path(new_data_dir, "Chan2021_Lung_step4.RDS"),
  #   column_id = "ident",
  #   cell_type_col = "Consensus_Cell_Type"
  # ),
  # list(
  #   path = file.path(new_data_dir, "Xing2021_Lung_step4.RDS"),
  #   column_id = "ident",
  #   cell_type_col = "Consensus_Cell_Type"
  # ),
  # list(
  #   path = file.path(new_data_dir, "Laughney2020_Lung_step4.RDS"),
  #   column_id = "ident",
  #   cell_type_col = "Consensus_Cell_Type"
  # )
  # list(
  #   path = file.path(new_data_dir, "Bischoff2021_Lung_step4.RDS"),
  #   column_id = "ident",
  #   cell_type_col = "Consensus_Cell_Type"
  # ),
  
  # H&N
  # list(
  #   path = file.path(new_data_dir, "HNSCC.Atlas.RData"),
  #   column_id = "Patient",
  #   cell_type_col = "Cell_Type"
  # ),
  list(
    path = file.path(new_data_dir, "Kurten_HNSC.RData"),
    column_id = "Patient",
    cell_type_col = "Cell_Type"
  ),
  list(
    path = file.path(new_data_dir, "Puram_HNSC.RData"),
    column_id = "Patient",
    cell_type_col = "Cell_Type"
  ),
  list(
    path = file.path(new_data_dir, "Bill_HNSC.RData"),
    column_id = "Patient",
    cell_type_col = "Cell_Type"
  ),
  list(
    path = file.path(new_data_dir, "Choi_HNSC.RData"),
    column_id = "Patient",
    cell_type_col = "Cell_Type"
  )
)

# --- Execution Loop ---
for (current_gene in genes) {
  for (ds in datasets_config) {
    
    current_path <- ds$path
    current_col  <- ds$column_id
    current_ct   <- ds$cell_type_col
    dataset_name <- tools::file_path_sans_ext(basename(current_path))
    
    summary_filename <- sprintf("%s_%s_grouped.rds", current_gene, dataset_name)
    summary_path <- file.path(results_path, dataset_name, summary_filename)
    if (file.exists(summary_path)) {
      next
    }
    message("\n" , paste(rep("=", 50), collapse = ""))
    message("RUNNING: Gene =", current_gene, " | Dataset =", basename(current_path))
    
    # 1. Start with the basic command
    job_name <- paste0(current_gene, "_", dataset_name)
    
    cmd <- sprintf(
      "qsub -q intel_all.q -S /bin/Rscript -V -cwd -N %s -o %s -e %s %s --gene %s --dataset_path %s --output_path %s",
      job_name,
      logs_path,
      logs_path,
      script_path,
      current_gene,
      current_path,
      results_path
    )
    
    if (!is.null(current_col)) cmd <- paste(cmd, "--column_id", current_col)
    if (!is.null(current_ct))  cmd <- paste(cmd, "--cell_type_col", current_ct)
    
    all_commands <- c(all_commands, cmd)
  }
}

writeLines(all_commands, bash_file_path)