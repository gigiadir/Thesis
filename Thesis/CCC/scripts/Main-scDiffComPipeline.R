# --- Master Runner Config: scDiffCom Analysis ---
# Thesis/CCC/scripts is a symlink to /local/.../Dropbox — often missing on compute nodes. The job
# script and qsub targets live under ~/Scripts on GPFS; copy edits there or sync from Thesis.
script_path <- normalizePath(path.expand("~/Scripts/scDiffComPipeline.R"), mustWork = FALSE)
preprocess_results_path <- "/gpfs0/bgu-ofircohen/users/gigiadir/CCC-PreProcess/results-RankGenes"
base_output_path <- "/gpfs0/bgu-ofircohen/users/gigiadir/CCC-scDiffCom/results/split-by-rank-genes"

# New data location
new_data_dir <- "/gpfs0/bgu-ofircohen/users/gigiadir/scObjects"

logs_path <- file.path(base_output_path, "logs")
genes_path <- file.path("/gpfs0/bgu-ofircohen/users/gigiadir/Thesis/CCC/outputs/data/Complexes.Oncogenes.OncoKB.Cosmic.NCG.rds")
if (!dir.exists(logs_path)) dir.create(logs_path, recursive = TRUE)

script_dir <- dirname(script_path)
wrapper_path <- normalizePath(path.expand("~/Scripts/run_with_scDiffComPipeline_env.sh"), mustWork = FALSE)
#genes <- c("AXL", "ERBB2", "EGFR", "HLA-A", "HLA-B", "HLA-C", "ESR1", "MKI67", "GATA3", "CTLA4", "CDH1", "CDH2", "CDH11", "CTNNB1", "SDC4", "THBS1", "VWF", "COL1A1", "COL1A2", "VEGFA", "IGF1", "IGF2", "CSF1", "CSF1R", "CSF2")
genes <- c("ABI1", "ACTB", "ACTG1", "APH1A", "ARID1A", "ARID1B", "ARID2", "AXL", "BARD1", "BAZ1A", "BCL6", "BLM", "BPTF", "BRCA1", "BRIP1", "BUB1B", "CARM1", "CCNC", "CD28", "CDC27", "CDH1", "CDH11", "CDH2", "CDK8", "CDKN2A", "CHD4", "COL1A1", "COL1A2", "CREBBP", "CSF1", "CSF1R", "CSF2", "CTLA4", "CTNNB1", "CUL3", "CUL4A", "CUL7", "CYFIP1", "DNMT3B", "EGFR", "EP300", "ERBB2", "ERCC2", "ERCC3", "ESR1", "FANCA", "FANCC", "FANCE", "FANCF", "FANCG", "FANCL", "GATA3", "GNA11", "GNAQ", "GNB1", "GPS2", "HDAC1", "HDAC2", "HDAC4", "HDAC7", "HLA-A", "HLA-B", "HLA-C", "IGF1", "IGF2", "LDB1", "LMO2", "LTB", "MAPK1", "MAPK3", "MDM2", "MED1", "MED12", "MKI67", "MLH1", "MNAT1", "NBN", "NCOA3", "NCOR1", "NCOR2", "NCSTN", "NDC80", "NDUFB9", "NPM1", "NUF2", "NUP98", "PARP1", "PBRM1", "PHC2", "PIGA", "POLR2A", "PSMB2", "RABEP1", "RAD21", "RAD50", "RAD51B", "RAD51C", "SDC4", "SIN3A", "SKP2", "SMARCA1", "SMARCA2", "SMARCA4", "SMARCB1", "SMARCD1", "SMARCE1", "SMC1A", "STAG1", "STAG2", "STAT1", "STAT2", "TBL1XR1", "TCEB1", "THBS1", "THRAP3", "TP53", "TRAK1", "VEGFA", "VHL", "VWF", "XRCC1", "XRCC2", "YY1")
# genes <- readRDS(genes_path)


bash_file_path <- file.path(script_dir, "submit_scDiffCom_jobs.sh")
all_commands <- c("#!/bin/bash", "")

# Updated dataset configuration pointing to your scObjects folder
datasets_config <- list(
  # Breast
  # list(path = file.path(new_data_dir, "Bassez2021_Breast_step4.RDS"), column_id = "ident", cell_type_col = "Consensus_Cell_Type"),
  # list(path = file.path(new_data_dir, "Qian2020_Breast_step4.RDS"), column_id = "ident", cell_type_col = "Consensus_Cell_Type"),
  # list(path = file.path(new_data_dir, "Wu2021_Breast_step4.RDS"), column_id = "ident", cell_type_col = "Consensus_Cell_Type")
  # 
  # Lung
  # list(path = file.path(new_data_dir, "Chan2021_Lung_step4.RDS"), column_id = "ident", cell_type_col = "Consensus_Cell_Type"),
  # list(path = file.path(new_data_dir, "Xing2021_Lung_step4.RDS"), column_id = "ident", cell_type_col = "Consensus_Cell_Type"),
  # list(path = file.path(new_data_dir, "Laughney2020_Lung_step4.RDS"), column_id = "ident", cell_type_col = "Consensus_Cell_Type")
  
  # H&N
  # list(path = file.path(new_data_dir, "HNSCC.Atlas.RData"), column_id = "Patient", cell_type_col = "Cell_Type")
  list(path = file.path(new_data_dir, "Kurten_HNSC.RData"), column_id = "Patient", cell_type_col = "Cell_Type"),
  list(path = file.path(new_data_dir, "Puram_HNSC.RData"), column_id = "Patient", cell_type_col = "Cell_Type"),
  list(path = file.path(new_data_dir, "Bill_HNSC.RData"), column_id = "Patient", cell_type_col = "Cell_Type"),
  list(path = file.path(new_data_dir, "Choi_HNSC.RData"), column_id = "Patient", cell_type_col = "Cell_Type")
)

# --- Execution Loop ---
for (current_gene in genes) {
  for (ds in datasets_config) {
    
    current_path <- ds$path
    current_col  <- ds$column_id
    current_ct   <- ds$cell_type_col
    dataset_name <- tools::file_path_sans_ext(basename(current_path))
  
    summary_filename <- sprintf("%s_%s_grouped.rds", current_gene, dataset_name)
    summary_path <- file.path(preprocess_results_path, dataset_name, summary_filename)
    
    output_filename <- sprintf("%s_%s_scDiffCom.rds", current_gene, dataset_name)
    output_path <- file.path(base_output_path, dataset_name, output_filename)
    
    # Check if the preprocess file actually exists
    if (!file.exists(summary_path)) {
      message("Skipping: Preprocess file not found at ", summary_path)
      next
    }
    
    if(file.exists(output_path)) {
      message("Skipping: Output exists on ", output_filename)
      next
      
    }
    
    job_name <- paste0("scDiff_", current_gene, "_", dataset_name)
    
    # qsub: bash runs wrapper (conda env + Rscript); no -V so the job does not inherit the submit shell env
    cmd <- sprintf(
      "qsub -q intel_all.q -S /bin/bash -cwd -N %s -o %s -e %s %s --gene %s --dataset_path %s --patient_summary_path %s --column_id %s --cell_type_col %s --iterations 10000",
      job_name,
      shQuote(logs_path, type = "sh"),
      shQuote(logs_path, type = "sh"),
      shQuote(wrapper_path, type = "sh"),
      shQuote(current_gene, type = "sh"),
      shQuote(current_path, type = "sh"),
      shQuote(summary_path, type = "sh"),
      shQuote(current_col, type = "sh"),
      shQuote(current_ct, type = "sh")
    )
    
    all_commands <- c(all_commands, cmd)
  }
}

writeLines(all_commands, bash_file_path)
message("\nDone! Created submit_scDiffCom_jobs.sh with ", length(all_commands) - 2, " commands.")