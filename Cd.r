# ============================================================================
# Single-Cell RNA-seq Analysis Pipeline for Testis Data
# ============================================================================

# 1. INITIAL SETUP -----------------------------------------------------------
# Create and initialize project environment
renv::init(project = "~/SCP_env", bare = TRUE, restart = TRUE)

# Install required packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SCP", "Seurat", "DoubletFinder", "devtools", "dplyr"))

# Load libraries
library(SCP)
library(Seurat)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(tibble)

# 2. DATA LOADING ------------------------------------------------------------
#' Load 10X data from directories matching pattern
load_10x_data <- function(pattern = "STE") {
  fileNames <- list.files(pattern = pattern)
  message("Found ", length(fileNames), " samples: ", paste(fileNames, collapse = ", "))
  
  # Load each sample into a named list
  samples <- list()
  for (sample in fileNames) {
    message("Loading sample: ", sample)
    samples[[sample]] <- Read10X(data.dir = sample)
  }
  return(samples)
}

# Load all samples
samples <- load_10x_data("STE")
list2env(samples, envir = .GlobalEnv)

# 3. SAMPLE NAMES ------------------------------------------------------------
sample_names <- c("C72", "C74", "C82", "CDH233", "CDH234", "CDH245", 
                  "CTL139", "CTL209", "Ctrl1", "Ctrl2", "Ctrl3",
                  "DEHP107", "DEHP109", "DEHP122", "FXT113", "FXT137", 
                  "FXT144", "MTX55", "MTX63", "MTX78", "STEKO135", 
                  "STEKO157", "STEKO182", "STEWT140", "STEWT143", 
                  "STEWT145", "WSCTL12", "WSCTL13", "WSCTL15", 
                  "WSCTL23", "WST12", "WST13", "WST15", "WST23")

# 4. CREATE SEURAT OBJECTS ---------------------------------------------------
#' Create Seurat objects for all samples
create_seurat_objects <- function(sample_names) {
  seurat_objects <- list()
  
  for (sample in sample_names) {
    if (exists(sample)) {
      message("Creating Seurat object for: ", sample)
      seurat_objects[[sample]] <- CreateSeuratObject(
        counts = get(sample),
        project = sample,
        min.cells = 3,
        min.features = 300
      )
    }
  }
  
  # Assign to global environment
  list2env(seurat_objects, envir = .GlobalEnv)
  return(seurat_objects)
}

seurat_objects <- create_seurat_objects(sample_names)

# 5. QC AND FILTERING --------------------------------------------------------
#' Add mitochondrial percentage and filter cells
process_seurat_object <- function(obj, mt_pattern = "^mt-") {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  obj <- subset(obj, subset = nFeature_RNA > 300 & nFeature_RNA < 10000 & percent.mt < 20)
  return(obj)
}

# Process all objects
for (name in names(seurat_objects)) {
  message("Processing: ", name)
  seurat_objects[[name]] <- process_seurat_object(seurat_objects[[name]])
}

# Update global environment
list2env(seurat_objects, envir = .GlobalEnv)

# 6. NORMALIZATION AND PCA ---------------------------------------------------
#' Normalize, find variable features, scale, and run PCA
run_pca_pipeline <- function(obj) {
  obj <- obj %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)
  return(obj)
}

# Apply to all objects
for (name in names(seurat_objects)) {
  message("Running PCA pipeline for: ", name)
  seurat_objects[[name]] <- run_pca_pipeline(seurat_objects[[name]])
}

list2env(seurat_objects, envir = .GlobalEnv)

# 7. DOUBLET DETECTION -------------------------------------------------------
#' Run DoubletFinder on a Seurat object
run_doubletfinder <- function(obj, PCs = 1:30) {
  # Parameter sweep
  sweep.res <- paramSweep_v3(obj, PCs = PCs, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # Get optimal pK
  pK_value <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  # Calculate doublet rate (based on 8 doublets per 1000 cells)
  DoubletRate <- ncol(obj) * 8 * 1e-6
  nExp_poi <- round(DoubletRate * length(obj@meta.data$orig.ident))
  
  # Run DoubletFinder
  obj <- doubletFinder_v3(obj, PCs = PCs, pN = 0.25, 
                          pK = pK_value, nExp = nExp_poi, 
                          reuse.pANN = FALSE, sct = FALSE)
  
  return(obj)
}

# Apply DoubletFinder to select objects
doublet_objects <- c("C72", "C74", "C82", "CDH233", "CDH234", "CDH245", 
                     "CTL139", "CTL209", "Ctrl1", "Ctrl2", "Ctrl3",
                     "DEHP107", "DEHP109", "DEHP122", "FXT113", "FXT137", 
                     "FXT144", "MTX55", "MTX63", "MTX78", "STEKO135", 
                     "STEKO157", "STEKO182", "STEWT140", "STEWT143", 
                     "STEWT145", "WSCTL12", "WSCTL13", "WSCTL15", 
                     "WSCTL23", "WST12", "WST13", "WST15", "WST23")

for (obj_name in doublet_objects) {
  if (exists(obj_name)) {
    message("Running DoubletFinder for: ", obj_name)
    assign(paste0(obj_name, ".DF"), run_doubletfinder(get(obj_name)))
  }
}

# 8. EXTRACT DOUBLET CLASSIFICATIONS -----------------------------------------
#' Extract doublet classifications
extract_doublet_class <- function(df_obj, pN = 0.25, pK = NULL, nExp = NULL) {
  # Find the classification column
  class_col <- grep("DF.classifications", colnames(df_obj@meta.data), value = TRUE)[1]
  if (length(class_col) > 0) {
    df_obj@meta.data$Doublet <- df_obj@meta.data[[class_col]]
  }
  return(df_obj)
}

# Apply to all DF objects
df_objects <- ls(pattern = "\\.DF$")
for (obj_name in df_objects) {
  assign(obj_name, extract_doublet_class(get(obj_name)))
}

# 9. MERGE ALL SAMPLES -------------------------------------------------------
#' Get all DF objects for merging
get_df_objects <- function() {
  df_names <- ls(pattern = "\\.DF$")
  df_list <- list()
  for (name in df_names) {
    df_list[[name]] <- get(name)
  }
  return(df_list)
}

df_list <- get_df_objects()

# Create cell IDs for merging
cell_ids <- gsub("\\.DF", "", names(df_list))

# Merge all objects
testis_total <- merge(
  df_list[[1]],
  y = df_list[-1],
  add.cell.ids = cell_ids,
  project = "testis"
)

# 10. CLEAN METADATA ---------------------------------------------------------
testis_total@meta.data <- testis_total@meta.data %>%
  dplyr::select(orig.ident, nCount_RNA, nFeature_RNA, percent.mt, Doublet)

# Remove doublets
testis_total.DF <- subset(testis_total, subset = Doublet == "Singlet")

# 11. INTEGRATION WITH HARMONY ------------------------------------------------
#' Run Harmony integration
run_harmony_integration <- function(obj, group_var = "orig.ident", vars_to_regress = "percent.mt") {
  obj <- obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
    ScaleData() %>%
    SCTransform(vars.to.regress = vars_to_regress)
  
  obj <- RunPCA(obj, assay = "SCT", npcs = 50)
  
  obj <- RunHarmony(obj, 
                    group.by.vars = group_var, 
                    reduction = "pca", 
                    assay.use = "SCT", 
                    reduction.save = "harmony")
  
  obj <- obj %>%
    RunUMAP(reduction = "harmony", dims = 1:50) %>%
    FindNeighbors(reduction = "harmony", dims = 1:50) %>%
    FindClusters(resolution = 0.6)
  
  return(obj)
}

testis_total.DF <- run_harmony_integration(testis_total.DF)

# 12. SAVE RESULTS -----------------------------------------------------------
saveRDS(testis_total, file = "testis_total.rds")
saveRDS(testis_total.DF, file = "testis_total.DF.rds")

# ============================================================================
# ALCOHOL EXPOSURE ANALYSIS
# ============================================================================

# 13. SUBSET ALCOHOL SAMPLES -------------------------------------------------
alcohol_samples <- c("WSCTL12", "WSCTL13", "WSCTL15", "WSCTL23", 
                     "CTL139", "CTL209", "WST12", "WST13", "WST15", "WST23")

alcohol <- subset(testis_total, idents = alcohol_samples)
alcohol.Singlet <- subset(alcohol, subset = Doublet == "Singlet")

# 14. INTEGRATION WITH SCP ---------------------------------------------------
alcohol.Singlet.harmony <- Integration_SCP(
  srtMerge = alcohol.Singlet, 
  batch = "orig.ident", 
  nHVF = 4000,
  integration_method = "Harmony"
)

alcohol.Singlet.harmony <- alcohol.Singlet.harmony %>%
  RunUMAP(dims = 1:30, reduction = "Harmony") %>%
  FindNeighbors(dims = 1:30, reduction = "Harmony") %>%
  FindClusters(resolution = 3)

# 15. CELL TYPE ANNOTATION FUNCTION ------------------------------------------
#' Annotate cell types based on cluster markers
annotate_cell_types <- function(obj, resolution_col = "RNA_snn_res.2") {
  obj@meta.data <- obj@meta.data %>%
    mutate(celltype = case_when(
      !!sym(resolution_col) == "50" ~ "Innate Lymph",
      !!sym(resolution_col) == "46" ~ "Macrophage",
      !!sym(resolution_col) == "39" ~ "Telocytes",
      !!sym(resolution_col) == "43" ~ "Leydig",
      !!sym(resolution_col) %in% c("0", "20", "45", "53") ~ "Sertoli",
      !!sym(resolution_col) %in% c("16", "23", "26", "27", "49", "18") ~ "Round STids",
      !!sym(resolution_col) %in% c("33", "35", "38", "40") ~ "SPG",
      !!sym(resolution_col) %in% c("44", "15", "21", "51", "7", "13", "22", "10", 
                                   "14", "1", "42", "36", "19", "4", "5", "6", "9") ~ "Elongating STids",
      !!sym(resolution_col) %in% c("28", "30", "31", "37", "12", "3", "2", "25", 
                                   "11", "34", "24", "17", "32", "47", "48", "52") ~ "Scytes",
      TRUE ~ as.character(!!sym(resolution_col))
    ))
  return(obj)
}

alcohol.Singlet.harmony <- annotate_cell_types(alcohol.Singlet.harmony)

# 16. ADD GROUP INFORMATION --------------------------------------------------
alcohol.Singlet.harmony@meta.data <- alcohol.Singlet.harmony@meta.data %>%
  mutate(group = case_when(
    orig.ident %in% c("CTL139", "CTL209", "WSCTL12", "WSCTL13", "WSCTL15", "WSCTL23") ~ "Normal",
    orig.ident %in% c("WST12", "WST13", "WST15", "WST23") ~ "Alcohol",
    TRUE ~ NA_character_
  ))

Idents(alcohol.Singlet.harmony) <- "celltype"

# 17. FIND ALL MARKERS -------------------------------------------------------
alcohol.Singlet.harmony <- RunDEtest(alcohol.Singlet.harmony, group_by = "celltype")

# 18. FILTER CELL TYPES ------------------------------------------------------
keep_celltypes <- c("Scytes", "Sertoli", "Elongating STids", "Round STids", 
                    "SPG", "Innate Lymph", "Leydig", "Telocytes", "Macrophage")

alcohol.Singlet.harmony.filter <- subset(alcohol.Singlet.harmony, 
                                         idents = keep_celltypes)

# 19. SUBSET SPECIFIC CELL TYPES ---------------------------------------------
alcohol.macrophage <- subset(alcohol.Singlet.harmony.filter, idents = "Macrophage")
alcohol.sertoli <- subset(alcohol.Singlet.harmony.filter, idents = "Sertoli")
alcohol.leydig <- subset(alcohol.Singlet.harmony.filter, idents = "Leydig")

# 20. MACROPHAGE SUBTYPING ---------------------------------------------------
process_cell_subtype <- function(obj, resolution = 3) {
  obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:45, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:45)
  obj <- FindClusters(obj, resolution = resolution)
  return(obj)
}

alcohol.macrophage <- process_cell_subtype(alcohol.macrophage)

# Annotate macrophage subtypes
alcohol.macrophage@meta.data <- alcohol.macrophage@meta.data %>%
  mutate(celltype = case_when(
    SCT_snn_res.1.8 %in% c("0", "4", "5", "6", "2") ~ "Cd74",
    SCT_snn_res.3 == "2" ~ "Rcan1",
    SCT_snn_res.3 %in% c("5", "8", "12") ~ "Lcn2",
    TRUE ~ as.character(SCT_snn_res.1.8)
  )) %>%
  mutate(celltype = case_when(
    celltype == "3" ~ "Cd74",
    celltype == "1" ~ "Rcan1",
    TRUE ~ as.character(celltype)
  ))

# 21. SERTOLI CELL SUBTYPING -------------------------------------------------
alcohol.sertoli <- process_cell_subtype(alcohol.sertoli, resolution = 0.4)

alcohol.sertoli@meta.data <- alcohol.sertoli@meta.data %>%
  mutate(celltype = case_when(
    SCT_snn_res.0.4 %in% c("3", "5") ~ "Sdc4",
    SCT_snn_res.0.4 == "2" ~ "Cypt4",
    SCT_snn_res.0.4 == "4" ~ "Elongating_sertoli",
    SCT_snn_res.0.4 %in% c("1", "0") ~ "Chchd10",
    SCT_snn_res.0.4 == "6" ~ "sertoli_unknown",
    TRUE ~ as.character(SCT_snn_res.0.4)
  ))

# 22. SPG SUBTYPING ----------------------------------------------------------
alcohol.SPG <- subset(alcohol.Singlet.harmony.filter, idents = "SPG")
alcohol.SPG <- process_cell_subtype(alcohol.SPG, resolution = 0.4)

alcohol.SPG@meta.data <- alcohol.SPG@meta.data %>%
  mutate(celltype = case_when(
    SCT_snn_res.0.4 %in% c("9", "0", "1", "5", "3", "8") ~ "Undiff",
    SCT_snn_res.0.4 == "6" ~ "Diff",
    TRUE ~ "Diff"
  ))

# 23. MERGE ANNOTATIONS ------------------------------------------------------
#' Merge annotations from subtype objects
merge_annotations <- function(main_obj, subtype_obj, id_var = "celltype") {
  main_obj@meta.data <- main_obj@meta.data %>%
    mutate(CellID = rownames(main_obj@meta.data))
  
  subtype_meta <- subtype_obj@meta.data %>%
    mutate(CellID = rownames(subtype_obj@meta.data)) %>%
    dplyr::select(celltype, CellID)
  
  main_obj@meta.data <- left_join(main_obj@meta.data, subtype_meta, by = "CellID")
  rownames(main_obj@meta.data) <- main_obj@meta.data$CellID
  
  return(main_obj)
}

alcohol.Singlet.harmony.filter <- merge_annotations(alcohol.Singlet.harmony.filter, alcohol.SPG)
alcohol.Singlet.harmony.filter <- merge_annotations(alcohol.Singlet.harmony.filter, alcohol.sertoli)
alcohol.Singlet.harmony.filter <- merge_annotations(alcohol.Singlet.harmony.filter, alcohol.macrophage)

# 24. FINALIZE ANNOTATIONS ---------------------------------------------------
alcohol.Singlet.harmony.filter@meta.data <- alcohol.Singlet.harmony.filter@meta.data %>%
  mutate(celltype = coalesce(celltype.y.y, celltype.x.x)) %>%
  mutate(celltype = case_when(
    celltype == "Undiff" ~ "Undifferentiated SSCs",
    celltype == "Diff" ~ "Differentiated SSCs",
    celltype == "Chchd10" ~ "Chchd10+ Sertoli",
    celltype == "Cypt4" ~ "Cypt4+ Sertoli",
    celltype == "Sdc4" ~ "Sdc4+ Sertoli",
    celltype == "Cd74" ~ "Cd74+ Macrophages",
    celltype == "Rcan1" ~ "Rcan1+ Macrophages",
    celltype == "Lcn2" ~ "Lcn2+ Macrophages",
    TRUE ~ as.character(celltype)
  )) %>%
  mutate(CellType = case_when(
    grepl("SSCs|SPG", celltype) ~ "SPG",
    grepl("Sertoli", celltype) ~ "Sertoli",
    grepl("Macrophages", celltype) ~ "Macrophages",
    TRUE ~ as.character(celltype)
  )) %>%
  dplyr::select(orig.ident, percent.mt, Doublet, nCount_RNA, nFeature_RNA, 
                celltype, group, CellType)

# 25. SUBSET FINAL CELL TYPES ------------------------------------------------
final_celltypes <- c("Scytes", "Cypt4+ Sertoli", "Elongating STids", "Round STids",
                     "Undifferentiated SSCs", "Innate Lymph", "Chchd10+ Sertoli",
                     "Leydig", "Sdc4+ Sertoli", "Telocytes", "Cd74+ Macrophages",
                     "Rcan1+ Macrophages", "Differentiated SSCs", "Lcn2+ Macrophages")

alcohol.final <- subset(alcohol.Singlet.harmony.filter, idents = final_celltypes)

# 26. SAVE RESULTS -----------------------------------------------------------
saveRDS(alcohol.Singlet.harmony.filter, file = "alcohol.Singlet.harmony.filter.rds")
saveRDS(alcohol.final, file = "alcohol.final.rds")
saveRDS(alcohol.SPG, file = "alcohol.SPG.rds")
saveRDS(alcohol.sertoli, file = "alcohol.sertoli.rds")
saveRDS(alcohol.macrophage, file = "alcohol.macrophage.rds")

# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

# 27. PLOTTING FUNCTIONS ----------------------------------------------------
plot_cell_types <- function(obj, group_var = "CellType", split_var = NULL) {
  p <- CellDimPlot(
    srt = obj, 
    group.by = group_var,
    reduction = "umap", 
    theme_use = "theme_blank",
    label = TRUE
  )
  
  if (!is.null(split_var)) {
    p <- p + facet_wrap(as.formula(paste("~", split_var)))
  }
  
  return(p)
}

plot_cell_proportions <- function(obj, stat_var = "group", group_var = "CellType") {
  CellStatPlot(
    obj, 
    stat.by = stat_var, 
    group.by = group_var, 
    stat_type = "percent", 
    position = "dodge", 
    label = TRUE
  )
}

# Generate plots
umap_plot <- plot_cell_types(alcohol.final)
proportion_plot <- plot_cell_proportions(alcohol.final)

# 28. DIFFERENTIAL EXPRESSION ANALYSIS ---------------------------------------
#' Run DE analysis for all cell types
run_de_analysis <- function(obj, group_var = "CellType") {
  obj <- RunDEtest(obj, group_by = group_var)
  
  # Filter significant markers
  de_filter <- dplyr::filter(
    obj@tools[[paste0("DEtest_", group_var)]]$AllMarkers_wilcox,
    p_val_adj < 0.05 & avg_log2FC > 1
  )
  
  # Get top markers
  de_top <- de_filter %>%
    group_by(gene) %>%
    top_n(1, avg_log2FC) %>%
    group_by(group1) %>%
    top_n(3, avg_log2FC)
  
  return(list(obj = obj, de_filter = de_filter, de_top = de_top))
}

de_results <- run_de_analysis(alcohol.final)
alcohol.final <- de_results$obj

# 29. HEATMAP VISUALIZATION -------------------------------------------------
GroupHeatmap(
  alcohol.final,
  features = de_results$de_top$gene,
  feature_split = de_results$de_top$group1,
  group.by = "CellType",
  heatmap_palette = "YlOrRd",
  add_dot = TRUE,
  add_bg = TRUE,
  show_row_names = TRUE
)

# ============================================================================
# ZMPSTE24 KO ANALYSIS
# ============================================================================

# 30. SUBSET ZMPSTE SAMPLES -------------------------------------------------
zmpste_samples <- c("STEKO135", "STEKO157", "STEKO182", 
                    "STEWT140", "STEWT143", "STEWT145")

ZMPSTE <- subset(testis_total, idents = zmpste_samples)
ZMPSTE.Singlet <- subset(ZMPSTE, subset = Doublet == "Singlet")

# 31. PROCESS ZMPSTE --------------------------------------------------------
ZMPSTE.Singlet <- NormalizeData(ZMPSTE.Singlet) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA()

ZMPSTE.Singlet <- FindNeighbors(ZMPSTE.Singlet, dims = 1:30) %>%
  FindClusters(resolution = 0.8) %>%
  RunUMAP(dims = 1:30)

# 32. ANNOTATE ZMPSTE -------------------------------------------------------
ZMPSTE.Singlet <- annotate_cell_types(ZMPSTE.Singlet, resolution_col = "RNA_snn_res.2")

ZMPSTE.Singlet@meta.data <- ZMPSTE.Singlet@meta.data %>%
  mutate(group = case_when(
    orig.ident %in% c("STEKO135", "STEKO157", "STEKO182") ~ "ZMPSTE KO",
    orig.ident %in% c("STEWT140", "STEWT143", "STEWT145") ~ "Normal",
    TRUE ~ NA_character_
  ))

# 33. FILTER AND SUBSET ZMPSTE ----------------------------------------------
zmpste_keep <- c("Scytes", "Sertoli", "Elongating STids", "Round STids", 
                 "SPG", "Innate Lymph", "Leydig", "Telocytes", "Macrophage")

ZMPSTE.filter <- subset(ZMPSTE.Singlet, idents = zmpste_keep)

# 34. PROCESS ZMPSTE SUBTYPES -----------------------------------------------
ZMPSTE.SPG <- subset(ZMPSTE.filter, idents = "SPG") %>% process_cell_subtype(resolution = 2)
ZMPSTE.Sertoli <- subset(ZMPSTE.filter, idents = "Sertoli") %>% process_cell_subtype(resolution = 0.15)
ZMPSTE.Macrophage <- subset(ZMPSTE.filter, idents = "Macrophage") %>% process_cell_subtype(resolution = 0.4)

# 35. ANNOTATE ZMPSTE SUBTYPES ----------------------------------------------
ZMPSTE.SPG@meta.data <- ZMPSTE.SPG@meta.data %>%
  mutate(celltype = case_when(
    RNA_snn_res.2 %in% c("8", "10", "0") ~ "Undifferentiated SSCs",
    RNA_snn_res.2 %in% c("15", "14", "1") ~ "Differentiated SSCs",
    TRUE ~ "PI"
  ))

ZMPSTE.Sertoli@meta.data <- ZMPSTE.Sertoli@meta.data %>%
  mutate(celltype = case_when(
    SCT_snn_res.0.15 == "2" ~ "Rhox5+ Sertoli",
    SCT_snn_res.0.15 == "1" ~ "Akap12+ Sertoli",
    SCT_snn_res.0.15 %in% c("0", "3") ~ "Defb36+ Sertoli",
    TRUE ~ as.character(SCT_snn_res.0.15)
  ))

ZMPSTE.Macrophage@meta.data <- ZMPSTE.Macrophage@meta.data %>%
  mutate(celltype = case_when(
    SCT_snn_res.0.4 == "2" ~ "Cd74+ macrophage",
    SCT_snn_res.0.4 == "1" ~ "Hmox1+ macrophage",
    SCT_snn_res.0.4 == "0" ~ "Crisp2+ macrophage",
    TRUE ~ as.character(SCT_snn_res.0.4)
  ))

# 36. MERGE ZMPSTE ANNOTATIONS ----------------------------------------------
ZMPSTE.filter <- merge_annotations(ZMPSTE.filter, ZMPSTE.SPG)
ZMPSTE.filter <- merge_annotations(ZMPSTE.filter, ZMPSTE.Macrophage)

# 37. SAVE ZMPSTE RESULTS ---------------------------------------------------
saveRDS(ZMPSTE.filter, file = "ZMPSTE.filter.rds")
saveRDS(ZMPSTE.SPG, file = "ZMPSTE.SPG.rds")
saveRDS(ZMPSTE.Sertoli, file = "ZMPSTE.Sertoli.rds")
saveRDS(ZMPSTE.Macrophage, file = "ZMPSTE.Macrophage.rds")

# ============================================================================
# FUNCTION TO PROCESS ALL DRUG TREATMENTS
# ============================================================================

#' Process a drug treatment dataset
#' @param treatment_name Name of treatment (e.g., "MTX", "FXT", "CDH", "DEHP", "CPA")
#' @param treatment_samples Vector of sample names for the treatment
#' @param control_samples Vector of control sample names
process_treatment <- function(treatment_name, treatment_samples, control_samples) {
  message("Processing ", treatment_name, " treatment...")
  
  # Subset data
  all_samples <- c(treatment_samples, control_samples)
  treatment_data <- subset(testis_total, idents = all_samples)
  treatment_data <- subset(treatment_data, subset = Doublet == "Singlet")
  
  # Integration
  integrated <- Integration_SCP(
    srtMerge = treatment_data, 
    batch = "orig.ident", 
    nHVF = 2000,
    integration_method = "Seurat"
  )
  
  # Add group information
  integrated@meta.data <- integrated@meta.data %>%
    mutate(group = case_when(
      orig.ident %in% treatment_samples ~ treatment_name,
      orig.ident %in% control_samples ~ "Normal",
      TRUE ~ NA_character_
    ))
  
  return(integrated)
}

# Control samples (can be reused)
control_samples <- c("CTL139", "CTL209", "Ctrl1", "Ctrl2", "Ctrl3")

# Process all treatments
MTX_samples <- c("MTX55", "MTX63", "MTX78")
FXT_samples <- c("FXT113", "FXT137", "FXT144")
CDH_samples <- c("CDH233", "CDH234", "CDH245")
DEHP_samples <- c("DEHP107", "DEHP109", "DEHP122")
CPA_samples <- c("C72", "C74", "C82")

MTX.combined <- process_treatment("MTX", MTX_samples, control_samples)
FXT.combined <- process_treatment("FXT", FXT_samples, control_samples)
CDH.combined <- process_treatment("CDH", CDH_samples, control_samples)
DEHP.combined <- process_treatment("DEHP", DEHP_samples, control_samples)
CPA.combined <- process_treatment("CPA", CPA_samples, control_samples)

# Save treatment results
saveRDS(MTX.combined, file = "MTX.combined.rds")
saveRDS(FXT.combined, file = "FXT.combined.rds")
saveRDS(CDH.combined, file = "CDH.combined.rds")
saveRDS(DEHP.combined, file = "DEHP.combined.rds")
saveRDS(CPA.combined, file = "CPA.combined.rds")

# ============================================================================
# PSEUDOTIME ANALYSIS FUNCTION
# ============================================================================

#' Run pseudotime analysis on SPG cells
run_pseudotime_analysis <- function(obj, group_var = "celltype", reduction = "umap") {
  obj <- RunSlingshot(
    srt = obj, 
    group.by = group_var, 
    reduction = reduction,
    start = "Undifferentiated SSCs"
  )
  
  obj <- RunDynamicFeatures(
    srt = obj, 
    lineages = c("Lineage1"),
    n_candidates = 2000
  )
  
  return(obj)
}

# Example usage (uncomment to run)
# alcohol.SPG <- run_pseudotime_analysis(alcohol.SPG)
# DynamicHeatmap(
#   srt = alcohol.SPG, lineages = c("Lineage1"),
#   use_fitted = TRUE, n_split = 3,
#   species = "Mus_musculus", db = "GO_BP",
#   heatmap_palette = "viridis", cell_annotation = "celltype"
# )

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

print_summary_stats <- function(obj, name) {
  cat("\n=== ", name, " ===\n")
  cat("Number of cells:", ncol(obj), "\n")
  cat("Number of genes:", nrow(obj), "\n")
  cat("Cell types:\n")
  print(table(Idents(obj)))
}

print_summary_stats(testis_total.DF, "All samples (after doublet removal)")
print_summary_stats(alcohol.final, "Alcohol treatment")
print_summary_stats(ZMPSTE.filter, "ZMPSTE24 KO")

cat("\nAnalysis complete!\n")
