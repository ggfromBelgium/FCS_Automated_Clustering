# ============================================================================
# SCRIPT D'ANALYSE CYTOMETRIE EN FLUX AUTOMATISÉ
# Gère les fichiers FCS contenant plusieurs échantillons (FileID / FolderID)
# Nouveauté: intergation d'une étape de diagnostic: normalisation nécéssaire? 
# Si oui (voir umap par condition et par fichier) alors CytoNorm
# Amélioration du prompt pour deepseek
# mise en cache des données de UMAP et tSNE pour le renomage et la fusion des clusters
# sécurisation de la clé API Deepseek
# ajout de la regle 7 de l etape 4 pour geré les marqueurs de linéages coexprimés abbérents
# modification du pipeline de détermination de k (nbr optimal de cluster)
#   -> ne pas hesiter a modifier les paramettre de cutreeDynamic en fonction des donnée
#     (method="hybrid/tree"; deepSplit = 2, minClusterSize = 5)
# remplacement du test chi carré par un autre test plus adapté
# ajout de test stat non apparié (pas correct si même patient pré-post)si plus de deux conditions avec la nomenclature Repondeurs_T0_Id, NonRepondeurs__T1_Id, ...
#
# les parametre du flowsom se trouve dans la SECTION 11: FONCTION PRINCIPALE D'ANALYSE
# ligne numéro 2425 environ
# ============================================================================
# Fuseau horaire de Bruxelles
Sys.setenv(TZ = "Europe/Brussels")

# ----- Récupération sécurisée de la clé API DeepSeek -----
api_key <- NULL

# 1) Essayer le fichier local (persistant après reboot)
key_file <- path.expand("~/.deepseek_api_key")
if (file.exists(key_file)) {
  api_key <- trimws(readLines(key_file, warn = FALSE))
  if (length(api_key) == 0 || nchar(api_key) == 0) api_key <- NULL
}

# 2) Sinon, lire depuis ~/.bashrc
if (is.null(api_key)) {
  bashrc_path <- path.expand("~/.bashrc")
  if (file.exists(bashrc_path)) {
    lines <- readLines(bashrc_path, warn = FALSE)
    key_line <- grep("^export DEEPSEEK_API_KEY=", lines, value = TRUE)
    if (length(key_line) > 0) {
      key_value <- sub("^export DEEPSEEK_API_KEY=", "", key_line[1])
      key_value <- gsub('"', '', key_value)
      key_value <- trimws(key_value)
      if (nchar(key_value) > 0) api_key <- key_value
    }
  }
}

# 3) Dernier recours : demande interactive + proposition de sauvegarde
if (is.null(api_key) || nchar(api_key) == 0) {
  cat("Aucune clé API DeepSeek trouvée.\n")
  cat("Entrez votre clé (ou laissez vide pour annotation manuelle) : ")
  api_key_interactive <- readline()
  if (nchar(trimws(api_key_interactive)) > 0) {
    api_key <- api_key_interactive
    cat("Sauvegarder cette clé dans ~/.deepseek_api_key pour les prochaines exécutions ? (y/n) : ")
    if (tolower(readline()) == "y") {
      writeLines(api_key, key_file)
      Sys.chmod(key_file, "600")
      cat("Clé sauvegardée dans", key_file, "\n")
    }
  } else {
    api_key <- ""
  }
}

Sys.setenv(DEEPSEEK_API_KEY = api_key)

# ============================================================================
# SECTION 1: CHARGEMENT DES LIBRARIES
# ============================================================================
suppressPackageStartupMessages({
  library(flowCore)
  library(FlowSOM)
  library(uwot)
  library(dplyr)
  library(ComplexHeatmap)
  library(Cairo)
  library(grid)
  library(ggpubr)
  library(cluster)
  library(viridis)
  library(ConsensusClusterPlus)
  library(ggplot2)
  library(reshape2)
  library(ggalluvial)
  library(data.table)
  library(httr)     
  library(jsonlite)
  library(scales)
  library(RColorBrewer)
  library(digest)
  library(keyring)
  library(dynamicTreeCut)
  library(fpc) 
})
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
Sys.setenv(TF_XLA_FLAGS = "--tf_xla_enable_xla_devices=false")
options(tensorflow.quiet = TRUE)
options(uwot.verbose = FALSE)
options(OutDec = ".")

# ============================================================================
# SECTION 2: FONCTIONS UTILITAIRES DE BASE
# ============================================================================
get_marker_name <- function(channel, marker_descriptions) {
  if (!is.na(marker_descriptions[channel]) && marker_descriptions[channel] != "")
    marker_descriptions[channel] else channel
}

get_file_name <- function(file_id, desc) {
  key <- paste0("File", file_id, "_Name")
  if (key %in% names(desc))
    gsub("\\.fcs$", "", desc[[key]], ignore.case = TRUE) else paste0("File_", file_id)
}

natural_sort <- function(x) {
  non_numeric <- gsub("[0-9]+", "", x)
  numeric <- as.numeric(gsub("[^0-9]+", "", x))
  if (all(is.na(numeric))) return(sort(x))
  x[order(non_numeric, numeric, na.last = TRUE)]
}

create_plot_window <- function(width = 8, height = 6) {
  graphics.off()
  if (.Platform$OS.type == "windows") windows(width, height)
  else if (capabilities("X11")) X11(width, height)
  else if (capabilities("aqua")) quartz(width, height)
  else tryCatch(x11(width, height), error = function(e) message("Impossible d'ouvrir une fenêtre graphique"))
}

open_pdf <- function(pdf_file) {
  if (file.exists(pdf_file)) {
    system2(if (.Platform$OS.type == "windows") "open" else "xdg-open",
            args = shQuote(pdf_file), wait = FALSE,
            stdout = if (.Platform$OS.type != "windows") "/dev/null" else FALSE,
            stderr = if (.Platform$OS.type != "windows") "/dev/null" else FALSE)
    cat("Ouverture du PDF :", pdf_file, "\n")
    return(TRUE)
  } else {
    warning("Le fichier PDF n'existe pas : ", pdf_file)
    return(FALSE)
  }
}

# ============================================================================
# SECTION 3: FONCTIONS D'APPEL API
# ============================================================================
call_deepseek <- function(prompt, api_key) {
  url <- "https://api.deepseek.com/v1/chat/completions"
  body <- list(
    model = "deepseek-flash", # deepseek-v4-pro, deepseek-flash
    messages = list(list(role = "user", content = prompt)),
    temperature = 0.2,
    max_tokens = 50000,
    thinking = list(type = "enabled") #disabled, enabled
  )
  response <- tryCatch({
    httr::POST(url,
               httr::add_headers("Authorization" = paste("Bearer", api_key),
                                 "Content-Type" = "application/json"),
               body = jsonlite::toJSON(body, auto_unbox = TRUE),
               encode = "raw")
  }, error = function(e) { message("❌ Erreur réseau / HTTP : ", e$message); return(NULL) })
  if (is.null(response)) return(NULL)
  if (httr::status_code(response) != 200) {
    message("❌ API DeepSeek a retourné une erreur : ", httr::status_code(response))
    message("   Réponse : ", httr::content(response, as = "text", encoding = "UTF-8"))
    return(NULL)
  }
  content_raw <- httr::content(response, as = "text", encoding = "UTF-8")
  parsed <- jsonlite::fromJSON(content_raw, simplifyVector = FALSE)
  if (!is.null(parsed$choices) && length(parsed$choices) > 0 &&
      !is.null(parsed$choices[[1]]$message$content)) {
    return(parsed$choices[[1]]$message$content)
  } else {
    message("❌ Structure de réponse inattendue : impossible d'extraire le contenu.")
    return(NULL)
  }
}

parse_json_with_fallback <- function(response_text, n_clusters) {
  json_str <- response_text
  json_str <- gsub("^```json\\s*", "", json_str)
  json_str <- gsub("^```\\s*", "", json_str)
  json_str <- gsub("\\s*```$", "", json_str)
  json_str <- trimws(json_str)
  
  first_brace <- regexpr("\\{", json_str)
  last_brace <- rev(gregexpr("\\}", json_str)[[1]])[1]
  if (first_brace == -1 || last_brace == -1) stop("Aucune accolade trouvée")
  json_str <- substr(json_str, first_brace, last_brace)
  
  parsed <- tryCatch(
    jsonlite::fromJSON(json_str, simplifyVector = FALSE),
    error = function(e) stop("JSON invalide : ", e$message)
  )
  
  cluster_keys <- grep("^cluster_\\d+$", names(parsed), value = TRUE, ignore.case = TRUE)
  if (length(cluster_keys) == 0) stop("Aucune clé 'cluster_X' trouvée")
  
  cluster_indices <- as.numeric(gsub("cluster_", "", cluster_keys, ignore.case = TRUE))
  ordered_keys <- cluster_keys[order(cluster_indices)]
  
  ann <- character(n_clusters)
  conf <- character(n_clusters)
  for (k in seq_len(min(length(ordered_keys), n_clusters))) {
    key <- ordered_keys[k]
    value <- parsed[[key]]
    if (is.list(value) && "annotation" %in% names(value)) {
      ann[k] <- as.character(value[["annotation"]])
      conf[k] <- if (!is.null(value[["confidence"]])) as.character(value[["confidence"]]) else "MEDIUM"
    } else if (is.character(value)) {
      ann[k] <- value
      conf[k] <- "MEDIUM"
    } else {
      ann[k] <- "Unknown"
      conf[k] <- "LOW"
    }
  }
  empty <- which(ann == "" | is.na(ann))
  ann[empty] <- "Unknown"
  conf[empty] <- "LOW"
  
  list(annotations = ann, confidence = conf)
}

# ============================================================================
# SECTION 4: ANNOTATION AUTOMATISÉE
# ============================================================================

annotate_clusters_improved <- function(cluster_medians, marker_descriptions, output_dir = NULL,
                                       api_key = NULL, population_description = NULL) {
  message("\n--- Début du processus d'annotation ---")
  
  clean_names <- sapply(colnames(cluster_medians), function(cn) get_marker_name(cn, marker_descriptions))
  colnames(cluster_medians) <- clean_names
  
  tech_cols <- c("Time", "FSC.A", "FSC.H", "FSC.W", "SSC.A", "SSC.H", "SSC.W", "Viability", "AF.A")
  bio_markers <- setdiff(colnames(cluster_medians), tech_cols)
  cluster_medians <- cluster_medians[, bio_markers, drop = FALSE]
  n_clusters <- nrow(cluster_medians)
  
  if (is.null(api_key)) api_key <- Sys.getenv("DEEPSEEK_API_KEY")
  if (is.null(api_key) || nchar(api_key) == 0) {
    message("⚠️  Aucune clé API DeepSeek. Annotation automatique impossible.")
    generic <- paste0("Cluster_", 1:n_clusters)
    return(list(
      annotations = generic,
      scores = rep(NA, n_clusters),
      details = lapply(generic, function(a) list(annotation = a)),
      source = "manual_fallback",
      confidence = rep("LOW", n_clusters),
      confidence_score = rep(0.0, n_clusters)
    ))
  }
  
  global_min <- min(cluster_medians, na.rm = TRUE)
  global_max <- max(cluster_medians, na.rm = TRUE)
  if (global_max == global_min) {
    cluster_medians_scaled <- matrix(0, nrow = n_clusters, ncol = ncol(cluster_medians))
    colnames(cluster_medians_scaled) <- colnames(cluster_medians)
  } else {
    cluster_medians_scaled <- round(((cluster_medians - global_min) / (global_max - global_min)) * 10, 1)
  }
  cluster_medians_scaled <- as.data.frame(cluster_medians_scaled)
  
  cluster_table <- paste0(
    "cluster\t", paste(colnames(cluster_medians_scaled), collapse = "\t"), "\n",
    paste(sapply(1:n_clusters, function(i) {
      paste(i, paste(cluster_medians_scaled[i, ], collapse = "\t"), sep = "\t")
    }), collapse = "\n")
  )
  
  gating_context <- if (!is.null(population_description) && nchar(trimws(population_description)) > 0) {
    paste0("All cells were pre-gated. Selected population: ", population_description, ".")
  } else {
    "All cells are from a pre-gated population (no specific description provided)."
  }
  
  prompt <- paste0(
    "You are an expert immunologist specializing in flow cytometry data analysis and immune cell phenotyping. ",
    "Your task is to annotate cytometry clusters based on marker expression profiles.\n\n",
    "═══════════════════════════════════════════════════════════════\n",
    "PHASE 1: CONTEXT & PRE-GATING\n",
    "═══════════════════════════════════════════════════════════════\n",
    gating_context, "\n\n",
    "Based on the pre-gating described above:\n",
    "- Determine which broad lineages are INCLUDED in this dataset.\n",
    "- Determine which populations are EXCLUDED and must NEVER appear in your annotations.\n",
    "- State these exclusions explicitly before proceeding.\n\n",
    
    "═══════════════════════════════════════════════════════════════\n",
    "PHASE 2: DATA SCALE & DYNAMIC THRESHOLD CALIBRATION\n",
    "═══════════════════════════════════════════════════════════════\n",
    "DATA SCALE: 0 = No expression, 10 = Maximum expression in this sample.\n\n",
    "IMPORTANT — Threshold determination strategy:\n",
    "Rather than using fixed cutoffs, you MUST determine positivity/negativity thresholds DYNAMICALLY for EACH marker by analyzing its distribution across ALL clusters:\n",
    "1. For each marker column, identify the RANGE (min → max) across all clusters.\n",
    "2. If a marker shows a BIMODAL distribution (clear gap between low and high values), use that gap as the threshold.\n",
    "3. If a marker is uniformly low across all clusters, it is not discriminatory — do not use it for annotation decisions.\n",
    "4. If a marker is uniformly high across all clusters, it is not discriminatory — do not use it for annotation decisions.\n",
    "5. As a fallback guideline (only when distribution analysis is ambiguous):\n",
    "   - 0–2.0 → Likely NEGATIVE\n",
    "   - 2.1–4.0 → LOW/DIM (context-dependent)\n",
    "   - 4.1–6.0 → INTERMEDIATE\n",
    "   - 6.1–8.0 → HIGH (clearly positive)\n",
    "   - 8.1–10 → VERY HIGH\n\n",
    "Before annotating, produce a short \"Marker Calibration Summary\" listing each marker, its observed range, and your chosen positive/negative threshold.\n\n",
    
    "═══════════════════════════════════════════════════════════════\n",
    "PHASE 3: UNIVERSAL REFERENCE PHENOTYPE DATABASE\n",
    "═══════════════════════════════════════════════════════════════\n",
    "Below is a comprehensive phenotype reference database. For each population, key markers are listed.\n",
    "You must MATCH clusters against this database using ONLY the markers present in the data.\n",
    "If a reference marker is NOT in the panel, IGNORE it (do not penalize or reward its absence).\n\n",
    
    "────────────────────────────\n",
    "A. T CELL LINEAGE\n",
    "   Lineage requirement: CD3+ (and/or TCRab+/TCRgd+ if available)\n",
    "────────────────────────────\n",
    "| Population | Expected Phenotype |\n",
    "|---|---|\n",
    "| CD4+ T cells (Th0/resting) | CD3+, CD4+, CD8a-, TCRgd-, Tbet_low, GATA3_low, RORgt-, Foxp3- |\n",
    "| CD4+ Th1 cells | CD3+, CD4+, CD8a-, Tbet_HIGH, GATA3_low, RORgt-, IFNg+ (if available) |\n",
    "| CD4+ Th2 cells | CD3+, CD4+, CD8a-, GATA3_HIGH, Tbet_low, RORgt-, IL4+ (if available) |\n",
    "| CD4+ Th17 cells | CD3+, CD4+, CD8a-, RORgt_HIGH, Tbet_low, CCR6+ (if available) |\n",
    "| CD4+ Th1/17 cells | CD3+, CD4+, CD8a-, RORgt+, Tbet+, dual expression |\n",
    "| CD4+ Th9 cells | CD3+, CD4+, PU.1+ (if available), IRF4+ (if available) |\n",
    "| CD4+ Th22 cells | CD3+, CD4+, AHR+ (if available), CCR10+ (if available) |\n",
    "| CD4+ Tfh cells | CD3+, CD4+, CXCR5+ (if available), PD1+ (if available), Bcl6+ (if available) |\n",
    "| CD4+ Tregs | CD3+, CD4+, CD8a-, Foxp3_HIGH, CD127_low/neg, CD25_HIGH (if available) |\n",
    "| CD4+ Tregs (activated/effector) | CD3+, CD4+, Foxp3_HIGH, CD127_low, GATA3+ or Tbet+ possible, ICOS+ (if available) |\n",
    "| CD4+ Th1-like Tregs | CD3+, CD4+, Foxp3_HIGH, Tbet+, CD127_low |\n",
    "| CD4+ Tfr cells | CD3+, CD4+, Foxp3_HIGH, CXCR5+ (if available) |\n",
    "| CD8+ T cells (resting/Tc0) | CD3+, CD8a+, CD4-, TCRgd- |\n",
    "| CD8+ Tc1 cells | CD3+, CD8a+, CD4-, Tbet_HIGH, Granzyme+ (if available) |\n",
    "| CD8+ Tc2 cells | CD3+, CD8a+, CD4-, GATA3_HIGH |\n",
    "| CD8+ Tc17 cells | CD3+, CD8a+, CD4-, RORgt_HIGH |\n",
    "| CD8+ regulatory T cells | CD3+, CD8a+, Foxp3+, CD127_low |\n",
    "| DN T cells (Double Negative) | CD3+, CD4-, CD8a-, TCRgd- |\n",
    "| DP T cells (Double Positive) | CD3+, CD4+, CD8a+ |\n",
    "| γδ T cells | CD3+, TCRgd_HIGH, CD4- usually, CD8a- usually |\n",
    "| γδ T cells (Vγ9Vδ2) | CD3+, TCRgd+, specific Vg/Vd staining (if available) |\n",
    "| MAIT cells | CD3+, TCRVa7.2+ (if available), CD161+ (if available), CD8+ or DN |\n",
    "| NKT cells (type I / iNKT) | CD3+, CD49b+ or NK1.1+ (if available), NKp46-/low, CD1d-tet+ (if available) |\n",
    "| NKT cells (type II) | CD3+, CD49b+ or NK1.1+, NKp46-, CD1d-tet- |\n",
    "| CD4+ NKT-like cells | CD3+, CD4+, CD49b_HIGH, NKp46_low/neg |\n",
    "| DN NKT-like cells | CD3+, CD4-, CD8a-, CD49b_HIGH |\n\n",
    
    "Subsets by differentiation state (apply as suffix if markers available):\n",
    "| Naive | CD44_low (or CD45RA+), CD62L_HIGH, CCR7+ (if available) |\n",
    "| TCM (Central Memory) | CD44_HIGH, CD62L+, CCR7+ (if available) |\n",
    "| TEM (Effector Memory) | CD44_HIGH, CD62L-, CCR7- (if available) |\n",
    "| TEMRA (Terminally differentiated) | CD44_HIGH, CD62L-, CD45RA+ (if available), KLRG1+ (if available) |\n",
    "| Effector | CD44_HIGH, CD62L-, cytokine producing |\n",
    "| Exhausted T cells | PD1+, Tim3+, LAG3+ (if available) |\n",
    "| Tissue-resident memory (Trm) | CD69+ (if available), CD103+ (if available) |\n\n",
    
    "────────────────────────────\n",
    "B. B CELL LINEAGE\n",
    "   Lineage requirement: B220+ (and/or CD19+, CD20+ if available), CD3-\n",
    "────────────────────────────\n",
    "| Population | Expected Phenotype |\n",
    "|---|---|\n",
    "| B cells (resting/naive) | B220+, CD3-, CD4-, CD8a-, IgD+ (if available), CD27- (if available) |\n",
    "| B cells (activated) | B220+, CD3-, CD69+ or GL7+ (if available) |\n",
    "| Germinal center B cells | B220+, GL7+ (if available), Fas+ (if available) |\n",
    "| Memory B cells | B220+, CD27+ (if available), IgD- (if available) |\n",
    "| B-1a cells | B220_low/+, CD5+ (if available), CD43+ (if available) |\n",
    "| B-1b cells | B220_low/+, CD5- (if available) |\n",
    "| Marginal zone B cells | B220+, CD21+ (if available), CD23- (if available) |\n",
    "| Plasmablasts | B220_low/neg, CD138+ (if available), BLIMP1+ (if available) |\n",
    "| Plasma cells | B220_neg, CD138+ (if available) |\n",
    "| Age-associated B cells (ABCs) | B220+, CD3-, Tbet+, CD11c+ (if available) |\n",
    "| Breg (regulatory B cells) | B220+, IL10+ (if available), CD1d+ (if available) |\n",
    "| Transitional B cells | B220+, AA4.1+ (if available), IgM_HIGH (if available) |\n\n",
    
    "────────────────────────────\n",
    "C. NK CELLS\n",
    "   Lineage requirement: CD3-, NKp46+ and/or NK1.1+ (if available), CD49b+\n",
    "────────────────────────────\n",
    "| Population | Expected Phenotype |\n",
    "|---|---|\n",
    "| NK cells | CD3-, NKp46+, CD49b+, B220-, Tbet+ |\n",
    "| NK cells (immature) | CD3-, NKp46+, CD49b_low/neg, CD27+ (if available) |\n",
    "| NK cells (mature) | CD3-, NKp46+, CD49b+, Tbet_HIGH, CD27- (if available) |\n",
    "| NK cells (cytokine-producing) | CD3-, NKp46+, CD49b+, IFNg+ (if available) |\n",
    "| NK cells (cytotoxic) | CD3-, NKp46+, CD49b+, Granzyme+ (if available), Perforin+ (if available) |\n",
    "| Uterine NK / tissue NK | CD3-, NKp46+, CD49a+ (if available), CD49b- |\n\n",
    
    "────────────────────────────\n",
    "D. INNATE LYMPHOID CELLS (ILCs)\n",
    "   Lineage requirement: CD3-, CD127+ (IL-7Ra), Lineage-negative\n",
    "────────────────────────────\n",
    "| Population | Expected Phenotype |\n",
    "|---|---|\n",
    "| ILC1 | CD3-, CD127+, Tbet+, NKp46+/-, RORgt-, GATA3_low, B220-, CD49b- or low, CD49a+ (if available) |\n",
    "| ILC2 | CD3-, CD127+, GATA3_HIGH, RORgt-, Tbet-, NKp46-, B220-, ICOS+ (if available), ST2+ (if available) |\n",
    "| ILC3 (NCR+) | CD3-, CD127+, RORgt_HIGH, NKp46+, Tbet_low |\n",
    "| ILC3 (NCR-) | CD3-, CD127+, RORgt_HIGH, NKp46-, CD4+/- |\n",
    "| ILC3 (LTi-like) | CD3-, CD127+, RORgt_HIGH, CD4+ (if available), CCR6+ (if available) |\n",
    "| ILCp (progenitor) | CD3-, CD127+, low TF expression, PLZF+ (if available) |\n",
    "| Ex-ILC3 / ILC1-like | CD3-, CD127+, Tbet+, RORgt_low (formerly RORgt+) |\n\n",
    
    "────────────────────────────\n",
    "E. DENDRITIC CELLS\n",
    "   Lineage requirement: CD3-, B220-/+ depending on subset\n",
    "────────────────────────────\n",
    "| Population | Expected Phenotype |\n",
    "|---|---|\n",
    "| pDC (plasmacytoid DC) | CD3-, B220+, CD4-/low, CD8a-, PDCA1+ (if available), SiglecH+ (if available), CD11c_low (if available) |\n",
    "| cDC1 | CD3-, B220-, CD8a+ or CD103+ (if available), XCR1+ (if available), CLEC9A+ (if available), CD11c+ (if available) |\n",
    "| cDC2 | CD3-, B220-, CD11b+ (if available), CD4+/-, SIRPa+ (if available), CD11c+ (if available) |\n",
    "| moDC (monocyte-derived DC) | CD3-, CD11c+ (if available), CD11b+ (if available), variable |\n",
    "| Migratory DC | CD3-, MHC-II_HIGH (if available), CCR7+ (if available) |\n\n",
    
    "────────────────────────────\n",
    "F. MYELOID / OTHER\n",
    "────────────────────────────\n",
    "| Population | Expected Phenotype |\n",
    "|---|---|\n",
    "| Classical monocytes | CD14+, CD16- (if available), CD11b+ (if available) |\n",
    "| Non-classical monocytes | CD14_low/neg, CD16+ (if available), CD11b+ (if available) |\n",
    "| Intermediate monocytes | CD14+, CD16+ (if available) |\n",
    "| Neutrophils | Ly6G+ (if available), CD11b+ (if available), CD3-, B220- |\n",
    "| Eosinophils | SiglecF+ (if available), CD11b+ (if available) |\n",
    "| Basophils | CD3-, B220-, FcεRI+ (if available), CD49b+ possible, CD200R3+ (if available) |\n",
    "| Mast cells | FcεRI+ (if available), c-Kit+ (if available) |\n",
    "| Macrophages | F4/80+ (if available), CD11b+ (if available), CD14+ possible |\n",
    "| MDSCs (myeloid-derived suppressor) | CD11b+ (if available), Gr1+ (if available), variable |\n\n",
    
    "────────────────────────────\n",
    "G. OTHER / RARE\n",
    "────────────────────────────\n",
    "| Population | Expected Phenotype |\n",
    "|---|---|\n",
    "| Hematopoietic stem/progenitor cells | Lin-, c-Kit+ (if available), Sca1+ (if available) |\n",
    "| Erythroid progenitors | Ter119+ (if available), CD71+ (if available) |\n",
    "| Megakaryocyte progenitors | CD41+ (if available) |\n",
    "| Unknown [marker]+ | Fallback when no confident match exists |\n\n",
    
    "═══════════════════════════════════════════════════════════════\n",
    "═══════════════════════════════════════════════════════════════\n",
    "RULE 1 — PANEL AWARENESS:\n",
    "- Identify ALL markers present in the data table. These are your ONLY available markers.\n",
    "- Any reference phenotype marker NOT in the panel is INVISIBLE: you cannot confirm or deny it.\n",
    "- A population is a candidate ONLY if its key MEASURABLE markers (those present in your panel) match.\n\n",
    "RULE 2 — HIERARCHICAL ANNOTATION (follow this order strictly):\n",
    "  Step 1: LINEAGE ASSIGNMENT\n",
    "    - Use lineage-defining markers (e.g., CD3 for T cells, B220/CD19 for B cells, NKp46/NK1.1 for NK, CD127 for ILCs, etc.)\n",
    "    - A cluster belongs to one primary lineage. Assign it.\n",
    "  Step 2: SUBSET REFINEMENT\n",
    "    - Within the assigned lineage, use transcription factors (Tbet, GATA3, RORgt, Foxp3, etc.), \n",
    "      co-receptors (CD4, CD8a, TCRgd), activation markers, and other available markers.\n",
    "  Step 3: DIFFERENTIATION STATE (if applicable markers are present)\n",
    "    - Determine Naive/TCM/TEM/TEMRA/Effector/Exhausted status if memory/activation markers are in the panel.\n\n",
    "RULE 3 — EXCLUSION ENFORCEMENT:\n",
    "- NEVER annotate a cluster as a population that was EXCLUDED by pre-gating.\n",
    "- Before finalizing each annotation, re-check it against the pre-gating exclusion list from Phase 1.\n\n",
    "RULE 4 — SELECTION CONSTRAINT:\n",
    "- You MUST select a cell identity from the Phase 3 reference database.\n",
    "- If no population matches satisfactorily, use: \"Unknown [most prominent marker(s)]+\"\n\n",
    "RULE 5 — AMBIGUITY HANDLING:\n",
    "- If two populations are equally plausible, list both as: \"Likely [Population A] or [Population B]\"\n",
    "- Assign LOW confidence in such cases.\n\n",
    "RULE 6 — BIOLOGICAL PLAUSIBILITY:\n",
    "- Consider whether your annotation makes biological sense in the context of the tissue/sample type.\n",
    "- Consider expected relative frequencies: if you annotate 8 clusters as the same rare cell type, reconsider.\n\n",
    "RULE 7 — RESOLVING AMBIGUOUS CO‑EXPRESSION & MIXED LINEAGE SIGNALS\n",
    "- When a cluster expresses mutually exclusive lineage markers\n",
    "  (e.g., CD3 AND B220 AND CD127 high), follow this strict\n",
    "  tie‑breaking procedure:\n\n",
    "  Step 1. IDENTIFY THE BI-MODALITY CHAMPION\n",
    "   - For each ambiguous marker, refer back to the Marker Calibration Summary (Phase 2).\n",
    "   - Pick the marker whose distribution across ALL clusters shows the CLEAREST bimodal\n",
    "     separation (largest gap between negative and positive peaks, highest fold‑change).\n",
    "   - This marker has the strongest “biological signal” in this specific sample.\n\n",
    "  Step 2. EVALUATE BIOLOGICAL SPECIFICITY\n",
    "   - Among the remaining candidates, favour the marker with the highest lineage\n",
    "     specificity in the given context (e.g., CD127 is highly specific for ILCs;\n",
    "   - B220 for B/plasmacytoid DC; CD3 for T cells). Prefer a marker that is\n",
    "     virtually never co‑expressed with the others in normal physiology over a\n",
    "     marker that can be promiscuous.\n\n",
    "  Step 3. ASSIGN THE LINEAGE BASED ON THE CHAMPION\n",
    "   - Let the marker selected in steps 1‑2 determine the primary lineage.\n",
    "   - Interpret the other intermediate markers as either:\n",
    "      - low-level non‑specific staining,\n",
    "      - activation‑induced upregulation (e.g., CD3 on activated NK cells is rare\n",
    "        but possible, while CD127 on T cells is common in memory subsets),\n",
    "      - or a true minor co‑expression (document it as an atypical feature).\n\n",
    "  Step 4. DOCUMENT THE DECISION EXPLICITLY\n",
    "   - In the annotation reasoning, state:\n",
    "     “Ambiguous co‑expression resolved by prioritising [marker X] because of its\n",
    "      clear bimodal separation (gap = Y) and high lineage specificity. The\n",
    "      intermediate [marker Z] value is treated as background/atypical co‑expression.”\n\n",
    
    "═══════════════════════════════════════════════════════════════\n",
    "PHASE 5: MARKER CALIBRATION & INITIAL ANNOTATION\n",
    "═══════════════════════════════════════════════════════════════\n",
    "Perform the following steps IN ORDER. Show your work.\n\n",
    "STEP 5.1 — Marker Calibration Summary:\n",
    "For each marker in the panel, report:\n",
    "| Marker | Min value | Max value | Observed distribution pattern | Chosen POS/NEG threshold |\n\n",
    "STEP 5.2 — Pre-gating Exclusion List:\n",
    "Based on Phase 1 context, list populations that CANNOT appear in annotations:\n",
    "- [Population 1]: excluded because [reason]\n",
    "- ...\n\n",
    "STEP 5.3 — Cluster-by-Cluster Annotation:\n",
    "For each cluster, produce:\n",
    "| Cluster | Marker Summary (POS/NEG) | Lineage | Subset Evidence | Candidate Annotation | Confidence |\n\n",
    
    "═══════════════════════════════════════════════════════════════\n",
    "PHASE 6: MANDATORY CONTRADICTION AUDIT (CRITICAL)\n",
    "═══════════════════════════════════════════════════════════════\n",
    "⚠️ DO NOT produce the final JSON output until this entire phase is completed for ALL clusters.\n\n",
    "You must now systematically re-read EVERY annotation against the RAW DATA and check for contradictions. This is a FULL AUDIT, not a spot-check.\n\n",
    "For EACH cluster, execute ALL of the following checks:\n\n",
    "────────────────────────────\n",
    "CHECK 1 — LINEAGE MARKER CONTRADICTION\n",
    "────────────────────────────\n",
    "Re-read the RAW expression values for lineage markers.\n",
    "| Your annotation | Required marker condition | Actual value in data | PASS/FAIL |\n\n",
    "Examples of contradictions to catch:\n",
    "- ❌ Annotated as T cell but CD3 < positive threshold\n",
    "- ❌ Annotated as B cell but B220/CD19 < positive threshold\n",
    "- ❌ Annotated as NK cell but CD3 > positive threshold (unless NKT)\n",
    "- ❌ Annotated as ILC but CD3 > positive threshold\n",
    "- ❌ Annotated as non-T cell but CD3 is clearly positive\n",
    "→ If FAIL: REVISE the annotation immediately. State the correction.\n\n",
    "────────────────────────────\n",
    "CHECK 2 — SUBSET MARKER CONTRADICTION\n",
    "────────────────────────────\n",
    "Re-read the RAW expression values for subset-defining markers.\n",
    "| Your annotation | Required subset marker | Expected expression | Actual value | PASS/FAIL |\n\n",
    "Examples of contradictions to catch:\n",
    "- ❌ Annotated as CD4+ T cell but CD4 < positive threshold\n",
    "- ❌ Annotated as CD8+ T cell but CD8a < positive threshold\n",
    "- ❌ Annotated as Treg but Foxp3 < positive threshold\n",
    "- ❌ Annotated as Th17/ILC3 but RORgt < positive threshold\n",
    "- ❌ Annotated as ILC2 but GATA3 is not the dominant transcription factor\n",
    "- ❌ Annotated as γδ T cell but TCRgd < positive threshold\n",
    "- ❌ Annotated as NK cell but NKp46 < positive threshold AND CD49b < positive threshold\n",
    "→ If FAIL: REVISE the annotation immediately. State the correction.\n\n",
    "────────────────────────────\n",
    "CHECK 3 — EXCLUSION MARKER CONTRADICTION\n",
    "────────────────────────────\n",
    "Re-read markers that should be NEGATIVE for the annotated population.\n",
    "| Your annotation | Should be negative for | Actual value | PASS/FAIL |\n\n",
    "Examples of contradictions to catch:\n",
    "- ❌ Annotated as CD4+ T cell but CD8a is clearly positive (unless DP T)\n",
    "- ❌ Annotated as B cell but CD3 is clearly positive\n",
    "- ❌ Annotated as NK cell but B220 is clearly positive (unless specific subset)\n",
    "- ❌ Annotated as CD8+ T cell but CD4 is clearly positive (unless DP T)\n",
    "→ If FAIL: REVISE or reclassify (e.g., consider DP T cells). State the correction.\n\n",
    "────────────────────────────\n",
    "CHECK 4 — PRE-GATING VIOLATION CHECK\n",
    "────────────────────────────\n",
    "Re-verify that NO annotation corresponds to a population excluded by pre-gating.\n",
    "| Your annotation | Is it excluded by pre-gating? | PASS/FAIL |\n",
    "→ If FAIL: REVISE the annotation immediately.\n\n",
    "────────────────────────────\n",
    "CHECK 5 — DUPLICATE / FREQUENCY PLAUSIBILITY CHECK\n",
    "────────────────────────────\n",
    "List all final annotations. Check for:\n",
    "- Are any two clusters given IDENTICAL annotations? \n",
    "  → If yes: Compare their marker profiles side-by-side. Can they be biologically distinguished? \n",
    "  → If they are truly indistinguishable with available markers, keep the same name and note \"similar profile.\"\n",
    "  → If they differ on key markers, revise one or both annotations.\n",
    "- Is any very rare population annotated for multiple clusters? \n",
    "  → Flag and verify.\n\n",
    "────────────────────────────\n",
    "CHECK 6 — GLOBAL CONSISTENCY REVIEW\n",
    "────────────────────────────\n",
    "Review all annotations as a SET:\n",
    "- Does the overall immune landscape make biological sense?\n",
    "- Are expected major populations represented (if their markers are in the panel)?\n",
    "- Are there any surprising absences or over-representations?\n",
    "- Note any concerns (this does not necessarily require changes, but flags for the user).\n\n",
    "────────────────────────────\n",
    "AUDIT SUMMARY:\n",
    "────────────────────────────\n",
    "After completing ALL checks for ALL clusters, produce:\n",
    "| Cluster | Original Annotation | Checks Passed | Checks Failed | Revised Annotation (if any) | Final Confidence |\n\n",
    "Total contradictions found: X\n",
    "Total revisions made: Y\n\n",
    
    "═══════════════════════════════════════════════════════════════\n",
    "PHASE 7: FINAL OUTPUT (Only after Phase 6 is 100% complete)\n",
    "═══════════════════════════════════════════════════════════════\n",
    "⚠️ ONLY produce this output AFTER completing the full Phase 6 audit for ALL clusters.\n\n",
    "**PART A — Detailed Reasoning Table:**\n",
    "| Cluster | Key Positive Markers | Key Negative Markers | Lineage | Annotation | Confidence | Audit Status |\n\n",
    "**PART B — Final JSON (strictly valid):**\n",
    "{\n",
    "  \"cluster_1\": {\"annotation\": \"Cell Type\", \"confidence\": \"HIGH/MEDIUM/LOW\"},\n",
    "  \"cluster_2\": {\"annotation\": \"Cell Type\", \"confidence\": \"HIGH/MEDIUM/LOW\"},\n",
    "  ...\n",
    "}\n\n",
    "**PART C — Audit Report Summary:**\n",
    "- Contradictions detected and corrected: [list]\n",
    "- Remaining ambiguities: [list]\n",
    "- Populations not detected in this dataset (expected but absent): [list]\n",
    "- Suggestions for panel improvement (optional): [list]\n\n",
    
    "--- MARKER DATA TABLE ---\n",
    cluster_table, "\n\n",
    "Now, execute all phases sequentially and produce the final output."
  )
  
  if (!is.null(output_dir)) {
    writeLines(prompt, file.path(output_dir, "last_prompt_sent.txt"))
  }
  
  message("📡 Envoi de la requête à l'IA...")
  response_text <- tryCatch({
    res <- call_deepseek(prompt, api_key)
    if (is.null(res) || res == "") stop("Réponse vide")
    res
  }, error = function(e) {
    message("❌ ERREUR API : ", e$message)
    return(NULL)
  })
  
  if (!is.null(output_dir)) {
    writeLines(if (!is.null(response_text)) response_text else "Échec de l'appel",
               file.path(output_dir, "IA_raw_response.txt"))
  }
  
  if (is.null(response_text)) {
    generic <- paste0("Cluster_", 1:n_clusters)
    return(list(
      annotations = generic,
      scores = rep(NA, n_clusters),
      details = lapply(generic, function(a) list(annotation = a)),
      source = "generic",
      confidence = rep("LOW", n_clusters),
      confidence_score = rep(0.0, n_clusters)
    ))
  }
  
  parsed <- tryCatch(
    parse_json_with_fallback(response_text, n_clusters),
    error = function(e) {
      message("⚠️ Échec du parsing : ", e$message)
      NULL
    }
  )
  
  if (is.null(parsed)) {
    generic <- paste0("Cluster_", 1:n_clusters)
    return(list(
      annotations = generic,
      scores = rep(NA, n_clusters),
      details = lapply(generic, function(a) list(annotation = a)),
      source = "generic",
      confidence = rep("LOW", n_clusters),
      confidence_score = rep(0.0, n_clusters)
    ))
  }
  
  annotations <- parsed$annotations
  conf_char <- parsed$confidence
  conf_score <- sapply(conf_char, function(c) switch(toupper(c), HIGH=0.9, MEDIUM=0.5, LOW=0.3, 0.3))
  
  message("✅ Annotations obtenues avec confiance :\n",
          paste(paste0("  Cluster ", 1:n_clusters, " : ", annotations, " [", conf_char, "]"), collapse = "\n"))
  
  list(
    annotations = annotations,
    scores = conf_score,
    details = lapply(seq_along(annotations), function(i)
      list(annotation = annotations[i], confidence = conf_char[i])),
    source = "AI",
    confidence = conf_char,
    confidence_score = conf_score
  )
}

# ============================================================================
# SECTION 5: FONCTIONS DE VISUALISATION
# ============================================================================
show_single_channel_histogram <- function(expr_data, channel, description, a, b, sample_name, title_suffix = "", threshold = NULL) {
  transformed_data <- asinh(a + b * expr_data[, channel])
  marker_name <- ifelse(!is.na(description) && description != "", description, channel)
  h <- hist(transformed_data, breaks = 200, plot = FALSE)
  plot(h,
       main = paste0(marker_name, " - ", sample_name, " ", title_suffix),
       xlab = "Intensité (transformée)", ylab = "Nombre d'événements",
       col = "lightblue", cex.main = 0.9, xaxt = "n")
  axis(1, at = axTicks(1), labels = round(axTicks(1), 2))
  axis(1, at = pretty(transformed_data, n = 20), labels = NA, tcl = -0.25)
  mtext(paste("a =", a, ", b =", b), side = 3, line = 0.2, cex = 0.6)
  if (!is.null(threshold)) {
    abline(v = asinh(a + b * threshold), col = "red", lwd = 2, lty = 2)
    legend("topright", legend = paste("Seuil:", round(threshold, 3)), col = "red", lwd = 2)
  }
}

show_fluorescence_histograms <- function(expr_data, channels, descriptions, transformation_list,
                                         filename = NULL, title_suffix = "", sample_name) {
  n_channels <- length(channels)
  n_cols <- 4
  n_rows <- ceiling(n_channels / n_cols)
  if (!is.null(filename)) {
    CairoPDF(filename, width = 15, height = 4 * n_rows)
  } else {
    create_plot_window(15, 4 * n_rows)
  }
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
  for (channel in channels) {
    a <- transformation_list[[channel]]$a
    b <- transformation_list[[channel]]$b
    marker_name <- get_marker_name(channel, descriptions)
    show_single_channel_histogram(expr_data, channel, marker_name, a, b, sample_name, title_suffix)
  }
  if (!is.null(filename)) {
    dev.off()
  } else {
    cat("Appuyez sur Entrée pour continuer...")
    readline()
    dev.off()
  }
}

create_dimred_plot <- function(df, x_col, y_col, color_col, sample_name, output_dir,
                               suffix = "", prefix = "UMAP", publication = FALSE,
                               separate_legend = FALSE) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  required <- c(x_col, y_col, color_col)
  if (!all(required %in% colnames(df))) {
    warning(paste("Colonnes manquantes :", paste(setdiff(required, colnames(df)), collapse=", ")))
    return(NULL)
  }
  custom_colors <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3",
    "#A65628", "#F781BF", "#00CED1", "#FFD700", "#1B9E77",
    "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
    "#006400", "#8B0000", "#00008B", "#FF4500", "#2E8B57",
    "#DC143C", "#00BFFF", "#32CD32", "#FF69B4", "#8A2BE2",
    "#FF8C00", "#00FA9A", "#4169E1", "#FF1493", "#20B2AA"
  )
  file_suffix <- if (suffix != "") paste0(prefix, suffix, ".pdf") else paste0(prefix, ".pdf")
  file <- file.path(output_dir, file_suffix)
  n_colors_needed <- length(unique(df[[color_col]]))
  if (n_colors_needed > length(custom_colors))
    custom_colors <- colorRampPalette(custom_colors)(n_colors_needed)
  
  # Graphique de base
  p <- ggplot(df, aes_string(x = x_col, y = y_col, color = color_col)) +
    geom_point(alpha = 0.4, size = if (publication) 0.5 else 0.6,
               stroke = if (publication) 0.1 else 0.1) +
    theme_minimal(base_size = if (publication) 11 else 12) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.ticks.length = unit(0.15, "cm"),
          legend.position = "right",
          legend.text = element_text(size = if (publication) 8 else 9),
          legend.title = element_text(size = if (publication) 9 else 10, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = if (publication) 12 else 14, face = "bold"),
          axis.title = element_text(size = if (publication) 10 else 11, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
    ggtitle(paste(prefix, "Visualization -", sample_name)) +
    xlab(paste(prefix, "Dimension 1")) + ylab(paste(prefix, "Dimension 2")) +
    scale_color_manual(values = custom_colors, name = "Clusters")
  
  if (separate_legend) {
    # --- PDF multipage : page 1 = graphique sans légende, page 2 = légende seule ---
    p_nolegend <- p + theme(legend.position = "none")
    
    # Extraction de la légende (nécessite ggpubr)
    if (!requireNamespace("ggpubr", quietly = TRUE)) {
      warning("ggpubr n'est pas disponible, la légende ne peut pas être séparée. Utilisation du mode normal.")
      separate_legend <- FALSE
    } else {
      legend_grob <- ggpubr::get_legend(p)
      
      # Page 2 : légende centrée
      p_legend <- ggplot() +
        theme_void() +
        annotation_custom(legend_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
      
      tryCatch({
        CairoPDF(file, width = 8.27, height = 11.69)  # A4 portrait
        print(p_nolegend)    # page 1
        grid::grid.newpage()
        print(p_legend)      # page 2
        dev.off()
        cat("✅ Graphique UMAP + légende séparée (A4, 2 pages) :", file, "\n")
      }, error = function(e) {
        cat("⚠️ Erreur création PDF multipage :", e$message, "\n")
        if (length(dev.list())) dev.off()
        # Fallback : graphique normal
        CairoPDF(file, width = if (publication) 10 else 12, height = if (publication) 8 else 8)
        print(p)
        dev.off()
      })
      return(file)
    }
  }
  
  # Mode normal (sans séparation) ou fallback
  tryCatch({
    CairoPDF(file, width = if (publication) 10 else 12, height = if (publication) 8 else 8)
    print(p)
    dev.off()
  }, error = function(e) {
    cat("⚠️ Erreur création PDF :", e$message, "\n")
    if (length(dev.list())) dev.off()
  })
  
  # Optionnel : fichier avec légende divisée en colonnes si trop de clusters
  if (n_colors_needed > 15 && !publication && !separate_legend) {
    file_split <- file.path(output_dir, gsub("\\.pdf$", "_SplitLegend.pdf", file_suffix))
    ncol_legend <- ceiling(n_colors_needed / 20)
    p2 <- p + theme(legend.position = "bottom",
                    legend.text = element_text(size = 7),
                    legend.key.size = unit(0.3, "cm")) +
      guides(color = guide_legend(ncol = ncol_legend, title.position = "top", title.hjust = 0.5))
    tryCatch({
      CairoPDF(file_split, width = 14, height = 10)
      print(p2)
      dev.off()
    }, error = function(e) if (length(dev.list())) dev.off())
  }
  return(file)
}

create_enhanced_umap_plot <- function(umap_df, cluster_labels, sample_name, output_dir, suffix = "",
                                      separate_legend = TRUE) {
  create_dimred_plot(umap_df, "UMAP1", "UMAP2", "Cluster", sample_name, output_dir,
                     suffix, "UMAP", publication = FALSE, separate_legend = separate_legend)
}

wrap_title <- function(title, max_chars = 20) {
  if (nchar(title) <= max_chars) return(title)
  words <- strsplit(title, " ")[[1]]
  line1 <- ""
  line2 <- ""
  for (w in words) {
    candidate <- if (nchar(line1) == 0) w else paste(line1, w)
    if (nchar(candidate) <= max_chars) {
      line1 <- candidate
    } else {
      line2 <- if (nchar(line2) == 0) w else paste(line2, w)
    }
  }
  if (nchar(line2) == 0) return(line1)
  paste0(line1, "\n", line2)
}


plot_fsc_ssc_by_cluster <- function(expr_data, cluster_ids, sample_name, output_dir) {
  if (!("FSC.A" %in% colnames(expr_data)) || !("SSC.A" %in% colnames(expr_data))) {
    warning("Les canaux FSC.A et/ou SSC.A ne sont pas présents dans les données.")
    return(NULL)
  }
  if (is.list(cluster_ids) || !is.vector(cluster_ids)) {
    cluster_ids <- as.numeric(unlist(cluster_ids))
  }
  pdf_file <- file.path(output_dir, "FSC_SSC_DotPlot_By_Cluster.pdf")
  CairoPDF(pdf_file, width = 10, height = 8)
  unique_clusters <- sort(unique(cluster_ids))
  n_clusters <- length(unique_clusters)
  par(mfrow = c(ceiling(sqrt(n_clusters)), ceiling(sqrt(n_clusters))),
      mar = c(4, 4, 3, 1))
  for (cluster in unique_clusters) {
    cluster_cells <- which(cluster_ids == cluster)
    other_cells <- which(cluster_ids != cluster)
    n_cells <- length(cluster_cells)
    plot(expr_data[other_cells, "FSC.A"], expr_data[other_cells, "SSC.A"],
         pch = 16, cex = 0.2, col = "gray",
         main = paste0("Cluster ", cluster, " (n=", n_cells, ")"),
         xlab = "FSC.A", ylab = "SSC.A",
         xlim = range(expr_data[, "FSC.A"]),
         ylim = range(expr_data[, "SSC.A"]))
    points(expr_data[cluster_cells, "FSC.A"], expr_data[cluster_cells, "SSC.A"],
           pch = 16, cex = 0.3, col = adjustcolor("red", alpha.f = 0.2))
  }
  dev.off()
  return(pdf_file)
}

plot_fsc_ssc_by_cluster_filtered <- function(expr_data, cluster_ids, cluster_labels, sample_name, output_dir) {
  if (!("FSC.A" %in% colnames(expr_data)) || !("SSC.A" %in% colnames(expr_data))) {
    warning("Les canaux FSC.A et/ou SSC.A ne sont pas présents dans les données.")
    return(NULL)
  }
  if (is.list(cluster_ids) || !is.vector(cluster_ids)) {
    cluster_ids <- as.numeric(unlist(cluster_ids))
  }
  pdf_file <- file.path(output_dir, "FSC_SSC_DotPlot_By_Cluster_Filtered.pdf")
  CairoPDF(pdf_file, width = 10, height = 8)
  cluster_id_to_label <- setNames(cluster_labels, 1:length(cluster_labels))
  unique_clusters <- sort(unique(cluster_ids))
  n_clusters <- length(unique_clusters)
  par(mfrow = c(ceiling(sqrt(n_clusters)), ceiling(sqrt(n_clusters))),
      mar = c(4, 4, 3, 1))
  for (cluster in unique_clusters) {
    cluster_cells <- which(cluster_ids == cluster)
    other_cells <- which(cluster_ids != cluster)
    n_cells <- length(cluster_cells)
    cluster_label <- cluster_id_to_label[as.character(cluster)]
    plot(expr_data[other_cells, "FSC.A"], expr_data[other_cells, "SSC.A"],
         pch = 16, cex = 0.2, col = "gray",
         main = wrap_title(paste0(cluster_label, " (n=", n_cells, ")")),
         xlab = "FSC.A", ylab = "SSC.A",
         xlim = range(expr_data[, "FSC.A"]),
         ylim = range(expr_data[, "SSC.A"]),
         cex.main = 0.8)
    points(expr_data[cluster_cells, "FSC.A"], expr_data[cluster_cells, "SSC.A"],
           pch = 16, cex = 0.3, col = adjustcolor("red", alpha.f = 0.2))
  }
  dev.off()
  return(pdf_file)
}

# ============================================================================
# SECTION 6: FONCTIONS DE GESTION INTERACTIVE DES CLUSTERS
# ============================================================================
rename_clusters_manually <- function(cluster_annotations, cluster_summary, sample_name) {
  cat("\n=== RENOMMAGE MANUEL DES CLUSTERS ===\n")
  cat("Vous allez renommer les clusters. Les visualisations (heatmap, UMAP) sont disponibles.\n")
  cat("Voici les annotations actuelles :\n\n")
  print(cluster_summary[, c("Cluster", "Annotation", "Count", "Percentage")])
  new_annotations <- cluster_annotations
  for (i in 1:length(cluster_annotations)) {
    cat(paste("\nCluster", i, ":", cluster_annotations[i]))
    cat(paste(" (", cluster_summary$Count[i], "cellules,", cluster_summary$Percentage[i], "%)"))
    cat("\nNouveau nom (Entrée pour garder '", cluster_annotations[i], "') : ")
    new_name <- readline()
    if (new_name != "") {
      new_annotations[i] <- new_name
      cat(paste("Cluster", i, "renommé :", new_annotations[i], "\n"))
    }
  }
  cat("\nRésumé des modifications :\n")
  comparison <- data.frame(
    Cluster = 1:length(cluster_annotations),
    Ancien = cluster_annotations,
    Nouveau = new_annotations
  )
  print(comparison)
  cat("\nConfirmez ces modifications? (o/n): ")
  confirm <- readline()
  if (tolower(confirm) %in% c("o", "oui", "y", "yes")) {
    return(new_annotations)
  } else {
    cat("Modifications annulées, conservation des noms précédents.\n")
    return(cluster_annotations)
  }
}

filter_clusters_interactive <- function(expr_sub, cluster_ids, cluster_annotations, markers_of_interest,
                                        marker_descriptions, sample_name, output_dir, umap_df, cluster_labels, annotation_results) {
  cat("\n=== FILTRAGE INTERACTIF DES CLUSTERS ===\n")
  unique_clusters <- sort(unique(cluster_ids))
  n_clusters_before <- length(unique_clusters)
  cluster_counts <- as.numeric(table(cluster_ids))
  total_cells_before <- sum(cluster_counts)
  cluster_stats <- data.frame(
    Cluster = unique_clusters,
    Annotation = cluster_annotations[unique_clusters],
    Count = cluster_counts,
    Percentage = round(100 * cluster_counts / total_cells_before, 2)
  )
  cat("\nClusters actuels:\n")
  cat("-----------------\n")
  print(cluster_stats)
  cat("\nEntrez les numéros des clusters à éliminer, séparés par des virgules (ex: 1,2,3) : ")
  clusters_input <- readline()
  clusters_to_remove <- as.numeric(unlist(strsplit(clusters_input, ",")))
  clusters_to_remove <- unique(clusters_to_remove[!is.na(clusters_to_remove)])
  clusters_to_remove <- intersect(clusters_to_remove, cluster_stats$Cluster)
  if (length(clusters_to_remove) == 0) {
    cat("❌ Aucun cluster valide sélectionné pour élimination. Aucun filtre effectué.\n")
    return(list(expr_sub = expr_sub, cluster_ids = cluster_ids, cluster_annotations = cluster_annotations,
                cluster_labels = cluster_labels, umap_df = umap_df, annotation_results = annotation_results,
                filter_info = NULL))
  }
  cat("\nClusters sélectionnés pour élimination:\n")
  total_cells_to_remove <- 0
  for (cl in clusters_to_remove) {
    idx <- which(cluster_stats$Cluster == cl)
    cell_count <- cluster_stats$Count[idx]
    total_cells_to_remove <- total_cells_to_remove + cell_count
    cat(sprintf(" Cluster %d: %s (%d cellules, %.2f%%)\n",
                cl, cluster_stats$Annotation[idx], cell_count, cluster_stats$Percentage[idx]))
  }
  cat(sprintf("\nTotal des cellules à éliminer: %d (%.2f%% du total)\n",
              total_cells_to_remove, 100 * total_cells_to_remove / total_cells_before))
  cat("\nConfirmez-vous l'élimination de ces clusters? (y/n): ")
  confirm_remove <- readline()
  if (tolower(confirm_remove) %in% c("y", "yes")) {
    keep_mask <- !(cluster_ids %in% clusters_to_remove)
    expr_sub_filtered <- expr_sub[keep_mask, ]
    cluster_ids_filtered <- cluster_ids[keep_mask]
    old_ids <- sort(unique(cluster_ids_filtered))
    mapping <- setNames(1:length(old_ids), old_ids)
    cluster_ids_renumbered <- sapply(cluster_ids_filtered, function(x) mapping[as.character(x)])
    new_annotations <- cluster_annotations[old_ids]
    new_labels <- paste0(1:length(new_annotations), "_", new_annotations)
    if (!is.null(annotation_results)) {
      new_scores <- annotation_results$scores[old_ids]
      new_details <- annotation_results$details[old_ids]
      annotation_results$annotations <- new_annotations
      annotation_results$scores <- new_scores
      annotation_results$details <- new_details
    }
    umap_df_filtered <- umap_df[keep_mask, ]
    umap_df_filtered$Cluster <- factor(new_labels[cluster_ids_renumbered], levels = new_labels)
    total_cells_after <- nrow(expr_sub_filtered)
    n_clusters_after <- length(unique(cluster_ids_renumbered))
    cat("\n🔍 VÉRIFICATION DU FILTRAGE :\n")
    cat("-----------------------------\n")
    cat(sprintf("Cellules avant filtrage: %d\n", total_cells_before))
    cat(sprintf("Cellules après filtrage: %d\n", total_cells_after))
    cat(sprintf("Clusters avant filtrage: %d\n", n_clusters_before))
    cat(sprintf("Clusters après filtrage: %d\n", n_clusters_after))
    cat(sprintf("Cellules éliminées: %d (%.2f%%)\n",
                total_cells_before - total_cells_after,
                100 * (total_cells_before - total_cells_after) / total_cells_before))
    filter_info <- list(
      removed_clusters = clusters_to_remove,
      removed_annotations = cluster_annotations[clusters_to_remove],
      total_cells_removed = total_cells_to_remove,
      total_cells_before = total_cells_before,
      total_cells_after = total_cells_after,
      n_clusters_before = n_clusters_before,
      n_clusters_after = n_clusters_after
    )
    filter_df <- data.frame(
      Removed_Clusters = paste(clusters_to_remove, collapse = ", "),
      Removed_Annotations = paste(cluster_annotations[clusters_to_remove], collapse = ", "),
      Total_Cells_Removed = total_cells_to_remove,
      Total_Cells_Before = total_cells_before,
      Total_Cells_After = total_cells_after,
      Clusters_Before = n_clusters_before,
      Clusters_After = n_clusters_after
    )
    filter_file <- file.path(output_dir, "Cluster_Filter_History.csv")
    if (file.exists(filter_file)) {
      existing_df <- read.csv(filter_file)
      filter_df <- rbind(existing_df, filter_df)
    }
    write.csv(filter_df, filter_file, row.names = FALSE)
    cat(sprintf("\n📋 Historique des filtrages sauvegardé : %s\n", filter_file))
    cat("\n✅ Filtrage terminé avec succès!\n")
    graphics.off()
    fsc_ssc_filtered_file <- plot_fsc_ssc_by_cluster_filtered(
      expr_sub_filtered, cluster_ids_renumbered, new_labels, sample_name, output_dir
    )
    return(list(expr_sub = expr_sub_filtered, cluster_ids = cluster_ids_renumbered,
                cluster_annotations = new_annotations, cluster_labels = new_labels,
                umap_df = umap_df_filtered, annotation_results = annotation_results,
                filter_info = filter_info, fsc_ssc_filtered_file = fsc_ssc_filtered_file))
  } else {
    cat("Filtrage annulé.\n")
    return(list(expr_sub = expr_sub, cluster_ids = cluster_ids,
                cluster_annotations = cluster_annotations, cluster_labels = cluster_labels,
                umap_df = umap_df, annotation_results = annotation_results,
                filter_info = NULL, fsc_ssc_filtered_file = NULL))
  }
}

merge_clusters_interactive <- function(expr_sub, cluster_ids, cluster_annotations, markers_of_interest,
                                       marker_descriptions, sample_name, output_dir, umap_df, cluster_labels,
                                       annotation_results, suffix = "") {
  cat("\n=== FUSION INTERACTIVE DES CLUSTERS ===\n")
  unique_clusters <- sort(unique(cluster_ids))
  n_clusters_before <- length(unique_clusters)
  cluster_counts <- as.numeric(table(cluster_ids))
  total_cells_before <- sum(cluster_counts)
  cluster_stats <- data.frame(
    Cluster = unique_clusters,
    Annotation = cluster_annotations[unique_clusters],
    Count = cluster_counts,
    Percentage = round(100 * cluster_counts / total_cells_before, 2)
  )
  cat("\nClusters actuels:\n")
  cat("-----------------\n")
  print(cluster_stats)
  cat("\nEntrez les numéros des clusters à fusionner, séparés par des virgules (ex: 1,2,3): ")
  clusters_input <- readline()
  clusters_to_merge <- as.numeric(unlist(strsplit(clusters_input, ",")))
  clusters_to_merge <- unique(clusters_to_merge[!is.na(clusters_to_merge)])
  clusters_to_merge <- intersect(clusters_to_merge, cluster_stats$Cluster)
  if (length(clusters_to_merge) < 2) {
    cat("❌ Vous devez sélectionner au moins 2 clusters pour fusionner. Aucune fusion effectuée.\n")
    return(list(cluster_ids = cluster_ids, cluster_annotations = cluster_annotations,
                cluster_labels = cluster_labels, umap_df = umap_df,
                annotation_results = annotation_results, merge_info = NULL,
                fsc_ssc_filtered_file = NULL))
  }
  cat("\nClusters sélectionnés pour fusion:\n")
  total_cells_to_merge <- 0
  for (cl in clusters_to_merge) {
    idx <- which(cluster_stats$Cluster == cl)
    cell_count <- cluster_stats$Count[idx]
    total_cells_to_merge <- total_cells_to_merge + cell_count
    cat(sprintf(" Cluster %d: %s (%d cellules, %.2f%%)\n",
                cl, cluster_stats$Annotation[idx], cell_count, cluster_stats$Percentage[idx]))
  }
  cat("\nEntrez le nom pour le nouveau cluster fusionné (ex: 'cDC1_merged', 'T_cells_combined'): ")
  new_cluster_name <- readline()
  if (new_cluster_name == "") {
    original_annotations <- cluster_stats$Annotation[cluster_stats$Cluster %in% clusters_to_merge]
    new_cluster_name <- paste(unique(original_annotations), collapse = "_")
    cat(sprintf("Nom par défaut utilisé: %s\n", new_cluster_name))
  }
  cat(sprintf("\nConfirmez-vous la fusion des clusters %s en un nouveau cluster '%s'? (y/n): ",
              paste(clusters_to_merge, collapse = ", "), new_cluster_name))
  confirm_merge <- readline()
  if (tolower(confirm_merge) %in% c("y", "yes")) {
    original_cluster_ids <- cluster_ids
    new_cluster_id <- min(clusters_to_merge)
    cluster_ids_updated <- cluster_ids
    for (cl in clusters_to_merge) {
      cluster_ids_updated[cluster_ids == cl] <- new_cluster_id
    }
    old_ids <- sort(unique(cluster_ids_updated))
    mapping <- setNames(1:length(old_ids), old_ids)
    cluster_ids_renumbered <- sapply(cluster_ids_updated, function(x) mapping[as.character(x)])
    new_annotations <- character(length(old_ids))
    new_labels <- character(length(old_ids))
    for (i in 1:length(old_ids)) {
      old_id <- old_ids[i]
      if (old_id == new_cluster_id) {
        new_annotations[i] <- new_cluster_name
        new_labels[i] <- paste0(i, "_", new_cluster_name)
      } else {
        orig_idx <- which(unique_clusters == old_id)
        new_annotations[i] <- cluster_annotations[orig_idx]
        new_labels[i] <- paste0(i, "_", new_annotations[i])
      }
    }
    if (!is.null(annotation_results)) {
      new_scores <- numeric(length(old_ids))
      new_details <- list()
      for (i in 1:length(old_ids)) {
        old_id <- old_ids[i]
        if (old_id == new_cluster_id) {
          scores_to_merge <- annotation_results$scores[clusters_to_merge]
          counts_to_merge <- cluster_counts[clusters_to_merge]
          weighted_score <- sum(scores_to_merge * counts_to_merge) / sum(counts_to_merge)
          new_scores[i] <- weighted_score
          new_details[[i]] <- list(
            annotation = new_cluster_name,
            score = weighted_score,
            details = paste("Fusion des clusters", paste(clusters_to_merge, collapse = ", ")),
            confidence_threshold = NA
          )
        } else {
          orig_idx <- which(unique_clusters == old_id)
          new_scores[i] <- annotation_results$scores[orig_idx]
          new_details[[i]] <- annotation_results$details[[orig_idx]]
        }
      }
      annotation_results$annotations <- new_annotations
      annotation_results$scores <- new_scores
      annotation_results$details <- new_details
    }
    umap_df_updated <- umap_df
    umap_df_updated$Cluster <- factor(new_labels[cluster_ids_renumbered], levels = new_labels)
    total_cells_after <- length(cluster_ids_renumbered)
    n_clusters_after <- length(unique(cluster_ids_renumbered))
    expected_clusters_after <- n_clusters_before - (length(clusters_to_merge) - 1)
    cat("\n🔍 VÉRIFICATION DE LA FUSION :\n")
    cat("-----------------------------\n")
    cat(sprintf("Clusters avant fusion: %d\n", n_clusters_before))
    cat(sprintf("Clusters après fusion: %d\n", n_clusters_after))
    if (n_clusters_after == expected_clusters_after) {
      cat("✅ Le nombre de clusters est correct\n")
    } else {
      cat(sprintf("⚠️ Attention: nombre de clusters différent de l'attendu (attendu: %d)\n", expected_clusters_after))
    }
    merge_info <- list(
      merged_clusters = clusters_to_merge,
      new_cluster_name = new_cluster_name,
      original_annotations = cluster_annotations[clusters_to_merge],
      total_cells_to_merge = total_cells_to_merge,
      total_cells_before = total_cells_before,
      total_cells_after = total_cells_after,
      n_clusters_before = n_clusters_before,
      n_clusters_after = n_clusters_after
    )
    cat(sprintf("\n✅ Fusion terminée avec succès!\n"))
    cat(sprintf(" Nouveau cluster: %s (%d cellules, %.2f%%)\n",
                new_cluster_name, total_cells_to_merge,
                round(100 * total_cells_to_merge / total_cells_after, 2)))
    merge_df <- data.frame(
      Merged_Clusters = paste(clusters_to_merge, collapse = ", "),
      New_Cluster_Name = new_cluster_name,
      Original_Annotations = paste(cluster_annotations[clusters_to_merge], collapse = ", "),
      Cells_in_New_Cluster = total_cells_to_merge,
      Total_Cells = total_cells_after,
      Clusters_Before = n_clusters_before,
      Clusters_After = n_clusters_after
    )
    if (suffix != "") {
      merge_file <- file.path(output_dir, paste0("Cluster_Merge_History", suffix, ".csv"))
    } else {
      merge_file <- file.path(output_dir, "Cluster_Merge_History.csv")
    }
    write.csv(merge_df, merge_file, row.names = FALSE)
    cat(sprintf("\n📋 Historique des fusions sauvegardé : %s\n", merge_file))
    fsc_ssc_filtered_file <- plot_fsc_ssc_by_cluster_filtered(
      expr_sub, cluster_ids_renumbered, new_labels, sample_name, output_dir
    )
    return(list(cluster_ids = cluster_ids_renumbered,
                cluster_annotations = new_annotations,
                cluster_labels = new_labels,
                umap_df = umap_df_updated,
                annotation_results = annotation_results,
                merge_info = merge_info,
                fsc_ssc_filtered_file = fsc_ssc_filtered_file))
  } else {
    cat("Fusion annulée.\n")
    return(list(cluster_ids = cluster_ids, cluster_annotations = cluster_annotations,
                cluster_labels = cluster_labels, umap_df = umap_df,
                annotation_results = annotation_results, merge_info = NULL,
                fsc_ssc_filtered_file = NULL))
  }
}

# ============================================================================
# SECTION 7: SAUVEGARDE DES MÉDIANES
# ============================================================================
save_cluster_medians_with_thresholds <- function(expr_sub, cluster_ids, marker_thresholds,
                                                 marker_descriptions, cluster_labels,
                                                 output_dir, sample_name) {
  cluster_medians <- aggregate(expr_sub, by = list(Cluster = cluster_ids), FUN = median)
  rownames(cluster_medians) <- cluster_medians$Cluster
  cluster_medians$Cluster <- NULL
  colnames(cluster_medians) <- sapply(colnames(cluster_medians),
                                      function(ch) get_marker_name(ch, marker_descriptions))
  medians_file <- file.path(output_dir, "cluster_medians.csv")
  write.csv(cluster_medians, medians_file, row.names = TRUE)
  return(list(medians = cluster_medians, thresholds = NULL))
}

# ============================================================================
# SECTION 8: SOUS-ÉCHANTILLONNAGE
# ============================================================================
perform_global_subsampling <- function(expr_data, max_cells = 200000) {
  total_cells <- nrow(expr_data)
  if (total_cells > max_cells) {
    set.seed(123)
    idx <- sample(1:total_cells, max_cells)
    expr_sub <- expr_data[idx, , drop = FALSE]
    cat(sprintf("Sous‑échantillonnage global : %d cellules sélectionnées sur %d.\n", max_cells, total_cells))
  } else {
    expr_sub <- expr_data
    cat(sprintf("Pas de sous‑échantillonnage : %d cellules disponibles.\n", total_cells))
  }
  return(expr_sub)
}

perform_balanced_subsampling <- function(expr_data, conditions, mode = "condition", n_cells_target = NULL) {
  if (!"FolderID" %in% colnames(expr_data) || !"FileID" %in% colnames(expr_data)) {
    warning("Colonnes FolderID ou FileID non trouvées. Retour des données originales.")
    return(expr_data)
  }
  
  cat("\n=== SOUS-ÉCHANTILLONNAGE ÉQUILIBRÉ ===\n")
  
  if (mode == "condition") {
    unique_conditions <- sort(unique(expr_data[, "FolderID"]))
    if (is.null(n_cells_target)) {
      cond_counts <- sapply(unique_conditions, function(cid) sum(expr_data[, "FolderID"] == cid))
      n_cells_target <- min(cond_counts)
      cat(sprintf("Taille cible par condition = %d (minimum disponible)\n", n_cells_target))
    }
    results <- list()
    for (cond_id in unique_conditions) {
      cond_name <- get_condition_name(cond_id, conditions)
      cond_data <- expr_data[expr_data[, "FolderID"] == cond_id, ]
      n_avail <- nrow(cond_data)
      if (n_avail >= n_cells_target) {
        set.seed(123 + cond_id)
        idx <- sample(1:n_avail, n_cells_target)
        cond_sampled <- cond_data[idx, ]
      } else {
        cat(sprintf("⚠️ Condition %s : seulement %d cellules disponibles (cible=%d). Toutes conservées.\n",
                    cond_name, n_avail, n_cells_target))
        cond_sampled <- cond_data
      }
      results[[as.character(cond_id)]] <- cond_sampled
      cat(sprintf("  %s : %d cellules échantillonnées (sur %d)\n", cond_name, nrow(cond_sampled), n_avail))
    }
    final_sampled <- do.call(rbind, results)
    
  } else if (mode == "file") {
    unique_files <- sort(unique(expr_data[, "FileID"]))
    if (is.null(n_cells_target)) {
      file_counts <- sapply(unique_files, function(fid) sum(expr_data[, "FileID"] == fid))
      n_cells_target <- min(file_counts)
      cat(sprintf("Taille cible par fichier = %d (minimum disponible)\n", n_cells_target))
    }
    results <- list()
    for (file_id in unique_files) {
      file_data <- expr_data[expr_data[, "FileID"] == file_id, ]
      n_avail <- nrow(file_data)
      file_name <- paste0("File_", file_id)
      if (n_avail >= n_cells_target) {
        set.seed(123 + file_id)
        idx <- sample(1:n_avail, n_cells_target)
        file_sampled <- file_data[idx, ]
      } else {
        cat(sprintf("⚠️ Fichier %s : seulement %d cellules disponibles (cible=%d). Toutes conservées.\n",
                    file_name, n_avail, n_cells_target))
        file_sampled <- file_data
      }
      results[[as.character(file_id)]] <- file_sampled
      cat(sprintf("  %s : %d cellules échantillonnées (sur %d)\n", file_name, nrow(file_sampled), n_avail))
    }
    final_sampled <- do.call(rbind, results)
  } else {
    stop("Mode de sous-échantillonnage inconnu.")
  }
  
  cat(sprintf("\nTotal final après sous-échantillonnage : %d cellules\n", nrow(final_sampled)))
  return(final_sampled)
}

# ============================================================================
# SECTION 9: FONCTIONS POUR ANALYSE PAR CONDITION
# ============================================================================
detect_conditions_from_metadata <- function(desc) {
  conditions <- list()
  folder_keys <- names(desc)[grep("^Folder[0-9]+_Name$", names(desc))]
  if (length(folder_keys) > 0) {
    for (key in folder_keys) {
      folder_id <- as.numeric(gsub("Folder([0-9]+)_Name", "\\1", key))
      folder_name <- desc[[key]]
      conditions[[as.character(folder_id)]] <- folder_name
    }
    cat("✅ Conditions détectées depuis les métadonnées :\n")
    for (id in sort(as.numeric(names(conditions)))) {
      cat(sprintf("  FolderID %s : %s\n", id, conditions[[as.character(id)]]))
    }
  } else {
    conditions[["1"]] <- "Condition_1"
    conditions[["2"]] <- "Condition_2"
    cat("⚠️  Aucune condition détectée dans les métadonnées. Utilisation des valeurs par défaut.\n")
  }
  return(conditions)
}

get_condition_name <- function(folder_id, conditions) {
  id_char <- as.character(folder_id)
  if (id_char %in% names(conditions)) {
    return(conditions[[id_char]])
  } else {
    return(paste("Condition", folder_id))
  }
}

create_proportion_heatmap_by_file <- function(prop_data, sample_name, output_dir, conditions) {
  prop_table <- table(prop_data$Cluster, prop_data$File)
  prop_percent <- prop.table(prop_table, margin = 2) * 100
  
  prop_long <- as.data.frame(prop_table)
  colnames(prop_long) <- c("Cluster", "File", "Count")
  prop_long$Percentage <- as.vector(prop_percent)
  
  file_condition <- unique(prop_data[, c("File", "Condition")])
  rownames(file_condition) <- file_condition$File
  
  file_order_df <- unique(prop_data[, c("File", "Condition")])
  file_order_df <- file_order_df[order(file_order_df$Condition, natural_sort(as.character(file_order_df$File))), ]
  file_levels <- file_order_df$File
  prop_long$File <- factor(prop_long$File, levels = file_levels)
  
  n_conditions <- length(unique(file_condition$Condition))
  condition_colors <- scales::hue_pal()(n_conditions)
  names(condition_colors) <- unique(file_condition$Condition)
  
  heatmap_file <- file.path(output_dir, "Proportion_Heatmap_by_File.pdf")
  
  CairoPDF(heatmap_file, width = max(10, length(file_levels) * 0.8),
           height = max(8, length(unique(prop_long$Cluster)) * 0.4))
  
  p <- ggplot(prop_long, aes(x = File, y = Cluster, fill = Percentage)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(Percentage >= 0.01, sprintf("%.1f%%\n(n=%d)", Percentage, Count), "")),
              size = 2.5, color = "orange", check_overlap = FALSE) +
    scale_fill_gradientn(colors = viridis(256), name = "Pourcentage (%)") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right",
      legend.key.height = unit(1.5, "cm"),
      panel.grid = element_blank()
    ) +
    ggtitle(paste("Proportions des Clusters par Fichier -", sample_name)) +
    xlab("Fichier") + ylab("Cluster")
  
  y_top <- length(unique(prop_long$Cluster)) + 0.5
  y_bottom <- length(unique(prop_long$Cluster)) + 1.0
  y_text   <- y_bottom + 0.4
  
  for (i in seq_along(file_levels)) {
    f <- file_levels[i]
    cond <- file_condition[f, "Condition"]
    p <- p + annotate("rect",
                      xmin = i - 0.4, xmax = i + 0.4,
                      ymin = y_top, ymax = y_bottom,
                      fill = condition_colors[cond], alpha = 0.7)
  }
  p <- p + annotate("text",
                    x = 1:length(file_levels),
                    y = y_text,
                    label = sapply(file_levels, function(f) file_condition[f, "Condition"]),
                    size = 2, fontface = "bold", angle = 0, hjust = 0.5)
  
  p <- p + coord_cartesian(ylim = c(0.5, y_bottom + 0.5), clip = "off") +
    theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10))
  
  print(p)
  dev.off()
  cat("✅ Heatmap des proportions par fichier générée :", heatmap_file, "\n")
}

create_stacked_bar_charts <- function(prop_long, output_dir, sample_name) {
  stacked_file <- file.path(output_dir, "Stacked_Bar_Charts_by_Condition.pdf")
  
  custom_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3",
                     "#A65628", "#F781BF", "#00CED1", "#FFD700", "#1B9E77",
                     "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
                     "#006400", "#8B0000", "#00008B", "#FF4500", "#2E8B57",
                     "#DC143C", "#00BFFF", "#32CD32", "#FF69B4", "#8A2BE2",
                     "#FF8C00", "#00FA9A", "#4169E1", "#FF1493", "#20B2AA")
  n_clusters <- length(unique(prop_long$Cluster))
  if (n_clusters > length(custom_colors))
    custom_colors <- colorRampPalette(custom_colors)(n_clusters)
  
  # Graphique 1 : pourcentages par condition (sans légende)
  p1 <- ggplot(prop_long, aes(x = Condition, y = Percentage_Condition, fill = Cluster)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", linewidth = 0.3) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
    ggtitle(paste("Répartition des populations par condition -", sample_name)) +
    xlab("Condition") + ylab("Pourcentage de la condition (%)") +
    scale_fill_manual(values = custom_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    geom_text(aes(label = ifelse(Percentage_Condition > 5, sprintf("%.1f%%", Percentage_Condition), "")),
              position = position_stack(vjust = 0.5), size = 3.5, color = "white", fontface = "bold")
  
  # Graphique 2 : comptages absolus par condition (sans légende)
  p2 <- ggplot(prop_long, aes(x = Condition, y = Count, fill = Cluster)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", linewidth = 0.3) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
    ggtitle(paste("Comptages absolus par condition -", sample_name)) +
    xlab("Condition") + ylab("Nombre de cellules") +
    scale_fill_manual(values = custom_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    geom_text(aes(label = ifelse(Count > max(Count) * 0.05, format(Count, big.mark = " "), "")),
              position = position_stack(vjust = 0.5), size = 3.2, color = "white", fontface = "bold")
  
  # Combinaison des deux graphiques sur une même page (côte à côte)
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    warning("ggpubr n'est pas disponible, la légende ne peut pas être séparée. Utilisation du mode original.")
    CairoPDF(stacked_file, width = 10, height = 8)
    print(p1 + theme(legend.position = "right"))
    print(p2 + theme(legend.position = "right"))
    dev.off()
    return()
  }
  
  combined <- ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = FALSE)
  
  # Création d'un graphique temporaire pour extraire la légende
  p_temp <- ggplot(prop_long, aes(x = Condition, y = Percentage_Condition, fill = Cluster)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = custom_colors, name = "Clusters")
  legend_grob <- ggpubr::get_legend(p_temp)
  
  # Page 2 : légende seule centrée
  p_legend <- ggplot() +
    theme_void() +
    annotation_custom(legend_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
  
  # PDF multipage format A4 paysage (11.69" × 8.27")
  tryCatch({
    CairoPDF(stacked_file, width = 11.69, height = 8.27)
    print(combined)      # page 1
    grid::grid.newpage()
    print(p_legend)      # page 2
    dev.off()
    cat("✅ Graphiques stacked bar + légende séparée (A4 paysage, 2 pages) :", stacked_file, "\n")
  }, error = function(e) {
    cat("⚠️ Erreur création PDF multipage :", e$message, "\n")
    if (length(dev.list())) dev.off()
    CairoPDF(stacked_file, width = 10, height = 8)
    print(p1 + theme(legend.position = "right"))
    print(p2 + theme(legend.position = "right"))
    dev.off()
  })
}

create_condition_dimred_plot <- function(df, x_col, y_col, sample_name, output_dir, prefix, conditions) {
  if (!"Condition" %in% colnames(df)) return(NULL)
  
  n_cond <- length(unique(df$Condition))
  color_palette <- scales::hue_pal()(n_cond)
  
  p <- ggplot(df, aes_string(x = x_col, y = y_col, color = "Condition")) +
    geom_point(alpha = 0.4, size = 0.4) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
    ggtitle(paste(prefix, "par Condition -", sample_name)) +
    xlab(paste(prefix, "Dimension 1")) + ylab(paste(prefix, "Dimension 2")) +
    scale_color_manual(values = color_palette)
  
  file <- file.path(output_dir, paste0(prefix, "_by_Condition.pdf"))
  CairoPDF(file, width = 10, height = 8)
  print(p)
  dev.off()
  cat(paste0("✅ ", prefix, " par condition généré :", file, "\n"))
}

create_facet_dimred_by_file <- function(dimred_df, prop_data, sample_name, output_dir, prefix, conditions) {
  if (!("File" %in% colnames(prop_data))) return(NULL)
  
  df <- dimred_df
  df$File <- prop_data$File
  df$Condition <- prop_data$Condition
  
  file_order_df <- unique(prop_data[, c("File", "Condition")])
  file_order_df <- file_order_df[order(file_order_df$Condition, natural_sort(as.character(file_order_df$File))), ]
  file_levels <- file_order_df$File
  df$File <- factor(df$File, levels = file_levels)
  
  file_counts <- table(df$File)
  file_labels <- paste0(file_levels, "\n(n=", file_counts[file_levels], ")")
  names(file_labels) <- file_levels
  df$FileLabel <- factor(df$File, levels = file_levels, labels = file_labels)
  
  p1 <- ggplot(df, aes_string(x = colnames(dimred_df)[1], y = colnames(dimred_df)[2])) +
    geom_point(alpha = 0.2, size = 0.3, color = "steelblue") +
    facet_wrap(~ FileLabel, ncol = 4) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "lightgray", color = "black"),
          strip.text = element_text(size = 8, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
    ggtitle(paste(prefix, "Facets par Fichier -", sample_name)) +
    xlab(paste(prefix, "Dimension 1")) + ylab(paste(prefix, "Dimension 2"))
  
  file1 <- file.path(output_dir, paste0(prefix, "_Facets_by_File.pdf"))
  CairoPDF(file1, width = 14, height = 10)
  print(p1)
  dev.off()
  cat(paste0("✅ ", prefix, " en facets par fichier généré :", file1, "\n"))
  
  n_cond <- length(unique(df$Condition))
  color_palette <- scales::hue_pal()(n_cond)
  
  p2 <- ggplot(df, aes_string(x = colnames(dimred_df)[1], y = colnames(dimred_df)[2], color = "Condition")) +
    geom_point(alpha = 0.2, size = 0.3) +
    facet_wrap(~ FileLabel, ncol = 4) +
    scale_color_manual(values = color_palette) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "lightgray", color = "black"),
          strip.text = element_text(size = 8, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "bottom") +
    ggtitle(paste(prefix, "Facets par Fichier (couleur par condition) -", sample_name)) +
    xlab(paste(prefix, "Dimension 1")) + ylab(paste(prefix, "Dimension 2"))
  
  file2 <- file.path(output_dir, paste0(prefix, "_Facets_by_File_with_Condition.pdf"))
  CairoPDF(file2, width = 14, height = 10)
  print(p2)
  dev.off()
  cat(paste0("✅ ", prefix, " en facets avec condition généré :", file2, "\n"))
}

perform_condition_analysis <- function(expr_sub, cluster_ids, cluster_labels,
                                       sample_name, output_dir, desc = NULL,
                                       umap_df = NULL) {
  
  if (!("FolderID" %in% colnames(expr_sub)) || !("FileID" %in% colnames(expr_sub))) {
    cat("[WARNING] Colonnes 'FolderID' ou 'FileID' absentes. Analyse par condition ignoree.\n")
    return(NULL)
  }
  
  cat("\n[INFO] Analyse par condition (donnees concatenees) ...\n")
  
  conditions <- detect_conditions_from_metadata(desc)
  
  file_mapping <- list()
  if (!is.null(desc)) {
    file_keys <- names(desc)[grep("^File[0-9]+_Name$", names(desc))]
    if (length(file_keys) > 0) {
      for (key in file_keys) {
        file_id <- as.numeric(gsub("File([0-9]+)_Name", "\\1", key))
        file_name <- desc[[key]]
        file_name <- gsub("\\.fcs$", "", file_name, ignore.case = TRUE)
        file_mapping[[as.character(file_id)]] <- file_name
      }
    }
  }
  
  folder_vec <- expr_sub[, "FolderID"]
  file_vec   <- expr_sub[, "FileID"]
  
  file_display <- sapply(file_vec, function(id) {
    id_char <- as.character(id)
    if (!is.na(id) && id_char %in% names(file_mapping)) {
      return(file_mapping[[id_char]])
    } else {
      return(paste0("File_", id_char))
    }
  })
  
  condition_vec <- sapply(folder_vec, function(id) get_condition_name(id, conditions))
  
  prop_data <- data.frame(
    Cluster   = cluster_labels[cluster_ids],
    FileID    = file_vec,
    File      = file_display,
    FolderID  = folder_vec,
    Condition = condition_vec,
    stringsAsFactors = FALSE
  )
  
  prop_table_cond <- table(prop_data$Cluster, prop_data$Condition)
  prop_percent_condition <- prop.table(prop_table_cond, margin = 2) * 100
  
  counts_csv <- file.path(output_dir, "Cluster_Counts_by_Condition.csv")
  write.csv(as.data.frame.matrix(prop_table_cond), file = counts_csv)
  write.csv(as.data.frame.matrix(round(prop_percent_condition, 2)),
            file = file.path(output_dir, "Cluster_Percent_by_Condition.csv"))
  
  prop_long_cond <- as.data.frame(prop_table_cond)
  colnames(prop_long_cond) <- c("Cluster", "Condition", "Count")
  prop_long_cond$Percentage_Condition <- as.vector(prop_percent_condition)
  prop_long_cond$Percentage_Total <- as.vector(prop.table(prop_table_cond) * 100)
  write.csv(prop_long_cond, file.path(output_dir, "Cluster_Quantification_by_Condition_Long.csv"), row.names = FALSE)
  
  unique_conditions <- unique(prop_data$Condition)
  
  counts_per_sample <- table(prop_data$Cluster, prop_data$File)
  props_per_sample <- prop.table(counts_per_sample, margin = 2)
  props_transformed <- asin(sqrt(props_per_sample))
  
  sample_conditions <- unique(prop_data[, c("File", "Condition")])
  stopifnot(all(colnames(props_transformed) %in% sample_conditions$File))
  cond_vector <- sample_conditions$Condition[
    match(colnames(props_transformed), sample_conditions$File)
  ]
  
  prop_sample_file <- file.path(output_dir, "Proportions_Par_Echantillon.csv")
  write.csv(as.data.frame.matrix(round(props_per_sample * 100, 2)), file = prop_sample_file)
  cat("[INFO] Proportions par echantillon sauvegardees dans", prop_sample_file, "\n")
  
  if (length(unique_conditions) == 2) {
    cat("[TEST] Analyse statistique : comparaison directe entre deux conditions (pas de design temporel detecte).\n")
    
    cond_levels <- unique(cond_vector)
    idx_cond1 <- which(cond_vector == cond_levels[1])
    idx_cond2 <- which(cond_vector == cond_levels[2])
    
    n_clusters <- nrow(props_transformed)
    p_values <- numeric(n_clusters)
    test_used <- character(n_clusters)
    mean_diff <- numeric(n_clusters)
    
    for (cl in seq_len(n_clusters)) {
      vals1 <- props_transformed[cl, idx_cond1]
      vals2 <- props_transformed[cl, idx_cond2]
      
      if (length(vals1) >= 2 && length(vals2) >= 2) {
        res <- t.test(vals2, vals1, var.equal = FALSE)
        p_values[cl] <- res$p.value
        test_used[cl] <- "Welch t-test"
        mean_diff[cl] <- mean(vals2) - mean(vals1)
      } else {
        p_values[cl] <- NA
        test_used[cl] <- "Effectif insuffisant (<2)"
        mean_diff[cl] <- NA
      }
    }
    
    adj_p_values <- p.adjust(p_values, method = "BH")
    significance <- ifelse(adj_p_values < 0.05 & !is.na(adj_p_values), "Oui", "Non")
    
    results_df <- data.frame(
      Cluster = rownames(props_transformed),
      Annotation = rownames(props_transformed),
      Moyenne_Transformee_Cond1 = rowMeans(props_transformed[, idx_cond1, drop = FALSE]),
      Moyenne_Transformee_Cond2 = rowMeans(props_transformed[, idx_cond2, drop = FALSE]),
      Difference = mean_diff,
      Test = test_used,
      P_brute = p_values,
      P_ajustee = adj_p_values,
      Significatif_FDR_5pct = significance,
      row.names = NULL
    )
    
    results_file <- file.path(output_dir, "Tests_Statistiques_2_Conditions.csv")
    header_fr <- c(
      "# =====================================================================",
      "# ANALYSE STATISTIQUE DES PROPORTIONS DE CLUSTERS ENTRE DEUX CONDITIONS",
      "# =====================================================================",
      paste0("# Condition 1 : ", cond_levels[1]),
      paste0("# Condition 2 : ", cond_levels[2]),
      "# Test utilise : Welch t-test non apparie"
    )
    writeLines(header_fr, results_file)
    suppressWarnings(
      write.table(results_df, file = results_file, sep = ",",
                  row.names = FALSE, col.names = TRUE, append = TRUE)
    )
    cat("[OK] Tests statistiques (2 conditions) sauvegardes :", results_file, "\n")
    
    sig <- which(adj_p_values < 0.05 & !is.na(adj_p_values))
    if (length(sig) > 0) {
      cat("   Clusters significativement differents (FDR < 0.05) :\n")
      for (i in sig) cat(sprintf("     %s (p-ajustee = %.4f)\n", rownames(props_transformed)[i], adj_p_values[i]))
    } else {
      cat("   Aucun cluster significativement different apres correction.\n")
    }
    
  } else if (length(unique_conditions) > 2) {
    cat("[INFO] Detection automatique du design experimental (plusieurs conditions detectees)...\n")
    
    parse_condition_design <- function(condition_names) {
      temps <- sapply(condition_names, function(x) {
        if (grepl("T1|t1|V1|v1|Visit1|visit1|Baseline|baseline|BL|bl", x)) return("T1")
        if (grepl("T2|t2|V2|v2|Visit2|visit2|FollowUp|followup|FU|fu", x)) return("T2")
        return(NA_character_)
      })
      groupes <- sapply(condition_names, function(x) {
        if (grepl("Non[_\\s-]?Repondeur|Non[_\\s-]?Rep|NR|Non[_\\s-]?Responder|Non[_\\s-]?Response", x, ignore.case = TRUE)) {
          return("NonRepondeur")
        }
        if (grepl("Repondeur|Rep|Responder|Response|R$|R[_\\s-]", x, ignore.case = TRUE)) {
          return("Repondeur")
        }
        return(NA_character_)
      })
      list(temps = temps, groupes = groupes)
    }
    
    parsed <- parse_condition_design(unique_conditions)
    
    if (any(is.na(parsed$temps)) || any(is.na(parsed$groupes))) {
      cat("[WARNING] Impossible de parser automatiquement les noms de conditions.\n")
      cat("    Attendu : noms contenant un indicateur temporel (T1/T2/V1/V2)\n")
      cat("    et un indicateur de groupe (Repondeur/NonRepondeur/Responder/etc.).\n")
      cat("    Exemples valides : 'Repondeurs_T1', 'NonRepondeurs_T2', 'R_T1', 'NR_T2'.\n")
      cat("[WARNING] Aucun test statistique automatique effectue.\n")
    } else {
      time_points <- sort(unique(parsed$temps))
      cat(sprintf("[INFO] Points temporels detectes : %s\n", paste(time_points, collapse = ", ")))
      
      for (tp in time_points) {
        idx_tp <- which(parsed$temps == tp)
        conds_tp <- unique_conditions[idx_tp]
        groups_tp <- parsed$groupes[idx_tp]
        
        cat(sprintf("\n--- Analyse entre groupes pour %s ---\n", tp))
        if (length(conds_tp) != 2) {
          cat(sprintf("[WARNING] %s : nombre de conditions != 2 (%d trouvees). Test ignore.\n", tp, length(conds_tp)))
          next
        }
        if (length(unique(groups_tp)) != 2) {
          cat(sprintf("[WARNING] %s : les conditions ne forment pas deux groupes distincts. Test ignore.\n", tp))
          next
        }
        
        mask_files <- sample_conditions$Condition %in% conds_tp
        files_tp <- sample_conditions$File[mask_files]
        if (length(files_tp) == 0) {
          cat(sprintf("[WARNING] %s : aucun echantillon trouve. Test ignore.\n", tp))
          next
        }
        
        sub_props <- props_transformed[, colnames(props_transformed) %in% files_tp, drop = FALSE]
        sub_cond_vector <- sample_conditions$Condition[match(colnames(sub_props), sample_conditions$File)]
        idx_cond1 <- which(sub_cond_vector == conds_tp[1])
        idx_cond2 <- which(sub_cond_vector == conds_tp[2])
        
        n_clusters <- nrow(sub_props)
        p_values <- numeric(n_clusters)
        test_used <- character(n_clusters)
        mean_diff <- numeric(n_clusters)
        
        for (cl in seq_len(n_clusters)) {
          vals1 <- sub_props[cl, idx_cond1]
          vals2 <- sub_props[cl, idx_cond2]
          if (length(vals1) >= 2 && length(vals2) >= 2) {
            res <- t.test(vals2, vals1, var.equal = FALSE)
            p_values[cl] <- res$p.value
            test_used[cl] <- "Welch t-test"
            mean_diff[cl] <- mean(vals2) - mean(vals1)
          } else {
            p_values[cl] <- NA
            test_used[cl] <- "Effectif insuffisant (<2)"
            mean_diff[cl] <- NA
          }
        }
        adj_p_values <- p.adjust(p_values, method = "BH")
        significance <- ifelse(adj_p_values < 0.05 & !is.na(adj_p_values), "Oui", "Non")
        
        results_df <- data.frame(
          Cluster = rownames(sub_props),
          Annotation = rownames(sub_props),
          Moyenne_Transformee_Cond1 = rowMeans(sub_props[, idx_cond1, drop = FALSE]),
          Moyenne_Transformee_Cond2 = rowMeans(sub_props[, idx_cond2, drop = FALSE]),
          Difference = mean_diff,
          Test = test_used,
          P_brute = p_values,
          P_ajustee = adj_p_values,
          Significatif_FDR_5pct = significance,
          row.names = NULL
        )
        
        results_file <- file.path(output_dir, sprintf("Tests_Statistiques_%s_GroupComp.csv", tp))
        write.csv(results_df, results_file, row.names = FALSE)
        cat(sprintf("[OK] Tests entre groupes pour %s sauvegardes : %s\n", tp, results_file))
      }
      
      groups <- sort(unique(parsed$groupes[!is.na(parsed$groupes)]))
      for (grp in groups) {
        cat(sprintf("\n--- Analyse temporelle (T1 vs T2) au sein du groupe : %s ---\n", grp))
        grp_conditions <- unique_conditions[parsed$groupes == grp]
        if (length(grp_conditions) != 2) {
          cat(sprintf("[WARNING] Groupe %s : %d conditions trouvees (attendu 2). Test ignore.\n", grp, length(grp_conditions)))
          next
        }
        t1_cond <- grp_conditions[grepl("T1|t1|V1|v1|Visit1|visit1|Baseline|baseline|BL|bl", grp_conditions, ignore.case = TRUE)]
        t2_cond <- grp_conditions[grepl("T2|t2|V2|v2|Visit2|visit2|FollowUp|followup|FU|fu", grp_conditions, ignore.case = TRUE)]
        if (length(t1_cond) == 0 || length(t2_cond) == 0) {
          cat(sprintf("[WARNING] Groupe %s : impossible d'identifier T1/T2 dans les conditions.\n", grp))
          next
        }
        files_t1 <- sample_conditions$File[sample_conditions$Condition == t1_cond]
        files_t2 <- sample_conditions$File[sample_conditions$Condition == t2_cond]
        if (length(files_t1) == 0 || length(files_t2) == 0) {
          cat(sprintf("[WARNING] Groupe %s : aucun echantillon trouve pour T1 ou T2.\n", grp))
          next
        }
        sub_props <- props_transformed[, c(files_t1, files_t2), drop = FALSE]
        idx_t1 <- which(colnames(sub_props) %in% files_t1)
        idx_t2 <- which(colnames(sub_props) %in% files_t2)
        
        n_clusters <- nrow(sub_props)
        p_values <- numeric(n_clusters)
        test_used <- character(n_clusters)
        mean_diff <- numeric(n_clusters)
        
        for (cl in seq_len(n_clusters)) {
          vals_t1 <- sub_props[cl, idx_t1]
          vals_t2 <- sub_props[cl, idx_t2]
          if (length(vals_t1) >= 2 && length(vals_t2) >= 2) {
            res <- t.test(vals_t2, vals_t1, var.equal = FALSE)
            p_values[cl] <- res$p.value
            test_used[cl] <- "Welch t-test (non apparie)"
            mean_diff[cl] <- mean(vals_t2) - mean(vals_t1)
          } else {
            p_values[cl] <- NA
            test_used[cl] <- "Effectif insuffisant (<2)"
            mean_diff[cl] <- NA
          }
        }
        adj_p_values <- p.adjust(p_values, method = "BH")
        significance <- ifelse(adj_p_values < 0.05 & !is.na(adj_p_values), "Oui", "Non")
        
        results_df <- data.frame(
          Cluster = rownames(sub_props),
          Annotation = rownames(sub_props),
          Moyenne_Transformee_T1 = rowMeans(sub_props[, idx_t1, drop = FALSE]),
          Moyenne_Transformee_T2 = rowMeans(sub_props[, idx_t2, drop = FALSE]),
          Difference_T2vsT1 = mean_diff,
          Test = test_used,
          P_brute = p_values,
          P_ajustee = adj_p_values,
          Significatif_FDR_5pct = significance,
          row.names = NULL
        )
        
        results_file <- file.path(output_dir, sprintf("Tests_Statistiques_%s_T1vsT2.csv", grp))
        write.csv(results_df, results_file, row.names = FALSE)
        cat(sprintf("[OK] Tests T1 vs T2 pour le groupe %s sauvegardes : %s\n", grp, results_file))
        
        sig <- which(adj_p_values < 0.05 & !is.na(adj_p_values))
        if (length(sig) > 0) {
          cat(sprintf("   Clusters significativement differents (T2 vs T1) dans %s :\n", grp))
          for (i in sig) cat(sprintf("     %s (p-ajustee = %.4f)\n", rownames(sub_props)[i], adj_p_values[i]))
        } else {
          cat(sprintf("   Aucun cluster significativement different dans %s apres correction.\n", grp))
        }
      }
    }
  } else {
    cat("[WARNING] Nombre de conditions < 2. Aucun test statistique.\n")
  }
  
  create_proportion_heatmap_by_file(prop_data, sample_name, output_dir, conditions)
  create_stacked_bar_charts(prop_long_cond, output_dir, sample_name)
  
  if (!is.null(umap_df)) {
    umap_cond <- umap_df
    umap_cond$Condition <- condition_vec
    create_condition_dimred_plot(umap_cond, "UMAP1", "UMAP2", sample_name, output_dir, "UMAP", conditions)
  }
  if (!is.null(umap_df)) {
    create_facet_dimred_by_file(umap_df, prop_data, sample_name, output_dir, "UMAP", conditions)
  }
  
  cat("✅ Analyse par condition terminée.\n")
  return(list(prop_table = prop_table_cond, prop_long = prop_long_cond))
}
# ============================================================================
# SECTION 10: FONCTIONS DE DÉTECTION ET NORMALISATION DES BATCH EFFECTS (CytoNorm)
# ============================================================================

check_batch_effects_and_normalize <- function(expr_data, file_ids, folder_ids,
                                              channels, output_dir, sample_name,
                                              transformation_params, desc = NULL,
                                              marker_descriptions) {
  
  if (!requireNamespace("CytoNorm", quietly = TRUE)) {
    cat("Le package CytoNorm est nécessaire. Installation...\n")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("CytoNorm")
    library(CytoNorm)
  } else {
    library(CytoNorm)
  }
  
  batch_dir <- file.path(output_dir, "Batch_Effect_Analysis")
  if (!dir.exists(batch_dir)) dir.create(batch_dir)
  
  unique_files <- sort(unique(file_ids))
  n_files <- length(unique_files)
  cat("Nombre de fichiers (échantillons) détectés :", n_files, "\n")
  
  file_counts <- table(file_ids)
  cat("Effectifs par fichier :\n")
  print(file_counts)
  
  n_cells_per_file <- 5000
  
  data_list <- list()
  sample_names <- character()
  for (fid in unique_files) {
    idx <- which(file_ids == fid)
    cells <- expr_data[idx, channels, drop = FALSE]
    if (nrow(cells) > n_cells_per_file) {
      set.seed(123)
      idx_sample <- sample(1:nrow(cells), n_cells_per_file)
      cells <- cells[idx_sample, ]
    }
    data_list[[as.character(fid)]] <- cells
    sample_names <- c(sample_names, paste0("F", fid))
  }
  names(data_list) <- sample_names
  
  cat("\nDéfinition des lots (batch) :\n")
  cat("1. Utiliser les dossiers (FolderID) comme lots\n")
  cat("2. Considérer chaque fichier comme un lot séparé\n")
  cat("3. Définir manuellement les lots\n")
  repeat {
    batch_choice <- readline("Votre choix (1/2/3) : ")
    if (batch_choice %in% c("1", "2", "3")) break
  }
  
  batch_vector <- rep(NA, length(sample_names))
  if (batch_choice == "1") {
    folder_for_file <- sapply(unique_files, function(fid) {
      folder_ids[which(file_ids == fid)[1]]
    })
    batch_vector <- folder_for_file
    cat("Lots basés sur FolderID :\n")
    print(data.frame(File = unique_files, Folder = batch_vector))
  } else if (batch_choice == "2") {
    batch_vector <- seq_along(sample_names)
    cat("Chaque fichier est son propre lot.\n")
  } else {
    cat("Entrez pour chaque fichier le numéro de lot (séparé par des virgules ou espaces) :\n")
    cat("Fichiers :", paste(unique_files, collapse = ", "), "\n")
    batch_input <- readline("Lots : ")
    batch_vector <- as.numeric(unlist(strsplit(batch_input, "[, ]+")))
    if (length(batch_vector) != length(sample_names)) {
      stop("Nombre de lots incorrect.")
    }
  }
  
  nClusters <- 15
  repeat {
    cat("\nSouhaitez-vous déterminer le nombre optimal de clusters pour CytoNorm via un waterfall plot? (y/n) : ")
    test_clusters <- tolower(readline())
    if (test_clusters %in% c("y", "yes", "n", "no")) break
  }
  
  if (test_clusters %in% c("y", "yes")) {
    cat("Entrez la plage de nombres de clusters à tester (ex: 3:30) [défaut: 3:30] : ")
    range_input <- readline()
    if (range_input == "") {
      k_range <- 3:30
    } else {
      k_range <- eval(parse(text = range_input))
      k_range <- sort(unique(as.integer(k_range)))
      k_range <- k_range[k_range >= 3]
      if (length(k_range) == 0) k_range <- 3:30
    }
    
    cat("Test des clusters :", paste(k_range, collapse = ", "), "\n")
    
    unique_batches <- sort(unique(batch_vector))
    batch_agg_files <- character(length(unique_batches))
    temp_agg_dir <- file.path(batch_dir, "agg_for_cv")
    if (!dir.exists(temp_agg_dir)) dir.create(temp_agg_dir)
    
    for (i in seq_along(unique_batches)) {
      b <- unique_batches[i]
      batch_sample_indices <- which(batch_vector == b)
      batch_cells_list <- data_list[batch_sample_indices]
      batch_cells <- do.call(rbind, batch_cells_list)
      if (nrow(batch_cells) > 10000) {
        set.seed(123 + i)
        batch_cells <- batch_cells[sample(1:nrow(batch_cells), 10000), ]
      }
      ff <- flowFrame(batch_cells)
      agg_file <- file.path(temp_agg_dir, paste0("batch_", b, ".fcs"))
      write.FCS(ff, agg_file)
      batch_agg_files[i] <- agg_file
    }
    
    cat("Préparation du modèle FlowSOM pour le calcul des CV...\n")
    fsom <- CytoNorm::prepareFlowSOM(
      files = batch_agg_files,
      colsToUse = channels,
      nCells = 20000,
      FlowSOM.params = list(xdim = 10, ydim = 10, nClus = max(k_range), scale = FALSE),
      transformList = NULL,
      seed = 123,
      verbose = TRUE
    )
    
    cat("Calcul du waterfall plot...\n")
    cv_results <- CytoNorm::testCV(
      fsom = fsom,
      cluster_values = k_range,
      plot = FALSE,
      verbose = TRUE,
      seed = 123
    )
    
    waterfall_file <- file.path(batch_dir, "Waterfall_Plot.pdf")
    pdf(waterfall_file, width = 10, height = 8)
    CytoNorm::PlotOverviewCV(
      fsom = fsom,
      cv_res = cv_results,
      max_cv = 1,
      show_cv = 0.5
    )
    dev.off()
    
    cat("Waterfall plot généré :", waterfall_file, "\n")
    if (interactive()) open_pdf(waterfall_file)
    
    repeat {
      cat("Entrez le nombre de clusters choisi (parmi ceux testés) : ")
      chosen_k <- as.numeric(readline())
      if (!is.na(chosen_k) && chosen_k %in% k_range) break
      cat("Valeur invalide. Veuillez choisir parmi :", paste(k_range, collapse = ", "), "\n")
    }
    nClusters <- chosen_k
    cat("Nombre de clusters sélectionné :", nClusters, "\n")
    
    unlink(temp_agg_dir, recursive = TRUE)
  } else {
    cat("Utilisation de la valeur par défaut : nClusters = 15\n")
  }
  
  cat("\nEntraînement du modèle CytoNorm...\n")
  
  temp_dir_train <- tempfile(pattern = "cytonorm_train_")
  dir.create(temp_dir_train)
  fcs_files_train <- character(length(data_list))
  
  for (i in seq_along(data_list)) {
    mat <- data_list[[i]]
    ff <- flowFrame(mat)
    file_path <- file.path(temp_dir_train, paste0("sample_", i, ".fcs"))
    write.FCS(ff, file_path)
    fcs_files_train[i] <- file_path
  }
  
  model <- CytoNorm::CytoNorm.train(
    files = fcs_files_train,
    labels = sample_names,
    channels = channels,
    transformList = NULL,
    seed = 123,
    quantileValues = c(0.01, 0.5, 0.99),
    nClusters = nClusters,
    goal = "mean",
    verbose = TRUE
  )
  
  unlink(temp_dir_train, recursive = TRUE)
  
  cat("\nPréparation des données complètes pour la normalisation...\n")
  full_data_list <- list()
  for (fid in unique_files) {
    idx <- which(file_ids == fid)
    full_data_list[[as.character(fid)]] <- expr_data[idx, channels, drop = FALSE]
  }
  names(full_data_list) <- sample_names
  
  cat("\nApplication de la normalisation à toutes les cellules...\n")
  
  temp_dir_norm <- tempfile(pattern = "cytonorm_norm_")
  dir.create(temp_dir_norm)
  fcs_files_norm <- character(length(full_data_list))
  
  for (i in seq_along(full_data_list)) {
    mat <- full_data_list[[i]]
    ff <- flowFrame(mat)
    file_path <- file.path(temp_dir_norm, paste0("sample_", i, ".fcs"))
    write.FCS(ff, file_path)
    fcs_files_norm[i] <- file_path
  }
  
  CytoNorm::CytoNorm.normalize(
    model = model,
    files = fcs_files_norm,
    labels = sample_names,
    transformList = NULL,
    transformList.reverse = NULL,
    verbose = TRUE,
    outputDir = temp_dir_norm,
    prefix = "norm_"
  )
  
  norm_data_list <- list()
  for (i in seq_along(fcs_files_norm)) {
    norm_file <- file.path(temp_dir_norm, paste0("norm_sample_", i, ".fcs"))
    if (!file.exists(norm_file)) {
      stop("Fichier normalisé introuvable : ", norm_file)
    }
    ff_norm <- read.FCS(norm_file)
    norm_data_list[[i]] <- exprs(ff_norm)
  }
  names(norm_data_list) <- sample_names
  
  unlink(temp_dir_norm, recursive = TRUE)
  
  expr_data_norm <- expr_data
  for (i in seq_along(sample_names)) {
    fid <- unique_files[i]
    idx <- which(file_ids == fid)
    expr_data_norm[idx, channels] <- norm_data_list[[i]]
  }
  
  cat("\n📊 Génération des graphiques de densité avant/après normalisation...\n")
  
  temp_density_dir <- file.path(batch_dir, "temp_density_plots")
  if (!dir.exists(temp_density_dir)) dir.create(temp_density_dir)
  
  original_data_for_density <- list()
  normalized_data_for_density <- list()
  
  for (i in seq_along(data_list)) {
    sample_id <- names(data_list)[i]
    mat_orig <- data_list[[i]]
    mat_norm <- norm_data_list[[i]]
    n_orig <- nrow(mat_orig)
    n_norm <- nrow(mat_norm)
    n_min <- min(n_orig, n_norm)
    
    if (n_min < 100) {
      warning("Échantillon ", sample_id, " a moins de 100 cellules après appariement.")
    }
    
    set.seed(123 + i)
    idx_orig <- sample(1:n_orig, n_min)
    idx_norm <- sample(1:n_norm, n_min)
    
    original_data_for_density[[sample_id]] <- mat_orig[idx_orig, , drop = FALSE]
    normalized_data_for_density[[sample_id]] <- mat_norm[idx_norm, , drop = FALSE]
  }
  
  original_files_for_density <- c()
  normalized_files_for_density <- c()
  
  for (sample_id in names(original_data_for_density)) {
    ff_orig <- flowFrame(original_data_for_density[[sample_id]])
    temp_file_orig <- file.path(temp_density_dir, paste0(sample_id, ".fcs"))
    write.FCS(ff_orig, temp_file_orig)
    original_files_for_density <- c(original_files_for_density, temp_file_orig)
    
    ff_norm <- flowFrame(normalized_data_for_density[[sample_id]])
    temp_file_norm <- file.path(temp_density_dir, paste0(sample_id, "_norm.fcs"))
    write.FCS(ff_norm, temp_file_norm)
    normalized_files_for_density <- c(normalized_files_for_density, temp_file_norm)
  }
  
  input_list <- list()
  for (batch_name in unique(batch_vector)) {
    batch_indices <- which(batch_vector == batch_name)
    input_list[[paste0("Batch_", batch_name, "_original")]] <- original_files_for_density[batch_indices]
    input_list[[paste0("Batch_", batch_name, "_norm")]] <- normalized_files_for_density[batch_indices]
  }
  
  if (!requireNamespace("PeacoQC", quietly = TRUE)) {
    cat("Le package PeacoQC est nécessaire pour les graphiques de densité. Installation...\n")
    BiocManager::install("PeacoQC")
    library(PeacoQC)
  } else {
    library(PeacoQC)
  }
  
  density_plots <- tryCatch({
    plotDensities(
      input = input_list,
      channels = channels,
      transformList = NULL,
      show_goal = FALSE,
      suffix = c(original = "_original", normalized = "_norm")
    )
  }, error = function(e) {
    cat("⚠️  Erreur lors de la génération des graphiques de densité :", e$message, "\n")
    NULL
  })
  
  if (!is.null(density_plots) && length(density_plots) > 0) {
    density_pdf <- file.path(batch_dir, "Density_Plots_Before_After.pdf")
    n_plots <- length(density_plots)
    plots_per_page <- 2
    n_pages <- ceiling(n_plots / plots_per_page)
    
    pdf(density_pdf, width = 8.27, height = 11.69)
    for (page in seq_len(n_pages)) {
      idx <- ((page - 1) * plots_per_page + 1):min(page * plots_per_page, n_plots)
      page_plots <- density_plots[idx]
      
      if (length(page_plots) == 1) {
        print(page_plots[[1]])
      } else {
        combined <- ggpubr::ggarrange(plotlist = page_plots, ncol = 2, nrow = 1)
        print(combined)
      }
    }
    dev.off()
    cat("✅ Graphiques de densité sauvegardés :", density_pdf, "\n")
    if (interactive()) open_pdf(density_pdf)
  }
  
  unlink(temp_density_dir, recursive = TRUE)
  cat("✅ Nettoyage des fichiers temporaires terminé.\n")
  
  if (exists("emdEvaluation") && exists("madEvaluation")) {
    cat("\n📊 Calcul des métriques EMD et MAD...\n")
    
    temp_eval_dir <- file.path(batch_dir, "temp_eval")
    if (!dir.exists(temp_eval_dir)) dir.create(temp_eval_dir)
    
    original_files_full <- c()
    for (i in seq_along(full_data_list)) {
      sample_id <- names(full_data_list)[i]
      ff <- flowFrame(full_data_list[[i]])
      temp_file <- file.path(temp_eval_dir, paste0("original_", sample_id, ".fcs"))
      write.FCS(ff, temp_file)
      original_files_full <- c(original_files_full, temp_file)
    }
    
    normalized_files_full <- c()
    for (i in seq_along(norm_data_list)) {
      sample_id <- names(norm_data_list)[i]
      ff <- flowFrame(norm_data_list[[i]])
      temp_file <- file.path(temp_eval_dir, paste0("norm_", sample_id, ".fcs"))
      write.FCS(ff, temp_file)
      normalized_files_full <- c(normalized_files_full, temp_file)
    }
    
    safe_emd_mad_evaluation <- function(files, channels, prefix, type = "emd") {
      all_data <- list()
      for (f in files) {
        ff <- read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
        expr <- exprs(ff)
        for (ch in channels) {
          if (ch %in% colnames(expr)) {
            q_low <- quantile(expr[, ch], probs = 0.005, na.rm = TRUE)
            q_high <- quantile(expr[, ch], probs = 0.995, na.rm = TRUE)
            expr[, ch] <- pmax(pmin(expr[, ch], q_high), q_low)
          }
        }
        all_data[[f]] <- expr
      }
      
      temp_dir <- tempfile(pattern = "adjusted_")
      dir.create(temp_dir)
      adjusted_files <- character(length(files))
      for (i in seq_along(files)) {
        ff <- flowFrame(all_data[[files[i]]])
        adj_file <- file.path(temp_dir, basename(files[i]))
        write.FCS(ff, adj_file)
        adjusted_files[i] <- adj_file
      }
      
      result <- NULL
      tryCatch({
        if (type == "emd") {
          result <- emdEvaluation(
            files = adjusted_files,
            channels = channels,
            transformList = NULL,
            prefix = prefix,
            return_all = FALSE
          )
        } else {
          result <- madEvaluation(
            files = adjusted_files,
            channels = channels,
            transformList = NULL,
            prefix = prefix,
            return_all = FALSE
          )
        }
      }, error = function(e) {
        cat("⚠️  Erreur lors du calcul (après ajustement) :", e$message, "\n")
      })
      
      unlink(temp_dir, recursive = TRUE)
      return(result)
    }
    
    emd_before <- NULL
    emd_after <- NULL
    tryCatch({
      emd_before <- safe_emd_mad_evaluation(
        files = original_files_full,
        channels = channels,
        prefix = "^original_",
        type = "emd"
      )
      emd_after <- safe_emd_mad_evaluation(
        files = normalized_files_full,
        channels = channels,
        prefix = "^norm_",
        type = "emd"
      )
    }, error = function(e) {
      cat("⚠️  Erreur lors du calcul de l'EMD :", e$message, "\n")
    })
    
    mad_before <- NULL
    mad_after <- NULL
    tryCatch({
      mad_before <- safe_emd_mad_evaluation(
        files = original_files_full,
        channels = channels,
        prefix = "^original_",
        type = "mad"
      )
      mad_after <- safe_emd_mad_evaluation(
        files = normalized_files_full,
        channels = channels,
        prefix = "^norm_",
        type = "mad"
      )
    }, error = function(e) {
      cat("⚠️  Erreur lors du calcul du MAD :", e$message, "\n")
    })
    
    unlink(temp_eval_dir, recursive = TRUE)
    
    if (!is.null(emd_before) && !is.null(emd_after) && !is.null(mad_before) && !is.null(mad_after)) {
      eval_pdf <- file.path(batch_dir, "EMD_MAD_Evaluation.pdf")
      pdf(eval_pdf, width = 10, height = 5)
      par(mfrow = c(1, 2))
      
      plot(emd_before, emd_after,
           xlab = "EMD avant normalisation",
           ylab = "EMD après normalisation",
           main = "Évaluation EMD",
           pch = 19, col = rgb(0,0,1,0.5))
      abline(0, 1, col = "red", lty = 2)
      text(emd_before, emd_after, labels = rownames(emd_before), cex=0.7, pos=4)
      
      plot(mad_before, mad_after,
           xlab = "MAD avant normalisation",
           ylab = "MAD après normalisation",
           main = "Évaluation MAD",
           pch = 19, col = rgb(0,0,1,0.5))
      abline(0, 1, col = "red", lty = 2)
      
      dev.off()
      cat("✅ Graphiques EMD/MAD sauvegardés :", eval_pdf, "\n")
      if (interactive()) open_pdf(eval_pdf)
    }
  }
  
  cat("✅ Normalisation des batch effects terminée.\n")
  return(expr_data_norm)
}

# ============================================================================
# SECTION 11: FONCTION PRINCIPALE D'ANALYSE
# ============================================================================
analyze_fcs_file <- function(fcs_file, output_dir_base,
                             doublet_gate_coords = NULL,
                             fsc_ssa_gate_coords = NULL,
                             transformation_params = NULL,
                             script_name = NULL,
                             n_max_cells = NULL,
                             log_temp_file = NULL) {
  
  analysis_metrics <- list()
  analysis_metrics$start_time <- Sys.time()
  analysis_metrics$seed <- 123
  
  analysis_metrics$flowsom_xdim <- 10
  analysis_metrics$flowsom_ydim <- 10
  analysis_metrics$flowsom_rlen <- 100
  
  analysis_metrics$umap_n_neighbors <- 40
  analysis_metrics$umap_min_dist <- 0.15
  analysis_metrics$umap_metric <- "euclidean"
  analysis_metrics$umap_n_epochs <- 300
  analysis_metrics$umap_n_threads <- 9
  
  sample_name <- tools::file_path_sans_ext(basename(fcs_file))
  output_dir <- file.path(output_dir_base, sample_name)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  if (!is.null(log_temp_file)) {
    log_file <- file.path(output_dir, paste0(sample_name, "_console.log"))
    while (sink.number() > 0) sink()
    file.copy(log_temp_file, log_file, overwrite = TRUE)
    unlink(log_temp_file)
    sink(log_file, split = TRUE, append = TRUE)
    cat("📄 Log console déplacé vers :", log_file, "\n")
  }
  
  cat("\n", paste(rep("-", 60), collapse = ""), "\n", sep = "")
  cat("ANALYSE DE L'ÉCHANTILLON:", sample_name, "\n")
  cat(paste(rep("-", 60), collapse = ""), "\n", sep = "")
  
  cat("📁 Étape 1/12: Chargement du fichier FCS...\n")
  fs <- read.FCS(fcs_file, transformation = FALSE, truncate_max_range = FALSE)
  analysis_metrics$initial_cell_count <- nrow(fs)
  cat("✅ Fichier FCS chargé (", analysis_metrics$initial_cell_count, "cellules)\n\n")
  
  desc <- description(fs)
  
  current_names <- colnames(exprs(fs))
  new_names <- gsub("-", ".", current_names)
  colnames(fs) <- new_names
  
  params <- parameters(fs)
  marker_descriptions <- pData(params)$desc
  names(marker_descriptions) <- new_names
  
  cat("🧠 Étape 2/12: Définition des paramètres de transformation...\n")
  non_fluorescent_channels <- c("FSC.A", "FSC.H", "FSC.W", "SSC.H", "SSC.W",
                                "SSC.B.A", "SSC.B.H", "Time", "AF.A", "FolderID", "FileID")
  all_channels <- colnames(fs)
  channels_to_transform <- setdiff(all_channels, non_fluorescent_channels)
  
  if (is.null(transformation_params)) {
    transformation_params <- list()
    for (channel in channels_to_transform) {
      transformation_params[[channel]] <- list(a = 0, b = 1/1000)
    }
    
    default_hist_file <- file.path(output_dir, "transformation_histograms_initial.pdf")
    show_fluorescence_histograms(exprs(fs), channels_to_transform,
                                 marker_descriptions, transformation_params,
                                 filename = default_hist_file,
                                 title_suffix = "(paramètres par défaut)", sample_name = sample_name)
    
    if (file.exists(default_hist_file)) open_pdf(default_hist_file)
    
    cat("\nLe PDF des histogrammes a été généré. Veuillez le visualiser.\n")
    cat("Voulez-vous modifier les paramètres de transformation? (y/n): ")
    modify_params <- readline()
    
    if (tolower(modify_params) %in% c("y", "yes")) {
      for (channel in channels_to_transform) {
        repeat {
          cat(paste("\nParamètres actuels pour", channel, ":\n"))
          cat(paste("a =", transformation_params[[channel]]$a, ", b =", transformation_params[[channel]]$b, "\n"))
          
          create_plot_window(8, 6)
          show_single_channel_histogram(exprs(fs), channel,
                                        marker_descriptions[channel],
                                        transformation_params[[channel]]$a,
                                        transformation_params[[channel]]$b,
                                        sample_name, "(paramètres actuels)")
          
          cat(paste("Voulez-vous modifier les paramètres pour", channel, "? (y/n): "))
          modify_channel <- readline()
          
          if (tolower(modify_channel) %in% c("y", "yes")) {
            cat(paste("Entrez la valeur de a pour", channel, "(décalage, par défaut 0): "))
            a_input <- readline()
            if (a_input == "") a_input <- "0"
            a <- as.numeric(a_input)
            
            cat(paste("Entrez la valeur de b pour", channel, "(facteur d'échelle, par défaut 1/1000): "))
            b_input <- readline()
            if (b_input == "") b_input <- "1/1000"
            b <- tryCatch(eval(parse(text = b_input)), error = function(e) as.numeric(b_input))
            if (is.na(b)) b <- 1/150
            
            show_single_channel_histogram(exprs(fs), channel,
                                          marker_descriptions[channel], a, b,
                                          sample_name, paste("(a =", a, ", b =", b, ")"))
            
            cat("Êtes-vous satisfait de ces paramètres? (y/n): ")
            satisfied <- readline()
            
            if (tolower(satisfied) %in% c("y", "yes")) {
              transformation_params[[channel]] <- list(a = a, b = b)
              dev.off()
              break
            }
          } else {
            dev.off()
            break
          }
          dev.off()
        }
      }
      
      final_hist_file <- file.path(output_dir, "transformation_histograms_final.pdf")
      show_fluorescence_histograms(exprs(fs), channels_to_transform,
                                   marker_descriptions, transformation_params,
                                   filename = final_hist_file,
                                   title_suffix = "(paramètres finaux)", sample_name = sample_name)
      if (file.exists(final_hist_file)) open_pdf(final_hist_file)
    }
  }
  
  analysis_metrics$transformation_params <- transformation_params
  cat("✅ Paramètres de transformation définis\n")
  
  cat("🔧 Étape 3/12: Transformation des données...\n")
  fs_trans <- fs
  if (length(channels_to_transform) > 0) {
    trans_list <- list()
    for (channel in channels_to_transform) {
      a <- transformation_params[[channel]]$a
      b <- transformation_params[[channel]]$b
      asinh_trans <- arcsinhTransform(transformationId = paste0("arcsinh_", channel), a = a, b = b, c = 0)
      trans_list[[channel]] <- asinh_trans
    }
    fs_trans <- transform(fs_trans, transformList(channels_to_transform, trans_list))
  }
  cat("✅ Transformation appliquée\n")
  
  cat("🔬 Étape 4/12: Sélection interactive de populations avant clustering...\n")
  
  fluoro_channels_available <- setdiff(channels_to_transform,
                                       c(non_fluorescent_channels))
  
  fluoro_table <- data.frame(
    num     = seq_along(fluoro_channels_available),
    channel = fluoro_channels_available,
    marker  = sapply(fluoro_channels_available,
                     function(ch) get_marker_name(ch, marker_descriptions)),
    stringsAsFactors = FALSE
  )
  
  current_mask <- rep(TRUE, nrow(fs_trans))
  gate_fluoro_channels_used <- character(0)
  population_description <- NULL
  gate_counter <- 0
  
  cat("\nVoulez-vous sélectionner une population précise avant le clustering? (y/n): ")
  do_gate <- readline()
  
  while (tolower(do_gate) %in% c("y", "yes", "o", "oui")) {
    gate_counter <- gate_counter + 1
    
    cat("\n--- Canaux fluorescents disponibles ---\n")
    cat(paste(rep("-", 50), collapse = ""), "\n")
    for (i in seq_len(nrow(fluoro_table))) {
      cat(sprintf("%2d. %-25s (Marqueur: %s)\n",
                  fluoro_table$num[i],
                  fluoro_table$channel[i],
                  fluoro_table$marker[i]))
    }
    cat(paste(rep("-", 50), collapse = ""), "\n")
    
    x_choice <- NA_integer_
    while (is.na(x_choice) || x_choice < 1 || x_choice > nrow(fluoro_table)) {
      cat("Entrez le numéro du canal pour l'axe X (dot plot): ")
      x_input <- readline()
      x_choice <- suppressWarnings(as.integer(x_input))
    }
    x_channel <- fluoro_table$channel[x_choice]
    x_label   <- fluoro_table$marker[x_choice]
    
    y_choice <- NA_integer_
    while (is.na(y_choice) || y_choice < 1 || y_choice > nrow(fluoro_table)) {
      cat("Entrez le numéro du canal pour l'axe Y (dot plot): ")
      y_input <- readline()
      y_choice <- suppressWarnings(as.integer(y_input))
    }
    y_channel <- fluoro_table$channel[y_choice]
    y_label   <- fluoro_table$marker[y_choice]
    
    n_before_gate <- sum(current_mask)
    cat(sprintf("\n📊 Cellules avant le gate %d : %d\n", gate_counter, n_before_gate))
    
    cells_in_current <- which(current_mask)
    num_points_gate <- min(20000, length(cells_in_current))
    plot_idx <- sample(cells_in_current, num_points_gate)
    plot_data_gate <- exprs(fs_trans)[plot_idx, c(x_channel, y_channel), drop = FALSE]
    
    create_plot_window(8, 6)
    plot(plot_data_gate[, x_channel], plot_data_gate[, y_channel],
         pch = 16, cex = 0.6,
         main = paste0(sample_name, " - Gate ", gate_counter,
                       ": ", x_label, " vs ", y_label,
                       "\nDessinez le gate sur le dot plot"),
         xlab = paste0(x_channel, " (", x_label, ")"),
         ylab = paste0(y_channel, " (", y_label, ")"),
         col  = rgb(0, 0, 1, 0.3),
         xlim = range(plot_data_gate[, x_channel]),
         ylim = range(plot_data_gate[, y_channel]))
    
    cat(sprintf("Veuillez tracer le gate %d sur le dot plot (%s vs %s).\n",
                gate_counter, x_label, y_label))
    gate_coords_g <- locator(type = "l", col = "red", lwd = 2)
    
    if (!is.null(gate_coords_g) && length(gate_coords_g$x) >= 3) {
      gate_mat <- cbind(gate_coords_g$x, gate_coords_g$y)
      colnames(gate_mat) <- c(x_channel, y_channel)
      gate_poly <- polygonGate(.gate = gate_mat,
                               filterId = paste0("gate_", gate_counter, "_", x_channel, "_", y_channel))
      
      mask_temp <- suppressWarnings(flowCore::filter(fs_trans[cells_in_current, ], gate_poly)@subSet)
      mask_temp[is.na(mask_temp)] <- FALSE
      
      new_mask <- rep(FALSE, nrow(fs_trans))
      new_mask[cells_in_current[mask_temp]] <- TRUE
      current_mask <- new_mask
      
      n_after_gate <- sum(current_mask)
      cat(sprintf("✅ Cellules après le gate %d : %d (retenues: %.1f%%)\n",
                  gate_counter, n_after_gate,
                  100 * n_after_gate / n_before_gate))
      
      gate_pdf <- file.path(output_dir,
                            paste0("Gate_", gate_counter, "_",
                                   x_channel, "_vs_", y_channel, ".pdf"))
      pdf(gate_pdf, width = 8, height = 6)
      plot(plot_data_gate[, x_channel], plot_data_gate[, y_channel],
           pch = 16, cex = 0.6,
           main = paste0(sample_name, " - Gate ", gate_counter,
                         ": ", x_label, " vs ", y_label),
           xlab = paste0(x_channel, " (", x_label, ")"),
           ylab = paste0(y_channel, " (", y_label, ")"),
           col  = rgb(0, 0, 1, 0.3))
      lines(c(gate_coords_g$x, gate_coords_g$x[1]),
            c(gate_coords_g$y, gate_coords_g$y[1]),
            col = "red", lwd = 2)
      dev.off()
      cat("✅ Dot plot du gate sauvegardé :", gate_pdf, "\n")
      
      gate_fluoro_channels_used <- unique(c(gate_fluoro_channels_used,
                                            x_channel, y_channel))
    } else {
      cat("⚠️ Aucun gate valide tracé. Utilisation de toutes les cellules du masque courant.\n")
    }
    dev.off()
    
    cat("\nVoulez-vous effectuer une autre sélection? (y/n): ")
    do_gate <- readline()
  }
  
  if (gate_counter > 0) {
    cat("\nVeuillez décrire la/les population(s) sélectionnée(s) par vos gates\n")
    cat("(ex: CD3+ CD8+ CD4-, CD56+ CD3-, etc.) : ")
    population_description <- readline()
    if (nchar(trimws(population_description)) == 0) population_description <- NULL
  }
  
  gated_cells_mask <- current_mask
  n_after_all_fluoro_gates <- sum(gated_cells_mask)
  cat(sprintf("\n✅ Sélection des populations terminée (%d gates réalisés, %d cellules retenues)\n",
              gate_counter, n_after_all_fluoro_gates))
  
  if (n_after_all_fluoro_gates == 0) {
    cat("❌ Aucune cellule après les gates fluorescents. Arrêt de l'analyse.\n")
    return(NULL)
  }
  
  fs_gated <- fs[gated_cells_mask, ]
  output_fcs_file <- file.path(output_dir, paste0(sample_name, "_target_population.fcs"))
  write.FCS(fs_gated, filename = output_fcs_file)
  cat("✅ Population cible sauvegardée (données brutes non transformées) :", output_fcs_file, "\n")
  
  analysis_metrics$after_gate_CD8_CD56  <- n_after_all_fluoro_gates
  analysis_metrics$after_gate_CD3_UnconvT <- n_after_all_fluoro_gates
  
  cat("🔬 Étape 5/12: Gating des doublets...\n")
  final_cells_mask <- gated_cells_mask
  
  if ("FSC.A" %in% colnames(fs_trans) && "FSC.H" %in% colnames(fs_trans)) {
    if (is.null(doublet_gate_coords)) {
      create_plot_window(8, 6)
      expr_subset_for_doublet_plot <- exprs(fs_trans)[final_cells_mask, ]
      num_points_doublet_plot <- min(20000, nrow(expr_subset_for_doublet_plot))
      doublet_plot_indices <- sample(1:nrow(expr_subset_for_doublet_plot), num_points_doublet_plot)
      doublet_plot_data <- expr_subset_for_doublet_plot[doublet_plot_indices, ]
      
      plot(doublet_plot_data[, "FSC.A"], doublet_plot_data[, "FSC.H"],
           pch = 16, cex = 0.5,
           main = paste(sample_name, "- Gating des Doublets (FSC.A vs FSC.H)"),
           xlab = "FSC.A", ylab = "FSC.H", col = rgb(0, 0, 1, 0.3))
      
      cat("Veuillez tracer le gate des singlets sur le graphique.\n")
      gate_coords <- locator(type = "l", col = "green", lwd = 2)
      
      if (!is.null(gate_coords) && length(gate_coords$x) >= 3) {
        doublet_gate_coords <- list(x = gate_coords$x, y = gate_coords$y)
        pdf(file.path(output_dir, "Doublet_Gating.pdf"), width = 8, height = 6)
        plot(doublet_plot_data[, "FSC.A"], doublet_plot_data[, "FSC.H"],
             pch = 16, cex = 0.5,
             main = paste(sample_name, "- Gating des Doublets (FSC.A vs FSC.H)"),
             xlab = "FSC.A", ylab = "FSC.H", col = rgb(0, 0, 1, 0.3))
        lines(c(gate_coords$x, gate_coords$x[1]), c(gate_coords$y, gate_coords$y[1]),
              col = "green", lwd = 2)
        dev.off()
      } else {
        doublet_gate_coords <- NULL
      }
      dev.off()
    }
    
    if (!is.null(doublet_gate_coords)) {
      doublet_gate <- polygonGate(.gate = data.frame("FSC.A" = doublet_gate_coords$x,
                                                     "FSC.H" = doublet_gate_coords$y),
                                  filterId = "doublet_gate")
      indices_gated_cells <- which(final_cells_mask)
      singlets_mask <- suppressWarnings(flowCore::filter(fs_trans[indices_gated_cells, ], doublet_gate)@subSet)
      
      updated_mask <- rep(FALSE, length(final_cells_mask))
      updated_mask[indices_gated_cells[singlets_mask]] <- TRUE
      final_cells_mask <- updated_mask
    }
  }
  
  analysis_metrics$after_doublet_removal <- sum(final_cells_mask)
  cat("✅ Gating des doublets terminé (", analysis_metrics$after_doublet_removal, "cellules retenues)\n")
  
  cat("🔬 Étape 6/12: Gating FSC-A vs SSC-A...\n")
  final_cells_mask_after_ssa <- final_cells_mask
  
  if ("FSC.A" %in% colnames(fs_trans) && "SSC.A" %in% colnames(fs_trans)) {
    create_plot_window(8, 6)
    expr_data <- exprs(fs_trans)[final_cells_mask, ]
    num_points_plot <- min(20000, nrow(expr_data))
    plot_indices <- sample(1:nrow(expr_data), num_points_plot)
    plot_data <- expr_data[plot_indices, c("FSC.A", "SSC.A")]
    
    plot(plot_data[, "FSC.A"], plot_data[, "SSC.A"],
         pch = "16", cex = 0.5,
         main = paste(sample_name, "- Gating FSC.A vs SSC.A"),
         xlab = "FSC.A", ylab = "SSC.A", col = rgb(0, 0, 1, 0.1))
    
    cat("Veuillez tracer le gate sur le graphique FSC-A vs SSC-A.\n")
    gate_coords <- locator(type = "l", col = "green", lwd = 2)
    
    if (!is.null(gate_coords) && length(gate_coords$x) >= 3) {
      fsc_ssa_gate_coords <- list(x = gate_coords$x, y = gate_coords$y)
      pdf(file.path(output_dir, "FSC_SSC_Gating.pdf"), width = 8, height = 6)
      plot(plot_data[, "FSC.A"], plot_data[, "SSC.A"],
           pch = ".", cex = 0.5,
           main = paste(sample_name, "- Gating FSC.A vs SSC.A"),
           xlab = "FSC.A", ylab = "SSC.A", col = rgb(0, 0, 1, 0.1))
      lines(c(gate_coords$x, gate_coords$x[1]), c(gate_coords$y, gate_coords$y[1]),
            col = "green", lwd = 2)
      dev.off()
    } else {
      fsc_ssa_gate_coords <- NULL
    }
    dev.off()
    
    if (!is.null(fsc_ssa_gate_coords)) {
      fsc_ssa_gate <- polygonGate(.gate = data.frame("FSC.A" = fsc_ssa_gate_coords$x,
                                                     "SSC.A" = fsc_ssa_gate_coords$y),
                                  filterId = "fsc_ssa_gate")
      indices_gated_cells <- which(final_cells_mask)
      ssa_mask <- suppressWarnings(flowCore::filter(fs_trans[indices_gated_cells, ], fsc_ssa_gate)@subSet)
      
      updated_mask <- rep(FALSE, length(final_cells_mask))
      updated_mask[indices_gated_cells[ssa_mask]] <- TRUE
      final_cells_mask_after_ssa <- updated_mask
    }
  }
  
  analysis_metrics$after_fsc_ssc_gating <- sum(final_cells_mask_after_ssa)
  cat("✅ Gating FSC-A vs SSC-A terminé (", analysis_metrics$after_fsc_ssc_gating, "cellules retenues)\n")
  
  channels_to_exclude <- gate_fluoro_channels_used
  
  if (length(channels_to_exclude) > 0) {
    cat("\n📋 Canaux automatiquement exclus du FlowSOM (utilisés pour les gates fluorescents) :\n")
    for (ch in channels_to_exclude) {
      cat(sprintf("  - %s (%s)\n", ch, get_marker_name(ch, marker_descriptions)))
    }
  }
  
  cat("\n🎯 EXCLUSION D'AUTRES CANAUX DE L'ANALYSE\n")
  available_channels <- setdiff(channels_to_transform, c(non_fluorescent_channels))
  available_channels <- setdiff(available_channels, channels_to_exclude)
  
  more_channels_to_exclude <- TRUE
  while (more_channels_to_exclude && length(available_channels) > 0) {
    cat("\nVoulez-vous exclure un canal de l'analyse? (y/n): ")
    exclude_another <- readline()
    if (tolower(exclude_another) %in% c("y", "yes")) {
      cat("\nCanaux disponibles pour exclusion:\n")
      cat(paste(rep("-", 50), collapse = ""), "\n")
      for (i in 1:length(available_channels)) {
        channel <- available_channels[i]
        marker_name <- get_marker_name(channel, marker_descriptions)
        cat(sprintf("%2d. %-25s (Marqueur: %s)\n", i, channel, marker_name))
      }
      cat(paste(rep("-", 50), collapse = ""), "\n")
      cat("Entrez le numéro du canal à exclure (ou 0 pour annuler): ")
      channel_choice <- as.integer(readline())
      if (!is.na(channel_choice) && channel_choice > 0 && channel_choice <= length(available_channels)) {
        channel_to_exclude <- available_channels[channel_choice]
        channels_to_exclude <- c(channels_to_exclude, channel_to_exclude)
        available_channels <- setdiff(available_channels, channel_to_exclude)
        marker_name <- get_marker_name(channel_to_exclude, marker_descriptions)
        cat(sprintf("✅ Canal '%s' (%s) ajouté à la liste des canaux exclus.\n",
                    channel_to_exclude, marker_name))
        if (length(available_channels) == 0) {
          cat("⚠️ Tous les canaux disponibles ont été exclus. Fin de la sélection.\n")
          more_channels_to_exclude <- FALSE
        }
      } else if (channel_choice == 0) {
        cat("Annulation de la sélection de canal.\n")
      } else {
        cat("❌ Sélection invalide. Veuillez choisir un numéro valide.\n")
      }
    } else {
      more_channels_to_exclude <- FALSE
      cat("Fin de la sélection des canaux à exclure.\n")
    }
  }
  
  if (length(channels_to_exclude) > 0) {
    cat("\n📋 RÉSUMÉ DES CANAUX EXCLUS DE L'ANALYSE:\n")
    for (channel in channels_to_exclude) {
      marker_name <- get_marker_name(channel, marker_descriptions)
      cat(sprintf(" - %s (%s)\n", channel, marker_name))
    }
  } else {
    cat("\n✅ Aucun canal supplémentaire exclu de l'analyse.\n")
  }
  analysis_metrics$excluded_channels <- channels_to_exclude
  
  cat("💻️ Étape 7/12: Sous-échantillonnage des données...\n")
  expr_data <- exprs(fs_trans)[final_cells_mask_after_ssa, ]
  total_cells_available <- nrow(expr_data)
  cat(sprintf("Nombre total de cellules disponibles après filtrage: %d\n", total_cells_available))
  
  has_metadata <- ("FolderID" %in% colnames(expr_data)) && ("FileID" %in% colnames(expr_data))
  
  if (has_metadata) {
    cat("\n📋 Aperçu des cellules disponibles :\n")
    if ("FolderID" %in% colnames(expr_data)) {
      cat("Par FolderID (condition) :\n")
      print(table(expr_data[, "FolderID"]))
    }
    if ("FileID" %in% colnames(expr_data)) {
      cat("\nPar FileID (fichier) :\n")
      print(table(expr_data[, "FileID"]))
    }
    
    cat("\nChoisissez le mode de sous-échantillonnage :\n")
    cat("  1. Global (échantillonnage uniforme sur toutes les cellules)\n")
    cat("  2. Équilibré par condition (même nombre de cellules par condition)\n")
    cat("  3. Équilibré par fichier (même nombre de cellules par fichier)\n")
    cat("Votre choix (1/2/3) [1] : ")
    choice <- readline()
    if (choice == "") choice <- "1"
  } else {
    choice <- "1"
    cat("Pas de métadonnées FolderID/FileID, sous-échantillonnage global uniquement.\n")
  }
  
  if (choice == "1") {
    if (is.null(n_max_cells)) {
      cat("Entrez le nombre maximum de cellules à conserver (défaut 200000) : ")
      user_input <- readline()
      if (user_input == "") {
        max_cells <- 200000
      } else {
        max_cells <- as.numeric(user_input)
        if (is.na(max_cells) || max_cells < 1000) {
          cat("⚠️ Valeur invalide. Utilisation de 200000 par défaut.\n")
          max_cells <- 200000
        }
      }
    } else {
      max_cells <- n_max_cells
    }
    expr_sub <- perform_global_subsampling(expr_data, max_cells)
    analysis_metrics$subsampling_method <- "global"
    
  } else if (choice == "2") {
    conditions_list <- detect_conditions_from_metadata(desc)
    cat("Entrez le nombre de cellules par condition (défaut = minimum disponible) : ")
    n_input <- readline()
    if (n_input == "") {
      n_target <- NULL
    } else {
      n_target <- as.numeric(n_input)
    }
    expr_sub <- perform_balanced_subsampling(expr_data, conditions_list, mode = "condition", n_cells_target = n_target)
    analysis_metrics$subsampling_method <- "balanced_by_condition"
    
  } else if (choice == "3") {
    conditions_list <- detect_conditions_from_metadata(desc)
    cat("Entrez le nombre de cellules par fichier (défaut = minimum disponible) : ")
    n_input <- readline()
    if (n_input == "") {
      n_target <- NULL
    } else {
      n_target <- as.numeric(n_input)
    }
    expr_sub <- perform_balanced_subsampling(expr_data, conditions_list, mode = "file", n_cells_target = n_target)
    analysis_metrics$subsampling_method <- "balanced_by_file"
    
  } else {
    stop("Choix invalide.")
  }
  
  analysis_metrics$subsampled_cell_count <- nrow(expr_sub)
  cat("✅ Sous-échantillonnage terminé (", analysis_metrics$subsampled_cell_count, "cellules pour l'analyse)\n")
  
  cat("📊 Génération des histogrammes de transformation pour les cellules sélectionnées...\n")
  selected_hist_file <- file.path(output_dir, "transformation_histograms_selected.pdf")
  show_fluorescence_histograms(expr_sub, channels_to_transform,
                               marker_descriptions, transformation_params,
                               filename = selected_hist_file,
                               title_suffix = "(cellules sélectionnées)", sample_name = sample_name)
  if (file.exists(selected_hist_file)) open_pdf(selected_hist_file)
  cat("✅ PDF des histogrammes sélectionnés généré :", selected_hist_file, "\n")
  
  markers_to_exclude <- c("FSC.A", "FSC.H", "FSC.W", "SSC.H", "SSC.W", "SSC.B.A", "SSC.B.H")
  markers_to_exclude <- c(markers_to_exclude, channels_to_exclude)
  markers_of_interest <- setdiff(channels_to_transform, markers_to_exclude)
  
  cat("\n🔍 Étape 7b/12 (Diagnostic) : Visualisation UMAP pour détecter les effets batch ...\n")
  
  if ("FolderID" %in% colnames(expr_sub) && "FileID" %in% colnames(expr_sub)) {
    
    mat_diag <- expr_sub[, markers_of_interest, drop = FALSE]
    mat_diag_scaled <- scale(mat_diag)
    
    set.seed(123)
    umap_diag <- uwot::umap(mat_diag_scaled,
                            n_neighbors = 25, min_dist = 0.25,
                            metric = "cosine", n_epochs = 200,
                            verbose = FALSE)
    
    diag_df <- data.frame(
      UMAP1 = umap_diag[, 1],
      UMAP2 = umap_diag[, 2],
      FolderID = expr_sub[, "FolderID"],
      FileID   = expr_sub[, "FileID"],
      stringsAsFactors = FALSE
    )
    
    conditions <- detect_conditions_from_metadata(desc)
    diag_df$Condition <- sapply(diag_df$FolderID,
                                function(id) get_condition_name(id, conditions))
    
    file_mapping <- list()
    if (!is.null(desc)) {
      file_keys <- names(desc)[grep("^File[0-9]+_Name$", names(desc))]
      if (length(file_keys) > 0) {
        for (key in file_keys) {
          file_id <- as.numeric(gsub("File([0-9]+)_Name", "\\1", key))
          file_name <- gsub("\\.fcs$", "", desc[[key]], ignore.case = TRUE)
          file_mapping[[as.character(file_id)]] <- file_name
        }
      }
    }
    diag_df$File <- sapply(diag_df$FileID, function(id) {
      id_char <- as.character(id)
      if (!is.na(id) && id_char %in% names(file_mapping)) {
        return(file_mapping[[id_char]])
      } else {
        return(paste0("File_", id_char))
      }
    })
    
    p_cond <- ggplot(diag_df, aes(x = UMAP1, y = UMAP2, color = Condition)) +
      geom_point(alpha = 0.2, size = 0.5) +
      theme_minimal(base_size = 12) +
      ggtitle(paste("UMAP diagnostique coloré par Condition -", sample_name)) +
      xlab("UMAP1") + ylab("UMAP2") +
      scale_color_discrete(name = "Condition") +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 8))
    
    n_files <- length(unique(diag_df$File))
    
    p_file <- ggplot(diag_df, aes(x = UMAP1, y = UMAP2, color = File)) +
      geom_point(alpha = 0.2, size = 0.5) +
      theme_minimal(base_size = 12) +
      ggtitle(paste("UMAP diagnostique coloré par Fichier -", sample_name)) +
      xlab("UMAP1") + ylab("UMAP2") +
      {if(n_files <= 8) scale_color_brewer(palette = "Set1", name = "Fichier")
        else scale_color_viridis_d(name = "Fichier")} +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
            legend.position = "bottom",
            legend.text = element_text(size = 6),
            legend.key.size = unit(0.2, "cm"))
    
    diag_pdf <- file.path(output_dir, "Batch_Diagnostic_UMAP.pdf")
    CairoPDF(diag_pdf, width = 14, height = 7)
    tryCatch({
      if (requireNamespace("ggpubr", quietly = TRUE)) {
        combined <- ggpubr::ggarrange(p_cond, p_file, ncol = 2)
        print(combined)
      } else {
        print(p_cond)
        print(p_file)
      }
    }, error = function(e) {
      print(p_cond)
      print(p_file)
    })
    dev.off()
    cat("✅ Graphique de diagnostic batch sauvegardé :", diag_pdf, "\n")
    if (interactive()) open_pdf(diag_pdf)
    
    cat("\n--- DIAGNOSTIC DES EFFETS BATCH ---\n")
    cat("Veuillez examiner le PDF 'Batch_Diagnostic_UMAP.pdf' dans le dossier de sortie.\n")
    cat("Une séparation nette par condition ou par fichier est-elle visible ?\n")
    cat("Si oui, un effet batch est probable et la normalisation (CytoNorm) est recommandée.\n")
    cat("Voulez-vous normaliser les batch effects ? (y/n) : ")
    do_batch <- tolower(readline())
    
    if (do_batch %in% c("y", "yes")) {
      expr_sub <- check_batch_effects_and_normalize(
        expr_data = expr_sub,
        file_ids = expr_sub[, "FileID"],
        folder_ids = expr_sub[, "FolderID"],
        channels = markers_of_interest,
        output_dir = output_dir,
        sample_name = sample_name,
        transformation_params = transformation_params,
        desc = desc,
        marker_descriptions = marker_descriptions
      )
      analysis_metrics$batch_normalized <- TRUE
    } else {
      analysis_metrics$batch_normalized <- FALSE
      cat("Normalisation des batch effects ignorée.\n")
    }
    
  } else {
    cat("Colonnes FolderID et/ou FileID absentes. Normalisation des batch effects non applicable.\n")
    analysis_metrics$batch_normalized <- FALSE
  }
  
  cat("🤖️ Étape 8/12: Analyse FlowSOM en cours...\n")
  flowSOM_res <- FlowSOM(input = expr_sub[, markers_of_interest],
                         seed = analysis_metrics$seed,
                         compensate = FALSE,
                         transform = FALSE,
                         scale = FALSE,
                         colsToUse = markers_of_interest,
                         xdim = analysis_metrics$flowsom_xdim,
                         ydim = analysis_metrics$flowsom_ydim,
                         rlen = analysis_metrics$flowsom_rlen)
  cat("✅ Analyse FlowSOM terminée\n")
  
  cat("💻 Étape 9/12: Détermination robuste du nombre optimal de clusters...\n")
  
  codes <- flowSOM_res$map$codes
  n_codes <- nrow(codes)
  n_clusters_max <- min(50, n_codes - 1)
  
  dist_codes <- dist(codes)
  hc <- hclust(dist_codes, method = "ward.D2")
  
  cat("🔍 Application de Dynamic Tree Cut (méthode hybrid)...\n")
  dynamic_clusters <- tryCatch(
    cutreeDynamic(
      dendro = hc,
      distM = as.matrix(dist_codes),   
      method = "hybrid",
      deepSplit = 2,
      minClusterSize = max(5, round(n_codes * 0.05)),
      pamRespectsDendro = TRUE,
      verbose = 2
    ),
    error = function(e) {
      message("Dynamic Tree Cut a échoué, repli sur un seul cluster.")
      rep(1, n_codes)
    }
  )
  k_dynamic <- length(unique(dynamic_clusters[dynamic_clusters > 0]))
  k_dynamic <- max(2, min(k_dynamic, n_clusters_max))
  
  sil_scores <- numeric(n_clusters_max - 1)
  for (k in 2:n_clusters_max) {
    km <- kmeans(codes, centers = k, nstart = 25, iter.max = 100)
    ss <- silhouette(km$cluster, dist_codes)
    sil_scores[k - 1] <- mean(ss[, 3])
  }
  optimal_k_sil <- which.max(sil_scores) + 1
  
  gap_stat <- NULL
  optimal_k_gap <- NA
  tryCatch({
    gap_stat <- clusGap(
      codes,
      FUNcluster = kmeans,
      nstart = 25,
      K.max = n_clusters_max,
      B = 50,
      verbose = FALSE,
      iter.max = 100
    )
    optimal_k_gap <- maxSE(
      gap_stat$Tab[, "gap"],
      gap_stat$Tab[, "SE.sim"],
      method = "Tibs2001SEmax"
    )
  }, error = function(e) {
    message("Gap statistic non disponible : ", e$message)
  })
  
  cat("\n🔍 Calcul du taux de stabilité des clusters pour k = 2 à", n_clusters_max, "...\n")
  cat("   (Bootstrap avec B = 200, seuil de stabilité = 0.75)\n")
  cat("   Cela peut prendre quelques instants...\n")
  
  stable_rate_scores    <- numeric(n_clusters_max - 1)
  stable_rate_scores_60 <- numeric(n_clusters_max - 1)   # AJOUT : seuil ≥0.60
  jaccard_mean_scores   <- numeric(n_clusters_max - 1)
  names(stable_rate_scores)    <- 2:n_clusters_max
  names(stable_rate_scores_60) <- 2:n_clusters_max       # AJOUT
  names(jaccard_mean_scores)   <- 2:n_clusters_max
  
  for (k_test in 2:n_clusters_max) {
    set.seed(123)
    codes_df <- as.data.frame(codes)
    msg_con <- file(tempfile(), open = "wt")
    sink(msg_con, type = "message")
    stab_test <- tryCatch(
      clusterboot(
        codes_df,
        B = 200,
        clustermethod = hclustCBI,
        method = "ward.D2",
        k = k_test,
        seed = 123,
        scaling = FALSE,
        iter.max = 200
      ),
      finally = {
        sink(type = "message")
        close(msg_con)
      }
    )
    try(sink(type = "message"), silent = TRUE)
    
    jaccard_per_cluster <- stab_test$bootmean
    if (is.null(jaccard_per_cluster)) {
      jaccard_per_cluster <- stab_test$bootbrd
    }
    
    if (!is.null(jaccard_per_cluster)) {
      stable_rate    <- mean(jaccard_per_cluster >= 0.75, na.rm = TRUE)
      stable_rate_60 <- mean(jaccard_per_cluster >= 0.60, na.rm = TRUE)  # AJOUT
      mean_jaccard   <- mean(jaccard_per_cluster, na.rm = TRUE)
    } else {
      stable_rate    <- NA
      stable_rate_60 <- NA                                                 # AJOUT
      mean_jaccard   <- NA
    }
    
    stable_rate_scores[as.character(k_test)]    <- stable_rate
    stable_rate_scores_60[as.character(k_test)] <- stable_rate_60         # AJOUT
    jaccard_mean_scores[as.character(k_test)]   <- mean_jaccard
    
    cat(sprintf("   k = %d : Taux stable = %.1f%% | Jaccard moyen = %.3f\n", 
                k_test, stable_rate * 100, mean_jaccard))
  }
  
  k_candidates <- which(stable_rate_scores >= 0.75) + 1
  if (length(k_candidates) > 0) {
    optimal_k_stable <- max(k_candidates)
  } else {
    k_candidates <- which(stable_rate_scores >= 0.60) + 1
    if (length(k_candidates) > 0) {
      optimal_k_stable <- max(k_candidates)
      cat("⚠️ Aucun k n'atteint 75% de stabilité. Sélection basée sur 60%.\n")
    } else {
      optimal_k_stable <- which.max(stable_rate_scores) + 1
      cat("⚠️ Stabilité faible globalement. Sélection du meilleur k disponible.\n")
    }
  }
  
  cluster_selection_file <- file.path(output_dir, "Cluster_Selection_Methods.pdf")
  CairoPDF(cluster_selection_file, width = 14, height = 10)
  par(mfrow = c(2, 3))
  
  plot(hc,
       main = "Dendrogramme (Ward) avec coupe Dynamic Tree Cut",
       xlab = "Codes SOM", sub = "", cex = 0.6)
  dynamic_colors <- rainbow(k_dynamic)[dynamic_clusters]
  rect.hclust(hc, k = k_dynamic, border = "red")
  
  plot(2:n_clusters_max, sil_scores, type = "b", pch = 19,
       xlab = "Nombre de clusters (k)", ylab = "Score de silhouette moyen",
       main = paste("Silhouette - maximum à k =", optimal_k_sil))
  points(optimal_k_sil, sil_scores[optimal_k_sil - 1], col = "red", pch = 19, cex = 2)
  abline(v = k_dynamic, col = "blue", lty = 2, lwd = 2)
  text(k_dynamic, max(sil_scores) * 0.9,
       labels = paste("Dynamic k =", k_dynamic), col = "blue")
  
  if (!is.null(gap_stat)) {
    gap_values <- gap_stat$Tab[, "gap"]
    se_values  <- gap_stat$Tab[, "SE.sim"]
    plot(1:n_clusters_max, gap_values, type = "b", pch = 19,
         xlab = "Nombre de clusters (k)", ylab = "Gap statistic",
         main = paste("Gap statistic - k choisi =", optimal_k_gap))
    if (!is.na(optimal_k_gap)) {
      points(optimal_k_gap, gap_values[optimal_k_gap], col = "red", pch = 19, cex = 2)
    }
    arrows(1:n_clusters_max, gap_values - se_values,
           1:n_clusters_max, gap_values + se_values,
           length = 0.05, angle = 90, code = 3)
  } else {
    plot.new()
    text(0.5, 0.5, "Gap statistic non disponible", cex = 1.2)
  }
  
  # --- Graphique taux de stabilité (modifié) ---
  k_axis <- 2:n_clusters_max
  
  # Lignes verticales semi-transparentes tous les 10 clusters
  v_lines <- seq(10, n_clusters_max, by = 10)
  v_lines <- v_lines[v_lines >= 2 & v_lines <= n_clusters_max]
  
  plot(k_axis, stable_rate_scores * 100,
       type = "b", pch = 19, col = "darkgreen",
       xlab = "Nombre de clusters (k)", ylab = "% de clusters stables",
       main = paste("Taux de stabilité - k suggéré =", optimal_k_stable),
       ylim = c(0, 105))
  
  # Lignes verticales semi-transparentes tous les 10 clusters
  if (length(v_lines) > 0) {
    for (vl in v_lines) {
      abline(v = vl, col = adjustcolor("grey40", alpha.f = 0.35), lty = 1, lwd = 1)
    }
  }
  
  # Lignes horizontales de seuil
  abline(h = 75, col = "darkgreen", lty = 2, lwd = 2)
  abline(h = 60, col = "darkorange", lty = 2, lwd = 2)
  
  # Seconde courbe : seuil ≥0.60
  lines(k_axis, stable_rate_scores_60 * 100,
        type = "b", pch = 17, col = "darkorange", lwd = 1.5)
  
  # Point optimal
  points(optimal_k_stable,
         stable_rate_scores[as.character(optimal_k_stable)] * 100,
         col = "red", pch = 19, cex = 2)
  
  # Légende
  legend("topright",
         legend = c("Jaccard ≥ 0.75",
                    "Jaccard ≥ 0.60)",
                    "Seuil 75%",
                    "Seuil 60%",
                    "k sélectionné"),
         col    = c("darkgreen", "darkorange", "darkgreen", "darkorange", "red"),
         lty    = c(1, 1, 2, 2, NA),
         pch    = c(19, 17, NA, NA, 19),
         lwd    = c(1.5, 1.5, 2, 2, NA),
         cex    = 0.8,
         bg     = "white")
  
  plot.new()
  text(0.5, 0.5,
       labels = paste(
         "Résumé des méthodes :\n\n",
         "Dynamic Tree Cut : k =", k_dynamic, "\n",
         "Silhouette maximale : k =", optimal_k_sil, "\n",
         if (!is.null(gap_stat)) paste("Gap statistic : k =", optimal_k_gap) else "", "\n",
         "Taux stabilité ≥75% : k =", optimal_k_stable, "\n"
       ),
       cex = 0.9, adj = c(0.5, 0.5))
  
  dev.off()
  cat("✅ Graphiques de sélection sauvegardés :", cluster_selection_file, "\n")
  if (interactive()) open_pdf(cluster_selection_file)
  
  cat("\n--- Méthodes de sélection automatique ---\n")
  cat("Dynamic Tree Cut (hybride, deepSplit=2) : k =", k_dynamic, "\n")
  cat("Silhouette moyenne (max)                : k =", optimal_k_sil, "\n")
  if (!is.null(gap_stat)) {
    cat("Gap statistic (Tibs2001SEmax)          : k =", optimal_k_gap, "\n")
  }
  cat("Taux stabilité ≥75% (suggéré)            : k =", optimal_k_stable, "\n")
  cat("Jaccard moyen max (référence)            : k =", which.max(jaccard_mean_scores) + 1, "\n")
  
  cat("\nVeuillez examiner le PDF 'Cluster_Selection_Methods.pdf'.\n")
  cat("Le taux de stabilité indique quelle proportion de vos clusters\n")
  cat("sont des populations cellulaires réellement reproductibles.\n\n")
  
  repeat {
    cat("Entrez le nombre de clusters souhaité (Entrée = accepter",
        optimal_k_stable, ", q = quitter) : ")
    user_input <- readline()
    if (tolower(user_input) == "q") {
      stop("Analyse interrompue par l'utilisateur.")
    } else if (user_input == "") {
      optimal_k <- optimal_k_stable
      break
    } else {
      user_k <- suppressWarnings(as.numeric(user_input))
      if (!is.na(user_k) && user_k >= 2 && user_k <= n_clusters_max) {
        optimal_k <- user_k
        break
      } else {
        cat("Valeur invalide. Veuillez entrer un nombre entre 2 et", n_clusters_max, "\n")
      }
    }
  }
  analysis_metrics$optimal_k <- optimal_k
  cat("✅ Nombre de clusters initial :", optimal_k, "\n")
  
  current_k <- optimal_k
  stability_result <- NULL
  
  repeat {
    cat("\n🧪 Évaluation détaillée de la stabilité pour k =", current_k, "(bootstrap, B = 600)...\n")
    set.seed(123)
    codes_df <- as.data.frame(codes)
    msg_con <- file(tempfile(), open = "wt")
    sink(msg_con, type = "message")
    stab <- tryCatch(
      clusterboot(
        codes_df,
        B = 600,
        clustermethod = hclustCBI,
        method = "ward.D2",
        k = current_k,
        seed = 123,
        scaling = FALSE,
        iter.max = 200
      ),
      finally = {
        sink(type = "message")
        close(msg_con)
      }
    )
    try(sink(type = "message"), silent = TRUE)
    
    jaccard_per_cluster <- stab$bootmean
    n_stable <- sum(jaccard_per_cluster >= 0.75, na.rm = TRUE)
    n_acceptable <- sum(jaccard_per_cluster >= 0.6, na.rm = TRUE)
    global_rate <- mean(jaccard_per_cluster >= 0.75, na.rm = TRUE)
    
    cat("\n--- Stabilité détaillée des clusters (Jaccard, B=500) ---\n")
    cat(sprintf("Seuil stable : ≥0.75 | Seuil acceptable : ≥0.60\n\n"))
    
    for (i in seq_along(jaccard_per_cluster)) {
      j <- jaccard_per_cluster[i]
      status <- if (j >= 0.75) "✅ Stable" else if (j >= 0.6) "⚠️ Acceptable" else "❌ Instable"
      cat(sprintf("Cluster %2d : Jaccard = %.3f %s\n", i, j, status))
    }
    
    cat(sprintf("\n--- Résumé ---\n"))
    cat(sprintf("Clusters stables (≥0.75)     : %d/%d (%.1f%%)\n", 
                n_stable, current_k, n_stable/current_k*100))
    cat(sprintf("Clusters acceptables (≥0.60) : %d/%d (%.1f%%)\n", 
                n_acceptable, current_k, n_acceptable/current_k*100))
    cat(sprintf("Jaccard moyen global         : %.3f\n", mean(jaccard_per_cluster)))
    
    stability_result <- stab
    
    if (global_rate < 0.5) {
      cat("\n🚨 ATTENTION : Moins de 50% de vos clusters sont stables.\n")
      cat("   Cela suggère un sur-clustering. Envisagez de diminuer k.\n")
    } else if (global_rate < 0.75) {
      cat("\n⚠️  Note : Seulement ", round(global_rate*100, 1), 
          "% des clusters sont stables. Certains clusters peuvent être artificiels.\n", sep="")
    }
    
    cat("\nVoulez-vous modifier le nombre de clusters ? (y/n) : ")
    answer <- tolower(readline())
    if (!(answer %in% c("y", "yes"))) break
    
    repeat {
      cat(sprintf("Entrez le nouveau nombre de clusters (entre 2 et %d) : ", n_clusters_max))
      new_k <- suppressWarnings(as.integer(readline()))
      if (!is.na(new_k) && new_k >= 2 && new_k <= n_clusters_max) {
        current_k <- new_k
        break
      }
      cat("Valeur invalide.\n")
    }
  }
  
  optimal_k <- current_k
  analysis_metrics$optimal_k <- optimal_k
  
  stab_file <- file.path(output_dir, "Cluster_Stability_Final.csv")
  stab_df <- data.frame(
    Cluster = seq_along(stability_result$bootmean),
    Jaccard = stability_result$bootmean,
    Stability = ifelse(stability_result$bootmean >= 0.75, "Stable",
                       ifelse(stability_result$bootmean >= 0.6, "Acceptable", "Instable"))
  )
  write.csv(stab_df, stab_file, row.names = FALSE)
  cat("📋 Résultats de stabilité détaillés sauvegardés dans", stab_file, "\n")
  cat("✅ Nombre de clusters retenu après validation :", optimal_k, "\n")
  
  cat("💻 Étape 10/12: Meta-clustering en cours...\n")
  metaclusters <- metaClustering_consensus(flowSOM_res$map$codes, k = optimal_k)
  cluster_ids <- metaclusters[flowSOM_res$map$mapping[,1]]
  
  if (length(cluster_ids) != nrow(expr_sub)) {
    stop("Erreur: cluster_ids et expr_sub n'ont pas la même longueur")
  }
  
  wrap_title <- function(title, max_chars = 20) {
    if (nchar(title) <= max_chars) return(title)
    words <- strsplit(title, " ")[[1]]
    line1 <- ""
    line2 <- ""
    for (w in words) {
      candidate <- if (nchar(line1) == 0) w else paste(line1, w)
      if (nchar(candidate) <= max_chars) {
        line1 <- candidate
      } else {
        line2 <- if (nchar(line2) == 0) w else paste(line2, w)
      }
    }
    if (nchar(line2) == 0) return(line1)
    paste0(line1, "\n", line2)
  }
  
  plot_fsc_ssc_by_cluster(expr_sub, cluster_ids, sample_name, output_dir)
  cat("✅ Meta-clustering terminé\n")
  
  cat("🔨 Étape 11/12: Annotation des clusters (DeepSeek)...\n")
  expr_sub_markers <- expr_sub[, markers_of_interest]
  unique_clusters <- sort(unique(cluster_ids))
  cluster_medians <- matrix(NA, nrow = length(unique_clusters), ncol = ncol(expr_sub_markers))
  rownames(cluster_medians) <- unique_clusters
  colnames(cluster_medians) <- colnames(expr_sub_markers)
  for (cluster in unique_clusters) {
    cluster_cells <- which(cluster_ids == cluster)
    if (length(cluster_cells) > 0) {
      cluster_medians[as.character(cluster), ] <- apply(expr_sub_markers[cluster_cells, , drop = FALSE],
                                                        2, median, na.rm = TRUE)
    }
  }
  
  annotation_results <- annotate_clusters_improved(
    cluster_medians = cluster_medians,
    marker_descriptions = marker_descriptions,
    output_dir = output_dir,
    population_description = population_description
  )
  cluster_annotations <- annotation_results$annotations
  cluster_labels <- paste0(1:length(cluster_annotations), "_", cluster_annotations)
  annotation_source <- annotation_results$source
  
  cluster_medians_data <- save_cluster_medians_with_thresholds(
    expr_sub, cluster_ids,
    setNames(rep(NA_real_, length(markers_of_interest)), markers_of_interest),
    marker_descriptions, cluster_labels, output_dir, sample_name
  )
  
  detailed_report_file <- generate_detailed_annotation_report(
    annotation_results, cluster_medians, marker_descriptions, output_dir, sample_name
  )
  cat("✅ Annotation des clusters terminée\n")
  cat("✅ Rapport d'annotation détaillé généré :", detailed_report_file, "\n")
  
  cat("📊 Génération du PDF FSC/SSC avec les clusters et leurs noms finaux...\n")
  fsc_ssc_filtered_file <- plot_fsc_ssc_by_cluster_filtered(
    expr_sub, cluster_ids, cluster_labels, sample_name, output_dir
  )
  cat("✅ PDF FSC/SSC généré :", fsc_ssc_filtered_file, "\n")
  
  cat("⏳ Étape 12/12: Génération des visualisations (UMAP et heatmap)...\n")
  expr_sub_markers_for_viz <- expr_sub[, markers_of_interest]
  write.csv(expr_sub_markers_for_viz, file.path(output_dir, "expression_data_for_viz.csv"), row.names = FALSE)
  
  cat("📊 Standardisation des données pour UMAP...\n")
  expr_sub_scaled <- scale(expr_sub_markers_for_viz)
  
  dimred_cache_file <- file.path(output_dir, "dimred_cache.rds")
  cached_dimred <- NULL
  if (file.exists(dimred_cache_file)) {
    cached_dimred <- tryCatch({
      tmp <- readRDS(dimred_cache_file)
      if (is.list(tmp) && all(c("umap", "data_hash") %in% names(tmp))) tmp
      else NULL
    }, error = function(e) NULL)
    if (!is.null(cached_dimred)) {
      cat("📦 Cache UMAP chargé depuis le disque.\n")
    } else {
      cat("⚠️ Cache existant corrompu, sera recréé.\n")
    }
  }
  
  current_data_hash <- digest::digest(expr_sub_markers_for_viz)
  if (!is.null(cached_dimred) && identical(cached_dimred$data_hash, current_data_hash)) {
    umap_res <- cached_dimred$umap
    cat("♻️ Réutilisation du cache UMAP\n")
  } else {
    cat("📊 Calcul de l'UMAP...\n")
    set.seed(123)
    umap_res <- uwot::umap(expr_sub_scaled,
                           n_neighbors = analysis_metrics$umap_n_neighbors,
                           min_dist = analysis_metrics$umap_min_dist,
                           metric = analysis_metrics$umap_metric,
                           n_epochs = analysis_metrics$umap_n_epochs,
                           n_threads = analysis_metrics$umap_n_threads,
                           verbose = TRUE)
    cached_dimred <- list(umap = umap_res, data_hash = current_data_hash)
    saveRDS(cached_dimred, dimred_cache_file)
    cat("💾 Cache UMAP sauvegardé.\n")
  }
  
  umap_df <- data.frame(UMAP1 = umap_res[, 1],
                        UMAP2 = umap_res[, 2],
                        Cluster = factor(cluster_labels[cluster_ids], levels = cluster_labels))
  
  cluster_medians_pre <- apply(expr_sub_markers_for_viz, 2, function(x) {
    tapply(x, cluster_ids, median, na.rm = TRUE)
  })
  rownames(cluster_medians_pre) <- cluster_labels[as.numeric(rownames(cluster_medians_pre))]
  marker_names_for_heatmap <- sapply(colnames(cluster_medians_pre), function(ch) {
    get_marker_name(ch, marker_descriptions)
  })
  colnames(cluster_medians_pre) <- marker_names_for_heatmap
  
  graphics.off()
  heatmap_file <- file.path(output_dir, "Cluster_Heatmap.pdf")
  cat("📊 Création du heatmap initial...\n")
  CairoPDF(heatmap_file, width = 10, height = 8)
  hm <- Heatmap(cluster_medians_pre,
                name = "Expression",
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                col = viridis(256),
                row_title = "Clusters",
                column_title = paste("Marker Expression -", sample_name),
                row_names_gp = gpar(fontsize = 5),
                column_names_gp = gpar(fontsize = 8))
  print(hm)
  dev.off()
  cat("✅ Heatmap initial créé avec succès\n")
  
  umap_files <- create_enhanced_umap_plot(umap_df, cluster_labels, sample_name, output_dir)
  
  open_pdf(heatmap_file)
  open_pdf(umap_files)
  
  cat("\n=== RÉSUMÉ DES ANNOTATIONS AUTOMATIQUES ===\n")
  annotation_summary <- table(cluster_annotations)
  for (annotation in names(annotation_summary)) {
    cat(paste(annotation, ":", annotation_summary[annotation], "clusters\n"))
  }
  cat("\n")
  
  manual_rename_done <- FALSE
  if (annotation_source != "AI") {
    cat("\n⚠️  L'annotation automatique n'a pas pu être réalisée.")
    cat("\n→ Lancement du renommage manuel des clusters (les heatmap et UMAP sont déjà affichés).\n")
    cluster_summary_for_rename <- data.frame(
      Cluster = 1:length(cluster_annotations),
      Annotation = cluster_annotations,
      Count = as.numeric(table(cluster_ids)),
      Percentage = round(100 * as.numeric(table(cluster_ids)) / sum(table(cluster_ids)), 2)
    )
    cluster_annotations <- rename_clusters_manually(cluster_annotations, cluster_summary_for_rename, sample_name)
    cluster_labels <- paste0(1:length(cluster_annotations), "_", cluster_annotations)
    if (!is.null(annotation_results)) annotation_results$annotations <- cluster_annotations
    manual_rename_done <- TRUE
  }
  
  cat("\n=== FILTRAGE DES CLUSTERS ===\n")
  cat("Voulez-vous éliminer des clusters après visualisation du heatmap et des visualisations? (y/n): ")
  filter_choice <- readline()
  if (tolower(filter_choice) %in% c("y", "yes")) {
    filter_result <- filter_clusters_interactive(
      expr_sub = expr_sub,
      cluster_ids = cluster_ids,
      cluster_annotations = cluster_annotations,
      markers_of_interest = markers_of_interest,
      marker_descriptions = marker_descriptions,
      sample_name = sample_name,
      output_dir = output_dir,
      umap_df = umap_df,
      cluster_labels = cluster_labels,
      annotation_results = annotation_results
    )
    if (!is.null(filter_result$filter_info)) {
      expr_sub <- filter_result$expr_sub
      cluster_ids <- filter_result$cluster_ids
      cluster_annotations <- filter_result$cluster_annotations
      cluster_labels <- filter_result$cluster_labels
      umap_df <- filter_result$umap_df
      annotation_results <- filter_result$annotation_results
      
      expr_sub_markers_for_viz <- expr_sub[, markers_of_interest]
      cluster_medians_pre <- apply(expr_sub_markers_for_viz, 2, function(x) {
        tapply(x, cluster_ids, median, na.rm = TRUE)
      })
      rownames(cluster_medians_pre) <- cluster_labels[as.numeric(rownames(cluster_medians_pre))]
      colnames(cluster_medians_pre) <- marker_names_for_heatmap
      
      graphics.off()
      heatmap_file_filtered <- file.path(output_dir, "Cluster_Heatmap_Filtered.pdf")
      CairoPDF(heatmap_file_filtered, width = 10, height = 8)
      hm_filtered <- Heatmap(cluster_medians_pre,
                             name = "Expression",
                             cluster_rows = TRUE,
                             cluster_columns = TRUE,
                             show_row_names = TRUE,
                             show_column_names = TRUE,
                             col = viridis(256),
                             row_title = "Clusters",
                             column_title = paste("Marker Expression -", sample_name, "(filtré)"),
                             row_names_gp = gpar(fontsize = 5),
                             column_names_gp = gpar(fontsize = 8))
      print(hm_filtered)
      dev.off()
      
      expr_sub_scaled <- scale(expr_sub_markers_for_viz)
      current_data_hash <- digest::digest(expr_sub_markers_for_viz)
      if (!is.null(cached_dimred) && identical(cached_dimred$data_hash, current_data_hash)) {
        umap_res <- cached_dimred$umap
        cat("♻️ Réutilisation du cache UMAP (post-filtrage).\n")
      } else {
        set.seed(123)
        umap_res <- uwot::umap(expr_sub_scaled,
                               n_neighbors = analysis_metrics$umap_n_neighbors,
                               min_dist = analysis_metrics$umap_min_dist,
                               metric = analysis_metrics$umap_metric,
                               n_epochs = analysis_metrics$umap_n_epochs,
                               n_threads = analysis_metrics$umap_n_threads,
                               verbose = TRUE)
        cached_dimred <- list(umap = umap_res, data_hash = current_data_hash)
        saveRDS(cached_dimred, dimred_cache_file)
        cat("💾 Nouveau cache UMAP sauvegardé (post-filtrage).\n")
      }
      
      umap_df <- data.frame(UMAP1 = umap_res[, 1],
                            UMAP2 = umap_res[, 2],
                            Cluster = factor(cluster_labels[cluster_ids], levels = cluster_labels))
      
      umap_files_filtered <- create_enhanced_umap_plot(umap_df, cluster_labels, sample_name, output_dir, suffix = "_Filtered")
      
      open_pdf(heatmap_file_filtered)
      open_pdf(umap_files_filtered)
      if (!is.null(filter_result$fsc_ssc_filtered_file)) open_pdf(filter_result$fsc_ssc_filtered_file)
      
      heatmap_file <- heatmap_file_filtered
      umap_files <- umap_files_filtered
      
      cat("\n✅ Filtrage des clusters terminé avec succès!\n")
    }
  }
  
  cat("\n🔄 Étape 12b/12: Fusion interactive des clusters...\n")
  merge_counter <- 0
  continue_merge <- TRUE
  while (continue_merge) {
    cat("\nVoulez-vous fusionner des clusters? (y/n): ")
    merge_choice <- readline()
    if (tolower(merge_choice) %in% c("y", "yes")) {
      merge_counter <- merge_counter + 1
      suffix <- if (merge_counter > 1) paste0("_merged", merge_counter) else "_merged1"
      cat(paste0("\n🔄 Cycle de fusion ", merge_counter, " en cours...\n"))
      
      merge_result <- merge_clusters_interactive(
        expr_sub = expr_sub,
        cluster_ids = cluster_ids,
        cluster_annotations = cluster_annotations,
        markers_of_interest = markers_of_interest,
        marker_descriptions = marker_descriptions,
        sample_name = sample_name,
        output_dir = output_dir,
        umap_df = umap_df,
        cluster_labels = cluster_labels,
        annotation_results = annotation_results,
        suffix = suffix
      )
      
      if (!is.null(merge_result$merge_info)) {
        cluster_ids <- merge_result$cluster_ids
        cluster_annotations <- merge_result$cluster_annotations
        cluster_labels <- merge_result$cluster_labels
        umap_df <- merge_result$umap_df
        annotation_results <- merge_result$annotation_results
        
        expr_sub_markers_for_viz <- expr_sub[, markers_of_interest]
        cluster_medians_pre <- apply(expr_sub_markers_for_viz, 2, function(x) {
          tapply(x, cluster_ids, median, na.rm = TRUE)
        })
        rownames(cluster_medians_pre) <- cluster_labels[as.numeric(rownames(cluster_medians_pre))]
        colnames(cluster_medians_pre) <- marker_names_for_heatmap
        
        graphics.off()
        heatmap_file_merged <- file.path(output_dir, paste0("Cluster_Heatmap", suffix, ".pdf"))
        CairoPDF(heatmap_file_merged, width = 10, height = 8)
        hm_merged <- Heatmap(cluster_medians_pre,
                             name = "Expression",
                             cluster_rows = TRUE,
                             cluster_columns = TRUE,
                             show_row_names = TRUE,
                             show_column_names = TRUE,
                             col = viridis(256),
                             row_title = "Clusters",
                             column_title = paste("Marker Expression -", sample_name, "(fusion", merge_counter, ")"),
                             row_names_gp = gpar(fontsize = 5),
                             column_names_gp = gpar(fontsize = 8))
        print(hm_merged)
        dev.off()
        
        expr_sub_scaled <- scale(expr_sub_markers_for_viz)
        current_data_hash <- digest::digest(expr_sub_markers_for_viz)
        if (!is.null(cached_dimred) && identical(cached_dimred$data_hash, current_data_hash)) {
          umap_res <- cached_dimred$umap
          cat("♻️ Réutilisation du cache UMAP (post-fusion).\n")
        } else {
          set.seed(123)
          umap_res <- uwot::umap(expr_sub_scaled,
                                 n_neighbors = analysis_metrics$umap_n_neighbors,
                                 min_dist = analysis_metrics$umap_min_dist,
                                 metric = analysis_metrics$umap_metric,
                                 n_epochs = analysis_metrics$umap_n_epochs,
                                 n_threads = analysis_metrics$umap_n_threads,
                                 verbose = TRUE)
          cached_dimred <- list(umap = umap_res, data_hash = current_data_hash)
          saveRDS(cached_dimred, dimred_cache_file)
          cat("💾 Cache UMAP sauvegardé (post-fusion).\n")
        }
        
        umap_df <- data.frame(UMAP1 = umap_res[, 1],
                              UMAP2 = umap_res[, 2],
                              Cluster = factor(cluster_labels[cluster_ids], levels = cluster_labels))
        
        umap_files_merged <- create_enhanced_umap_plot(umap_df, cluster_labels, sample_name, output_dir, suffix = suffix)
        
        open_pdf(heatmap_file_merged)
        open_pdf(umap_files_merged)
        if (!is.null(merge_result$fsc_ssc_filtered_file)) open_pdf(merge_result$fsc_ssc_filtered_file)
        
        cat("\n✅ Fusion ", merge_counter, " terminée avec succès!\n")
        
        cat("\nVoulez-vous encore fusionner des clusters? (y/n): ")
        continue_input <- readline()
        continue_merge <- tolower(continue_input) %in% c("y", "yes")
      } else {
        cat("❌ Aucune fusion effectuée.\n")
        continue_merge <- FALSE
      }
    } else {
      continue_merge <- FALSE
    }
  }
  cat("✅ Fusion interactive terminée\n")
  
  cat("🧰 Étape 12c/12: Proposition de renommage manuel...\n")
  if (!manual_rename_done) {
    cat("Voulez-vous renommer les clusters après visualisation du heatmap et des visualisations? (o/n): ")
    rename_after_viz <- readline()
  } else {
    cat("Un premier renommage a déjà été effectué. Voulez-vous modifier à nouveau les noms? (o/n): ")
    rename_after_viz <- readline()
  }
  if (tolower(rename_after_viz) %in% c("o", "oui", "y", "yes")) {
    cluster_summary_for_rename <- data.frame(
      Cluster = 1:length(cluster_annotations),
      Annotation = cluster_annotations,
      Count = as.numeric(table(cluster_ids)),
      Percentage = round(100 * as.numeric(table(cluster_ids)) / sum(table(cluster_ids)), 2),
      Confidence_Score = annotation_results$scores
    )
    new_annotations <- rename_clusters_manually(cluster_annotations, cluster_summary_for_rename, sample_name)
    cluster_annotations <- new_annotations
    cluster_labels <- paste0(1:length(cluster_annotations), "_", cluster_annotations)
    if (!is.null(annotation_results)) annotation_results$annotations <- cluster_annotations
    
    expr_sub_markers_for_viz <- expr_sub[, markers_of_interest]
    cluster_medians_pre_updated <- apply(expr_sub_markers_for_viz, 2, function(x) {
      tapply(x, cluster_ids, median, na.rm = TRUE)
    })
    rownames(cluster_medians_pre_updated) <- cluster_labels[as.numeric(rownames(cluster_medians_pre_updated))]
    colnames(cluster_medians_pre_updated) <- marker_names_for_heatmap
    
    graphics.off()
    heatmap_file_updated <- file.path(output_dir, "Cluster_Heatmap_Updated.pdf")
    CairoPDF(heatmap_file_updated, width = 10, height = 8)
    print(Heatmap(cluster_medians_pre_updated,
                  name = "Expression",
                  cluster_rows = TRUE,
                  cluster_columns = TRUE,
                  show_row_names = TRUE,
                  show_column_names = TRUE,
                  col = viridis(256),
                  row_title = "Clusters",
                  column_title = paste("Marker Expression -", sample_name),
                  row_names_gp = gpar(fontsize = 5),
                  column_names_gp = gpar(fontsize = 8)))
    dev.off()
    
    expr_sub_scaled <- scale(expr_sub_markers_for_viz)
    current_data_hash <- digest::digest(expr_sub_markers_for_viz)
    if (!is.null(cached_dimred) && identical(cached_dimred$data_hash, current_data_hash)) {
      umap_res <- cached_dimred$umap
      cat("♻️ Réutilisation du cache UMAP (post-renommage).\n")
    } else {
      set.seed(123)
      umap_res <- uwot::umap(expr_sub_scaled,
                             n_neighbors = analysis_metrics$umap_n_neighbors,
                             min_dist = analysis_metrics$umap_min_dist,
                             metric = analysis_metrics$umap_metric,
                             n_epochs = analysis_metrics$umap_n_epochs,
                             n_threads = analysis_metrics$umap_n_threads,
                             verbose = TRUE)
      cached_dimred <- list(umap = umap_res, data_hash = current_data_hash)
      saveRDS(cached_dimred, dimred_cache_file)
      cat("💾 Cache UMAP sauvegardé (post-renommage).\n")
    }
    
    umap_df <- data.frame(UMAP1 = umap_res[, 1],
                          UMAP2 = umap_res[, 2],
                          Cluster = factor(cluster_labels[cluster_ids], levels = cluster_labels))
    
    umap_files_updated <- create_enhanced_umap_plot(umap_df, cluster_labels, sample_name, output_dir, suffix = "_Updated")
    
    open_pdf(heatmap_file_updated)
    open_pdf(umap_files_updated)
    
    cat("📊 Régénération du PDF FSC/SSC avec les nouveaux noms de clusters...\n")
    fsc_ssc_filtered_file <- plot_fsc_ssc_by_cluster_filtered(
      expr_sub, cluster_ids, cluster_labels, sample_name, output_dir
    )
    open_pdf(fsc_ssc_filtered_file)
    cat("✅ PDF FSC/SSC mis à jour :", fsc_ssc_filtered_file, "\n")
    
    heatmap_file <- heatmap_file_updated
    umap_files <- umap_files_updated
    
    cluster_medians_data <- save_cluster_medians_with_thresholds(
      expr_sub, cluster_ids,
      setNames(rep(NA_real_, length(markers_of_interest)), markers_of_interest),
      marker_descriptions, cluster_labels, output_dir, sample_name
    )
  }
  
  annotation_mode <- if (annotation_source == "AI") "Automatique" else "Manuel"
  analysis_metrics$annotation_mode <- annotation_mode
  
  cat("\n🔍 Étape 13/12: Génération du graphique nœud → métacluster ...\n")
  
  filtered_node_mapping <- FlowSOM:::MapDataToCodes(
    flowSOM_res$map$codes, expr_sub[, markers_of_interest]
  )
  node_ids_filtered <- filtered_node_mapping[, 1]
  
  node_to_cluster <- tapply(cluster_ids, node_ids_filtered, function(x) {
    as.numeric(names(which.max(table(x))))
  })
  node_to_cluster_df <- data.frame(
    node          = as.numeric(names(node_to_cluster)),
    cluster       = as.numeric(node_to_cluster),
    stringsAsFactors = FALSE
  )
  node_to_cluster_df$cluster_label <- cluster_labels[node_to_cluster_df$cluster]
  node_to_cluster_df <- node_to_cluster_df[order(node_to_cluster_df$node), ]
  
  n_digits <- nchar(as.character(max(node_to_cluster_df$node)))
  node_to_cluster_df$node_label <- paste0(
    "Node_",
    formatC(node_to_cluster_df$node, width = n_digits, flag = "0")
  )
  
  # Ordre décroissant des nœuds (1 en haut)
  node_levels <- paste0(
    "Node_",
    formatC(sort(unique(node_to_cluster_df$node), decreasing = TRUE),
            width = n_digits, flag = "0")
  )
  node_to_cluster_df$node_label <- factor(
    node_to_cluster_df$node_label,
    levels = node_levels
  )
  
  # Ordre des clusters par numéro
  cluster_levels <- cluster_labels[order(as.numeric(
    gsub("^(\\d+)_.*", "\\1", cluster_labels)
  ))]
  node_to_cluster_df$cluster_label <- factor(
    node_to_cluster_df$cluster_label,
    levels = cluster_levels
  )
  
  n_cls <- length(cluster_levels)
  custom_colors_alluvial <- c(
    "#E41A1C","#377EB8","#4DAF4A","#FF7F00","#984EA3",
    "#A65628","#F781BF","#00CED1","#FFD700","#1B9E77",
    "#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
    "#006400","#8B0000","#00008B","#FF4500","#2E8B57"
  )
  if (n_cls > length(custom_colors_alluvial))
    custom_colors_alluvial <- colorRampPalette(custom_colors_alluvial)(n_cls)
  cls_color_map <- setNames(
    custom_colors_alluvial[seq_along(cluster_levels)],
    cluster_levels
  )
  
  # Nœuds multiples de 10 pour les lignes majeures et étiquettes Y
  all_nodes    <- sort(unique(node_to_cluster_df$node))
  n_nodes      <- length(all_nodes)
  major_nodes  <- all_nodes[all_nodes %% 10 == 0]
  major_labels <- paste0("Node_", formatC(major_nodes, width = n_digits, flag = "0"))
  major_positions <- which(node_levels %in% major_labels)
  
  # Nombre de nœuds par cluster pour annotation (optionnel, très réduit)
  nodes_per_cluster <- table(node_to_cluster_df$cluster_label)[cluster_levels]
  
  p_node_cluster <- ggplot(
    node_to_cluster_df,
    aes(x = cluster_label, y = node_label, color = cluster_label)
  ) +
    # Lignes majeures tous les 10 nœuds (seulement)
    geom_hline(
      yintercept = major_positions,
      color = "grey70", linewidth = 0.25
    ) +
    # Points plus petits
    geom_point(size = 1.0, shape = 16, alpha = 0.6) +
    # Annotation du nombre de nœuds par cluster (en petit, en haut)
    annotate(
      "text",
      x     = seq_along(cluster_levels),
      y     = n_nodes + 0.5,
      label = paste0("n=", as.numeric(nodes_per_cluster)),
      size  = 2.0, color = "grey30", hjust = 0.5
    ) +
    scale_color_manual(values = cls_color_map, guide = "none") +
    scale_x_discrete(position = "bottom", expand = c(0, 0.5)) +
    scale_y_discrete(
      breaks = major_labels,
      labels = formatC(major_nodes, width = n_digits, flag = "0"),
      expand = c(0, 0)    # pas d'espace supplémentaire vertical
    ) +
    theme_minimal(base_size = 8) +
    theme(
      axis.text.x = element_text(
        angle = 45, hjust = 1, vjust = 1, size = 5, face = "bold"
      ),
      axis.text.y = element_text(size = 5, family = "mono"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 7, face = "bold"),
      axis.ticks.y = element_line(color = "grey40", linewidth = 0.2),
      axis.ticks.length.y = unit(0.1, "cm"),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste("Affectation des nœuds FlowSOM aux métaclusters -", sample_name),
      y = "Numéro de nœud"
    )
  
  # Dimensions A3 paysage (16.54" × 11.69") – convient pour 400 nœuds et 30 clusters
  alluvial_file <- file.path(output_dir, "Alluvial_Node_to_Cluster.pdf")
  tryCatch({
    CairoPDF(alluvial_file, width = 16.54, height = 11.69)
    print(p_node_cluster)
    dev.off()
    cat("✅ Graphique nœud→cluster format A3 sauvegardé :", alluvial_file, "\n")
  }, error = function(e) {
    if (length(dev.list()) > 0) dev.off()
    cat("❌ Erreur graphique :", e$message, "\n")
  })
  if (file.exists(alluvial_file) && file.info(alluvial_file)$size > 0) open_pdf(alluvial_file)
  
  if ("FolderID" %in% colnames(expr_sub) && "FileID" %in% colnames(expr_sub)) {
    perform_condition_analysis(
      expr_sub = expr_sub,
      cluster_ids = cluster_ids,
      cluster_labels = cluster_labels,
      sample_name = sample_name,
      output_dir = output_dir,
      desc = desc,
      umap_df = umap_df
    )
  }
  
  cluster_summary <- data.frame(
    Cluster = 1:length(cluster_annotations),
    Cluster_Label = cluster_labels,
    Count = as.numeric(table(cluster_ids)),
    Percentage = round(100 * as.numeric(table(cluster_ids)) / sum(table(cluster_ids)), 2),
    Annotation = cluster_annotations,
    Confidence_Score = if (!is.null(annotation_results) && length(annotation_results$scores) == length(cluster_annotations)) {
      annotation_results$scores
    } else {
      rep(NA, length(cluster_annotations))
    }
  )
  
  annotation_mode_code <- if (annotation_mode == "Manuel") 2 else 1
  annotation_report_file <- generate_annotation_report(
    cluster_summary, annotation_mode_code, output_dir, sample_name
  )
  cluster_summary_file <- file.path(output_dir, "Cluster_Summary.csv")
  write.csv(cluster_summary, cluster_summary_file, row.names = FALSE)
  
  if (length(channels_to_exclude) > 0) {
    excluded_channels_info <- data.frame(
      Channel = channels_to_exclude,
      Marker_Name = sapply(channels_to_exclude, function(x) get_marker_name(x, marker_descriptions)),
      Reason = "Excluded by user choice"
    )
    write.csv(excluded_channels_info, file.path(output_dir, "Excluded_Channels.csv"), row.names = FALSE)
  }
  
  analysis_metrics$end_time <- Sys.time()
  analysis_metrics$total_time <- difftime(analysis_metrics$end_time, analysis_metrics$start_time, units = "secs")
  format_time <- function(seconds) {
    hours <- floor(seconds / 3600)
    minutes <- floor((seconds - hours * 3600) / 60)
    seconds <- round(seconds - hours * 3600 - minutes * 60, 2)
    sprintf("%02d:%02d:%05.2f", hours, minutes, seconds)
  }
  analysis_metrics$total_time_formatted <- format_time(as.numeric(analysis_metrics$total_time))
  analysis_metrics$script_name <- script_name
  analysis_metrics$analysis_date <- format(analysis_metrics$start_time, "%Y-%m-%d", tz = "Europe/Brussels")
  analysis_metrics$analysis_time <- format(analysis_metrics$start_time, "%H:%M:%S", tz = "Europe/Brussels")
  
  metrics_list <- list()
  metrics_list[["script_name"]] <- analysis_metrics$script_name
  metrics_list[["analysis_date"]] <- analysis_metrics$analysis_date
  metrics_list[["analysis_time"]] <- analysis_metrics$analysis_time
  metrics_list[["total_time"]] <- analysis_metrics$total_time_formatted
  metrics_list[["seed"]] <- analysis_metrics$seed
  metrics_list[["initial_cell_count"]] <- analysis_metrics$initial_cell_count
  metrics_list[["after_gate_CD8_CD56"]] <- analysis_metrics$after_gate_CD8_CD56
  metrics_list[["after_gate_CD3_UnconvT"]] <- analysis_metrics$after_gate_CD3_UnconvT
  metrics_list[["after_doublet_removal"]] <- analysis_metrics$after_doublet_removal
  metrics_list[["after_fsc_ssc_gating"]] <- analysis_metrics$after_fsc_ssc_gating
  metrics_list[["subsampled_cell_count"]] <- analysis_metrics$subsampled_cell_count
  metrics_list[["subsampling_method"]] <- analysis_metrics$subsampling_method
  metrics_list[["batch_normalized"]] <- ifelse(is.null(analysis_metrics$batch_normalized), FALSE, analysis_metrics$batch_normalized)
  metrics_list[["excluded_channels_count"]] <- length(channels_to_exclude)
  if (length(channels_to_exclude) > 0) metrics_list[["excluded_channels"]] <- paste(channels_to_exclude, collapse = "; ")
  metrics_list[["flowsom_xdim"]] <- analysis_metrics$flowsom_xdim
  metrics_list[["flowsom_ydim"]] <- analysis_metrics$flowsom_ydim
  metrics_list[["flowsom_rlen"]] <- analysis_metrics$flowsom_rlen
  metrics_list[["umap_n_neighbors"]] <- analysis_metrics$umap_n_neighbors
  metrics_list[["umap_min_dist"]] <- analysis_metrics$umap_min_dist
  metrics_list[["umap_metric"]] <- analysis_metrics$umap_metric
  metrics_list[["umap_n_epochs"]] <- analysis_metrics$umap_n_epochs
  metrics_list[["umap_n_threads"]] <- analysis_metrics$umap_n_threads
  metrics_list[["optimal_k"]] <- analysis_metrics$optimal_k
  metrics_list[["annotation_mode"]] <- analysis_metrics$annotation_mode
  
  metrics_df <- data.frame(Metric = names(metrics_list), Value = unlist(metrics_list))
  write.csv(metrics_df, file.path(output_dir, "analysis_metrics.csv"), row.names = FALSE)
  
  cat("✅ Analyse terminée avec succès!\n")
  cat("✅ Tous les résultats ont été sauvegardés dans :", output_dir, "\n")
  cat("=== FIN DE L'ANALYSE ===\n")
  
  if (requireNamespace("renv", quietly = TRUE)) {
    renv_snapshot_path <- file.path(output_dir,
                                    paste0("renv_", format(Sys.time(), "%Y-%m-%d_%Hh%M"), ".lock"))
    renv::snapshot(lockfile = renv_snapshot_path, type = "all")
    cat("✅ Copie du renv.lock horodatée sauvegardée :", renv_snapshot_path, "\n")
  } else {
    cat("⚠️ Package renv non disponible. Pas de snapshot horodaté.\n")
  }
  
  session_info_text <- capture.output({
    cat("Date :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("R version :", paste(R.version$major, R.version$minor, sep = "."), "\n")
    cat("OS :", sessionInfo()$running, "\n\n")
    sessionInfo()
  })
  writeLines(session_info_text, file.path(output_dir, "session_info.txt"))
  cat("✅ Session info sauvegardée :", file.path(output_dir, "session_info.txt"), "\n")
  
  return(list(
    live_dead_gate = NULL,
    doublet_gate_coords = doublet_gate_coords,
    fsc_ssa_gate_coords = fsc_ssa_gate_coords,
    transformation_params = transformation_params,
    optimal_k = optimal_k,
    expr_sub = expr_sub,
    cluster_ids = cluster_ids,
    cluster_labels = cluster_labels,
    markers_of_interest = markers_of_interest,
    sample_name = sample_name,
    output_dir = output_dir,
    umap_res = umap_res,
    marker_descriptions = marker_descriptions,
    fs_final = fs_trans[final_cells_mask_after_ssa, ],
    n_cells_final = nrow(expr_sub),
    cluster_annotations = cluster_annotations,
    marker_thresholds = NULL,
    cluster_medians_data = cluster_medians_data,
    analysis_metrics = analysis_metrics,
    annotation_details = annotation_results$details,
    excluded_channels = channels_to_exclude,
    fsc_ssc_filtered_file = fsc_ssc_filtered_file,
    annotation_results = annotation_results,
    umap_df = umap_df,
    desc = desc
  ))
}

# ============================================================================
# SECTION 12: FONCTIONS DE RAPPORT
# ============================================================================
generate_detailed_annotation_report <- function(annotation_results, cluster_medians,
                                                marker_descriptions, output_dir, sample_name) {
  report_file <- file.path(output_dir, "Detailed_Annotation_Report.txt")
  sink(report_file)
  cat("=== RAPPORT DÉTAILLÉ D'ANNOTATION DES CLUSTERS ===\n\n")
  cat("Échantillon :", sample_name, "\n")
  cat("Date :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("RÉSUMÉ DES ANNOTATIONS :\n")
  cat("=======================\n")
  annotation_summary <- table(annotation_results$annotations)
  for (annotation in names(annotation_summary)) {
    cat(paste(annotation, ":", annotation_summary[annotation], "clusters\n"))
  }
  cat("\nDÉTAIL PAR CLUSTER :\n")
  cat("====================\n")
  for (i in 1:length(annotation_results$annotations)) {
    cat(paste("\n--- Cluster", i, "---\n"))
    cat(paste("Annotation :", annotation_results$annotations[i], "\n"))
    cat("Expression des marqueurs :\n")
    for (marker in colnames(cluster_medians)) {
      if (marker %in% colnames(cluster_medians)) {
        expr_value <- cluster_medians[i, marker]
        marker_name <- get_marker_name(marker, marker_descriptions)
        cat(paste(" ", marker_name, ":", round(expr_value, 3), "\n"))
      }
    }
  }
  sink()
  return(report_file)
}

generate_annotation_report <- function(cluster_summary, annotation_mode, output_dir, sample_name) {
  report_file <- file.path(output_dir, "Annotation_Report.txt")
  sink(report_file)
  cat("=== RAPPORT D'ANNOTATION DES CLUSTERS ===\n\n")
  cat("Échantillon :", sample_name, "\n")
  cat("Date :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Mode d'annotation :", c("Automatique", "Manuel")[annotation_mode], "\n\n")
  cat("RÉSUMÉ DES CLUSTERS :\n")
  cat("=====================\n")
  print(cluster_summary)
  cat("\nDÉTAIL PAR CLUSTER :\n")
  cat("====================\n")
  for (i in 1:nrow(cluster_summary)) {
    cat(paste("\nCluster", cluster_summary$Cluster[i], ":\n"))
    cat(paste(" Label :", cluster_summary$Cluster_Label[i], "\n"))
    cat(paste(" Cellules :", cluster_summary$Count[i], "\n"))
    cat(paste(" Pourcentage :", cluster_summary$Percentage[i], "%\n"))
    cat(paste(" Annotation :", cluster_summary$Annotation[i], "\n"))
  }
  sink()
  return(report_file)
}

# ============================================================================
# SECTION 13: PROGRAMME PRINCIPAL
# ============================================================================
main <- function() {
  temp_log <- tempfile(pattern = "console_log_", fileext = ".txt")
  sink(temp_log, split = TRUE)
  
  cat("=== Démarrage du script FCS_Automated_Clustering (version concaténée) ===\n\n")
  cat("Entrez le nom du script utilisé: ")
  script_name <- readline()
  
  if (Sys.getenv("DEEPSEEK_API_KEY") == "") {
    cat("Aucune clé API DeepSeek trouvée dans l'environnement.\n")
    cat("Entrez votre clé API (ou laissez vide pour utiliser le mode manuel plus tard) : ")
    api_key <- readline()
    if (api_key != "") Sys.setenv(DEEPSEEK_API_KEY = api_key)
  }
  
  cat("Entrez le chemin complet du fichier FCS à analyser: ")
  fcs_file <- readline()
  if (!file.exists(fcs_file)) stop("❌ Le fichier spécifié n'existe pas.")
  
  output_dir_base <- file.path(dirname(fcs_file), "outputs")
  cat("✅ Dossier de sortie:", output_dir_base, "\n")
  
  cat("\n", paste(rep("=", 50), collapse = ""), "\n", sep = "")
  cat("DÉBUT DE L'ANALYSE PRINCIPALE\n")
  cat(paste(rep("=", 50), collapse = ""), "\n", sep = "")
  
  analysis_results <- analyze_fcs_file(
    fcs_file, output_dir_base, script_name = script_name,
    log_temp_file = temp_log
  )
  
  if (is.null(analysis_results)) {
    if (sink.number() > 0) sink()
    stop("❌ L'analyse a échoué.")
  }
  
  while (sink.number() > 0) sink()
  
  cat("\n", paste(rep("=", 50), collapse = ""), "\n", sep = "")
  cat("ANALYSE TERMINÉE AVEC SUCCÈS!\n")
  cat(paste(rep("=", 50), collapse = ""), "\n", sep = "")
  
  return(analysis_results)
}

results <- main()
