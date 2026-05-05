# global.R — loaded once, shared across sessions

library(shiny)
library(ggiraph)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# ── Data paths ────────────────────────────────────────────────────────────────


# ── Load all supplemental tables and data ─────────────────────────────────────
# Local deployment
# DATA_DIR <- "/Users/johanna/GitHub/mass_spec/analysis/do_cr_paper/shiny-app"
# tbl_s1_annotations  <- read_csv(file.path(DATA_DIR, "data/Table_S1_CompoundAnnotations.csv"),  show_col_types = FALSE)
# tbl_s2_qtls         <- readRDS(file.path(DATA_DIR, "data/Table_S9_QTLs.Rds"))
# tbl_s4_diet_age     <- readRDS(file.path(DATA_DIR, "data/Table_S2_Diet_Age_Model.Rds"))
# tbl_s6_pll_aging    <- readRDS(file.path(DATA_DIR, "data/Table_S4_PLL_Aging_Dying_Model.Rds"))
# all_data <- readRDS(file.path(DATA_DIR, "data/shiny-app-compressed-data.Rds"))

# App deployment
tbl_s1_annotations <- read_csv("data/Table_S1_CompoundAnnotations.csv", show_col_types = FALSE)
tbl_s2_qtls <- readRDS("data/Table_S9_QTLs.Rds")
tbl_s4_diet_age <- readRDS("data/Table_S2_Diet_Age_Model.Rds")
tbl_s6_pll_aging <- readRDS("data/Table_S4_PLL_Aging_Dying_Model.Rds")
all_data <- readRDS("data/shiny-app-compressed-data.Rds")

# Refactor diet names
all_data <- all_data %>%
  dplyr::mutate(diet = dplyr::recode(diet,
    "1D" = "IF-1D", "2D" = "IF-2D", "20" = "CR-20", "40" = "CR-40"
  ))

# ── Compound name list for selector ──────────────────────────────────────────
compound_choices <- sort(unique(tbl_s1_annotations$name_use))

# ── S4 term choices derived from data (after Fasting remapping) ───────────────
.tmp_terms <- tbl_s4_diet_age$Term
.tmp_terms[grepl("FastingFast20", .tmp_terms)] <- "Fasting-Diet20"
.tmp_terms[grepl("FastingFast40", .tmp_terms)] <- "Fasting-Diet40"
.tmp_terms[.tmp_terms == "(Intercept)"] <- "Intercept"
.tmp_terms[.tmp_terms == "Diet1D"] <- "IF-1D"
.tmp_terms[.tmp_terms == "Diet2D"] <- "IF-2D"
.tmp_terms[.tmp_terms == "Diet20"] <- "CR-20"
.tmp_terms[.tmp_terms == "Diet40"] <- "CR-40"
s4_term_choices <- sort(unique(.tmp_terms))
.canonical <- c("AgeYear1", "AgeYear2", "AgeYear3", "IF-1D", "IF-2D", "CR-20", "CR-40")
s4_term_choices <- c(
  intersect(.canonical, s4_term_choices),
  sort(setdiff(s4_term_choices, .canonical))
)
rm(.tmp_terms, .canonical)
s4_n_platforms <- length(unique(tbl_s4_diet_age$Platform))
s4_compound_choices <- sort(unique(tbl_s4_diet_age$Compound))

# ── S6 term choices and filters ───────────────────────────────────────────────
.recode_diet_terms <- function(x) {
  x[x == "Diet1D"] <- "IF-1D"
  x[x == "Diet2D"] <- "IF-2D"
  x[x == "Diet20"] <- "CR-20"
  x[x == "Diet40"] <- "CR-40"
  x
}

.s6_terms_aging <- tbl_s6_pll_aging$Term[tbl_s6_pll_aging$Phase == "Aging"]
.s6_terms_aging[.s6_terms_aging == "(Intercept)"] <- "Intercept"
.s6_terms_aging <- .recode_diet_terms(.s6_terms_aging)
s6_aging_term_choices <- sort(unique(.s6_terms_aging))
s6_aging_term_choices <- c(
  intersect("PLL_Aging", s6_aging_term_choices),
  sort(setdiff(s6_aging_term_choices, "PLL_Aging"))
)

.s6_terms_dying <- tbl_s6_pll_aging$Term[tbl_s6_pll_aging$Phase == "Dying"]
.s6_terms_dying[.s6_terms_dying == "(Intercept)"] <- "Intercept"
.s6_terms_dying <- .recode_diet_terms(.s6_terms_dying)
s6_dying_term_choices <- sort(unique(.s6_terms_dying))
s6_dying_term_choices <- c(
  intersect("PLL_Dying", s6_dying_term_choices),
  sort(setdiff(s6_dying_term_choices, "PLL_Dying"))
)

s6_normalization_choices <- sort(unique(tbl_s6_pll_aging$Normalization))
s6_n_platforms <- length(unique(tbl_s6_pll_aging$Platform))
rm(.s6_terms_aging, .s6_terms_dying)
# ── PLL trajectory compound list ─────────────────────────────────────────────
pll_compound_choices <- sort(unique(
  tbl_s6_pll_aging$Compound[tbl_s6_pll_aging$Term == "PLL_Aging"]
))

qtl_compound_choices <- sort(unique(tbl_s2_qtls$name_use))
qtl_n_platforms <- length(unique(tbl_s2_qtls$modality))
qtl_n_ages <- length(unique(substr(
  tbl_s2_qtls$trait,
  nchar(tbl_s2_qtls$trait) - 4,
  nchar(tbl_s2_qtls$trait)
)))

# ── Download table metadata ───────────────────────────────────────────────────
download_meta <- list(
  list(
    id    = "s4",
    label = "Regression Results",
    file  = "Table_S2_Diet_Age_Model.Rds",
    desc  = "Linear model results for diet and age effects on plasma metabolite, lipid, and protein features"
  ),
  list(
    id    = "s6",
    label = "PLL Trajectories",
    file  = "Table_S4_PLL_Aging_Dying_Model.Rds",
    desc  = "Linear model results capturing Aging and Dying slopes from PLL-standardized timepoints"
  ),
  list(
    id    = "s2",
    label = "QTL Results",
    file  = "Table_S9_QTLs.Rds",
    desc  = "Quantitative trait loci (QTL) mapping results for all metabolite, lipid, and protein features"
  )
)
