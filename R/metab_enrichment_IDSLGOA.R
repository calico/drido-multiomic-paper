# This function is derived from code and keys published for https://goa.idsl.site/goa/
# Priyanka Mahajan, Oliver Fiehn, and Dinesh Barupal. "IDSL. GOA: gene ontology analysis for interpreting metabolomic datasets." Scientific Reports 14, no. 1 (2024): 1299.(Link).

docr_inchikey_enrichment <- function(input_inchi_vector,
                                     go_names_filepath = "go_names.RData",
                                     go_inchi_key_filepath = "all_go_child_cpd_list_met_ik14.RData",
                                     input_universe = NULL) {
  # load R objects
  load(go_names_filepath)
  load(go_inchi_key_filepath)

  all_cpds <- unique(unlist(all_go_child_cpd_list_met_ik14))

  # Optional - user provided universe
  if (!is.null(input_universe)) {
    input_universe <- unique(as.character(sapply(
      input_universe,
      function(ik) {
        strsplit(ik, "-")[[1]][1]
      }
    )))

    # all compounds (universe)
    all_cpds <- intersect(all_cpds, input_universe)

    # universe-filtered compounds set lists
    all_go_child_cpd_list_met_ik14 <- lapply(
      all_go_child_cpd_list_met_ik14,
      function(cpd_list) {
        return(intersect(cpd_list, input_universe))
      }
    )

    # remove set lists that are reduced to length 0
    all_go_child_cpd_list_met_ik14 <- all_go_child_cpd_list_met_ik14[sapply(all_go_child_cpd_list_met_ik14, length) > 0]
  }

  ### INPUT LIST -----
  ik_vec <- unique(as.character(sapply(input_inchi_vector, function(ik) {
    suppressWarnings(strsplit(ik, "-")[[1]][1])
  })))

  xn <- length(all_cpds) # universe list size
  xk <- length(ik_vec[ik_vec %in% all_cpds]) # input list size

  # run GO if we have at least four compounds input list
  if (xk > 3) {
    go_analysis_list <- do.call(
      rbind,
      lapply(
        1:length(all_go_child_cpd_list_met_ik14),
        function(x) {
          # overlap size
          xm <- sum(all_go_child_cpd_list_met_ik14[[x]] %in% ik_vec)
          # values, to identify pathways that capture redundant information
          xm_values <- intersect(all_go_child_cpd_list_met_ik14[[x]], ik_vec) %>%
            sort() %>%
            paste0(., collapse = "|")
          # pathway list size
          pl <- length(all_go_child_cpd_list_met_ik14[[x]])
          # fold enrichment
          fe <- (xm / xk) / (pl / xn)

          res <- c(
            xm, pl, xk, xn,
            phyper(xm - 1, pl, xn - pl, xk, lower.tail = F),
            fe,
            xm_values
          )
          return(res)
        }
      )
    )

    go_res_df <- data.frame(go_analysis_list,
      stringsAsFactors = F
    ) %>%
      {
        setNames(., c(
          "Overlap", "Process", "Input", "All",
          "Pvalue", "Fold Enrichment", "Overlap_Compounds"
        ))
      } %>%
      dplyr::mutate(across(
        c(
          "Overlap", "Process", "Input", "All",
          "Pvalue", "Fold Enrichment"
        ),
        as.numeric
      )) %>%
      dplyr::mutate(
        GO_term = names(all_go_child_cpd_list_met_ik14),
        SetSizeRatio = Process / xn
      ) %>%
      dplyr::filter(Overlap > 3)
    go_res_df$GO_name <- as.character(go_names[go_res_df$GO_term])

    return(go_res_df)
  }
  return(NULL)
}
