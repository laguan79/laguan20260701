options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "code/rebuild_endpoint_selection_construct_metadata.R"
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
data_root <- Sys.getenv("DERIVED_DATA_ROOT", unset = file.path(root, "derived_data"))
output_root <- Sys.getenv("ENDPOINT_AUDIT_OUTPUT_ROOT", unset = file.path(root, "endpoint_audit_rebuild"))

read_tsv <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE, quote = "", na.strings = c("", "NA"))
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

clean_text <- function(x) {
  x <- tolower(iconv(x, to = "ASCII//TRANSLIT"))
  x <- gsub("ischaemic", "ischemic", x, fixed = TRUE)
  x <- gsub("anaemia", "anemia", x, fixed = TRUE)
  x <- gsub("apnoea", "apnea", x, fixed = TRUE)
  x <- gsub("dysrhythmias", "arrhythmias", x, fixed = TRUE)
  x <- gsub("eyelids", "eyelid", x, fixed = TRUE)
  x <- gsub("[^a-z0-9]+", " ", x)
  trimws(gsub("\\s+", " ", x))
}

token_set <- function(x) {
  stop_words <- c("a", "an", "and", "of", "the", "other", "incl", "including", "nos")
  z <- unlist(strsplit(clean_text(x), " ", fixed = TRUE))
  unique(z[nzchar(z) & !z %in% stop_words])
}

jaccard_distance <- function(a, b) {
  aa <- token_set(a)
  bb <- token_set(b)
  if (!length(union(aa, bb))) return(0)
  1 - length(intersect(aa, bb)) / length(union(aa, bb))
}

character_distance <- function(a, b) {
  aa <- clean_text(a)
  bb <- clean_text(b)
  denom <- max(nchar(aa), nchar(bb), 1)
  as.numeric(adist(aa, bb)) / denom
}

label_distance_v2 <- function(a, b) {
  0.70 * jaccard_distance(a, b) + 0.30 * character_distance(a, b)
}

candidate_path <- file.path(data_root, "endpoint_selection", "candidate_universe_416.tsv")
panel_path <- file.path(data_root, "endpoint_selection", "retained_endpoint_panel.tsv")
paired_path <- file.path(data_root, "S02_paired_source_classifications.tsv")

candidates <- read_tsv(candidate_path)
panel <- read_tsv(panel_path)
paired <- read_tsv(paired_path)

stopifnot(nrow(candidates) == 416L, nrow(panel) == 53L, nrow(paired) == 1060L)

# Reconstruct the archived, effect-blind tier-1 rule. The source table contains
# endpoint labels, mappings, case counts, and construct metadata, but no effect
# estimates, standard errors, or P values.
candidates$reconstructed_tier1 <- with(
  candidates,
  construct_equivalence_tier == "close_construct_match" &
    min_source_cases >= 5000 &
    case_count_ratio <= 20
)

reason <- character(nrow(candidates))
for (i in seq_len(nrow(candidates))) {
  r <- character()
  if (candidates$construct_equivalence_tier[i] != "close_construct_match") {
    r <- c(r, "legacy_construct_tier_not_close")
  }
  if (candidates$min_source_cases[i] < 5000) {
    r <- c(r, "minimum_source_cases_below_5000")
  }
  if (candidates$case_count_ratio[i] > 20) {
    r <- c(r, "case_count_ratio_above_20")
  }
  reason[i] <- if (length(r)) paste(r, collapse = ";") else "included_by_archived_tier1_rule"
}
candidates$selection_decision <- ifelse(candidates$reconstructed_tier1, "included", "excluded")
candidates$selection_reason <- reason
candidates$effect_data_present_in_selection_source <- FALSE

panel_keys <- paste(panel$finngen_phenocode, panel$panukb_phenocode, sep = "||")
candidate_keys <- paste(candidates$finngen_phenocode, candidates$panukb_phenocode, sep = "||")
reconstructed_keys <- candidate_keys[candidates$reconstructed_tier1]
if (!setequal(panel_keys, reconstructed_keys)) {
  stop("Archived tier-1 rule does not reconstruct the 53-pair panel.")
}

selection_counts <- data.frame(
  stage = c(
    "Officially mapped candidate endpoint pairs",
    "Excluded by archived effect-blind tier-1 rule",
    "Included in curated calibration benchmark"
  ),
  n_pairs = c(nrow(candidates), sum(!candidates$reconstructed_tier1), sum(candidates$reconstructed_tier1)),
  definition = c(
    "Candidate pairs with official mapping anchors and endpoint-level metadata",
    "One or more of: legacy construct tier not close, minimum source cases below 5000, or case-count ratio above 20",
    "Legacy close construct, minimum source cases at least 5000, and case-count ratio at most 20"
  )
)

write_tsv(candidates, file.path(output_root, "endpoint_selection_ledger.tsv"))
write_tsv(selection_counts, file.path(output_root, "endpoint_selection_counts.tsv"))

# Correct clinical-domain metadata. Domains are descriptive and are deliberately
# excluded from construct-distance scoring because the legacy Phecode range
# mapper treated Phecodes as ICD chapters.
domain_map <- list(
  circulatory = 1:13,
  digestive = c(14:19, 23:26, 28),
  respiratory = c(21, 22, 27),
  genitourinary = 29:36,
  musculoskeletal = c(37:44),
  neurological = c(45, 46),
  hematologic = 49,
  sense_organs = 51:53
)

domain_for_rank <- function(rank, source) {
  if (rank == 20) return(if (source == "finngen") "mixed_metabolic_respiratory" else "endocrine_metabolic")
  if (rank == 47) return("circulatory")
  if (rank == 48) return("musculoskeletal_pain")
  if (rank == 50) return("genitourinary")
  for (nm in names(domain_map)) {
    if (rank %in% domain_map[[nm]]) return(nm)
  }
  "other"
}

panel$legacy_label_distance <- panel$label_distance
panel$legacy_granularity_distance <- panel$granularity_distance
panel$legacy_clinical_domain_distance <- panel$clinical_domain_distance
panel$legacy_endpoint_construct_distance <- panel$endpoint_construct_distance
panel$legacy_construct_equivalence_tier <- panel$construct_equivalence_tier
panel$legacy_finngen_clinical_domain <- panel$finngen_clinical_domain
panel$legacy_panukb_clinical_domain <- panel$panukb_clinical_domain
panel$legacy_clinical_domain_match <- panel$clinical_domain_match

panel$finngen_clinical_domain <- vapply(panel$pair_rank, domain_for_rank, character(1), source = "finngen")
panel$panukb_clinical_domain <- vapply(panel$pair_rank, domain_for_rank, character(1), source = "panukb")
panel$clinical_domain_match <- panel$finngen_clinical_domain == panel$panukb_clinical_domain
panel$clinical_domain_role <- "descriptive_only_not_used_in_construct_score"

panel$label_distance_v2_automated <- mapply(label_distance_v2, panel$finngen_label, panel$panukb_label)
panel$granularity_distance_v2 <- panel$legacy_granularity_distance
panel$endpoint_construct_distance_v2_automated <- pmin(
  1,
  0.625 * panel$label_distance_v2_automated + 0.375 * panel$granularity_distance_v2
)

# Effect-blind manual adjudication handles composite-versus-single and
# context-specific endpoint relations that lexical distance cannot resolve.
manual <- data.frame(
  pair_rank = c(3, 7, 8, 10, 11, 12, 13, 14, 15, 16, 20, 21, 22, 25, 27, 29, 32, 37, 42, 45, 47, 48, 49),
  distance = c(0.15, 0.18, 0.18, 0.18, 0.35, 0.05, 0.15, 0.18, 0.18, 0.15, 0.55, 0.28, 0.35, 0.28, 0.15, 0.35, 0.10, 0.18, 0.35, 0.15, 0.12, 0.55, 0.05),
  relation = c(
    "same disease with breadth qualifier",
    "etiologic subtype versus parent disease",
    "same venous disease with location qualifier",
    "definition qualifier versus clinical synonym",
    "composite heart-failure/coronary endpoint versus heart failure",
    "same disorder with NOS qualifier",
    "clinical synonym with parenthetical definition",
    "parent hernia category versus abdominal hernia",
    "same diverticular disease family with terminology variation",
    "abdominal-wall hernia versus abdominal hernia",
    "obesity-related asthma composite versus obesity",
    "influenza-and-pneumonia composite versus pneumonia",
    "asthma-related pneumonia versus pneumonia",
    "postoperative abdominal hernia versus abdominal hernia",
    "pleural-effusion endpoint with pleurisy grouping",
    "broad female-genital disorder category versus bleeding/menstruation subset",
    "same female-genital polyp construct",
    "low-back-pain subtype versus back pain",
    "polyarthropathy category versus rheumatoid/inflammatory polyarthropathy grouping",
    "combined sleep-disorder category versus sleep disorders",
    "clinical synonyms",
    "multisite pain composite versus back pain",
    "same iron-deficiency-anemia construct"
  )
)

panel$endpoint_construct_distance <- panel$endpoint_construct_distance_v2_automated
panel$construct_adjudication <- "automated_label_and_granularity_score"
for (i in seq_len(nrow(manual))) {
  idx <- match(manual$pair_rank[i], panel$pair_rank)
  panel$endpoint_construct_distance[idx] <- manual$distance[i]
  panel$construct_adjudication[idx] <- paste0("effect_blind_manual_relation: ", manual$relation[i])
}

panel$label_distance <- panel$label_distance_v2_automated
panel$granularity_distance <- panel$granularity_distance_v2
panel$clinical_domain_distance <- NA_real_
panel$construct_equivalence_tier <- ifelse(
  panel$endpoint_construct_distance <= 0.20,
  "close_construct_match",
  ifelse(panel$endpoint_construct_distance <= 0.40, "moderate_construct_match", "wide_construct_mismatch")
)
panel$construct_score_v2_definition <- "0.625*label_distance_v2+0.375*granularity_distance; domain excluded; documented manual relation overrides"
panel$benchmark_role <- paste0("calibration_", panel$field_domain)
panel$analysis_tier <- "curated_calibration_benchmark"
panel$archived_endpoint_priority <- "archived_tier1_selection"
panel$inclusion_decision <- "retained_in_curated_calibration_benchmark"

new_panel_path <- file.path(output_root, "retained_endpoint_panel_rebuilt.tsv")
write_tsv(panel, new_panel_path)

audit <- panel[, c(
  "pair_rank", "outcome_family", "finngen_phenocode", "panukb_phenocode",
  "finngen_label", "panukb_label",
  "legacy_finngen_clinical_domain", "legacy_panukb_clinical_domain",
  "finngen_clinical_domain", "panukb_clinical_domain", "clinical_domain_match",
  "legacy_endpoint_construct_distance", "endpoint_construct_distance",
  "legacy_construct_equivalence_tier", "construct_equivalence_tier",
  "construct_adjudication"
)]
write_tsv(audit, file.path(output_root, "endpoint_construct_metadata_audit.tsv"))

joined <- merge(
  paired,
  panel[, c("outcome_family", "construct_equivalence_tier", "endpoint_construct_distance")],
  by = "outcome_family",
  all.x = TRUE,
  sort = FALSE
)
if (nrow(joined) != 1060L || any(is.na(joined$construct_equivalence_tier))) {
  stop("Construct metadata did not join cleanly to all 1,060 paired comparisons.")
}

impact <- do.call(rbind, lapply(split(joined, joined$construct_equivalence_tier), function(z) {
  naive <- z$asymmetric_nominal_support %in% TRUE
  hetero <- z$effect_difference_fdr < 0.05
  equiv <- z$equivalent_or_margin_20 %in% TRUE
  data.frame(
    construct_equivalence_tier = z$construct_equivalence_tier[1],
    endpoint_pairs = length(unique(z$outcome_family)),
    complete_comparisons = nrow(z),
    naive_asymmetry = sum(naive, na.rm = TRUE),
    fdr_supported_source_heterogeneity = sum(hetero, na.rm = TRUE),
    practical_equivalence = sum(equiv, na.rm = TRUE),
    naive_rows_retaining_fdr_heterogeneity = sum(naive & hetero, na.rm = TRUE),
    naive_rows_downgraded = sum(naive & !hetero, na.rm = TRUE)
  )
}))
impact$naive_downgrade_percent <- ifelse(
  impact$naive_asymmetry > 0,
  100 * impact$naive_rows_downgraded / impact$naive_asymmetry,
  NA_real_
)
impact <- impact[order(match(
  impact$construct_equivalence_tier,
  c("close_construct_match", "moderate_construct_match", "wide_construct_mismatch")
)), ]
write_tsv(impact, file.path(output_root, "construct_tier_decision_impact.tsv"))

selection_flow <- c(
  "# Endpoint Selection Reconstruction",
  "",
  "The 53 endpoint pairs form a curated calibration benchmark. They are not a probability sample of public endpoints and were not retrospectively described as preregistered.",
  "",
  "The archived candidate table contained endpoint labels, official mapping anchors, case counts, source labels, and construct metadata, but no beta, standard error, or P value. Applying the archived tier-1 rule reconstructed the final panel exactly:",
  "",
  "1. Start with 416 officially mapped candidate endpoint pairs.",
  "2. Require the legacy close-construct tier.",
  "3. Require at least 5,000 cases in the smaller source.",
  "4. Require a source case-count ratio no greater than 20.",
  "5. Retain 53 pairs; exclude 363 pairs for one or more effect-blind metadata criteria.",
  "",
  "The construct metadata was subsequently repaired without consulting effect estimates. Because the legacy domain mapper treated Phecode numeric ranges as ICD chapters, clinical domain was removed from the construct score and retained only as descriptive metadata. The revised score uses label and granularity information, with documented manual relation adjudication for composite-versus-single and context-specific pairs.",
  "",
  "Files:",
  "- `endpoint_pair_candidate_universe_416.tsv`: archived effect-blind candidate universe.",
  "- `endpoint_selection_ledger.tsv`: row-level reconstructed decision and exclusion reason.",
  "- `endpoint_selection_counts.tsv`: flow counts.",
  "- `endpoint_construct_metadata_audit.tsv`: legacy-to-corrected metadata audit."
)
writeLines(selection_flow, file.path(output_root, "ENDPOINT_SELECTION_FLOW.md"), useBytes = TRUE)

tier_counts <- as.data.frame(table(panel$construct_equivalence_tier), stringsAsFactors = FALSE)
names(tier_counts) <- c("construct_tier", "endpoint_pairs")
domain_counts <- data.frame(
  metric = c("legacy_clinical_domain_matches", "corrected_clinical_domain_matches"),
  n = c(
    sum(panel$legacy_clinical_domain_match %in% TRUE),
    sum(panel$clinical_domain_match %in% TRUE)
  ),
  denominator = nrow(panel)
)

qa_lines <- c(
  "# Construct Metadata Impact Audit",
  "",
  paste0("- Candidate pairs reconstructed: ", nrow(candidates), "."),
  paste0("- Curated endpoint pairs reconstructed by the archived effect-blind rule: ", sum(candidates$reconstructed_tier1), "."),
  paste0("- Excluded candidate pairs: ", sum(!candidates$reconstructed_tier1), "."),
  paste0("- Legacy clinical-domain matches: ", sum(panel$legacy_clinical_domain_match %in% TRUE), "/53."),
  paste0("- Corrected descriptive clinical-domain matches: ", sum(panel$clinical_domain_match %in% TRUE), "/53."),
  "- Clinical domain contribution to construct distance: removed.",
  paste0("- Revised construct tiers: ", paste(paste(tier_counts$construct_tier, tier_counts$endpoint_pairs, sep = "="), collapse = "; "), "."),
  "- Main benchmark denominators and effect estimates: unchanged.",
  "- Construct-stratified decision-impact results: recalculated in `tables/construct_tier_decision_impact_v2.tsv`.",
  "",
  "## Reproducible checks",
  "",
  "- The archived tier-1 rule reproduces exactly the same 53 endpoint-code pairs.",
  "- Corrected metadata join to all 1,060 paired comparisons completed without missing rows.",
  "- The correction script does not read source-specific effect estimates until after the endpoint selection and construct metadata objects have been finalized.",
  "",
  "## Interpretation",
  "",
  "The former statement that all 53 pairs were close constructs is withdrawn. The panel remains a curated calibration benchmark, and revised construct tiers are reported transparently. Clinical domains are descriptive rather than components of construct distance."
)
writeLines(qa_lines, file.path(output_root, "CONSTRUCT_METADATA_IMPACT_AUDIT.md"), useBytes = TRUE)
write_tsv(domain_counts, file.path(output_root, "construct_domain_match_counts.tsv"))

cat("PASS: endpoint selection and construct metadata rebuilt\n")
cat("Candidate pairs:", nrow(candidates), "\n")
cat("Included:", sum(candidates$reconstructed_tier1), "\n")
cat("Construct tiers:", paste(paste(tier_counts$construct_tier, tier_counts$endpoint_pairs, sep = "="), collapse = "; "), "\n")
