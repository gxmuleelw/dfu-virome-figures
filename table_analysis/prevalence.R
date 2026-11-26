prevalence_filter <- function(counts, detection_count, prevalence_cut, label) {
  n_prevalence_input <- ncol(counts)
  threshold <- ceiling(ncol(counts) * prevalence_cut)
  keep <- rowSums(counts >= detection_count) >= threshold
  filter <- data.frame(scenario = label, vOTU = rownames(counts),
                       detection_count = detection_count, prevalence_cut = prevalence_cut,
                       n_prevalence_input = n_prevalence_input,
                       n_meeting_threshold = rowSums(counts >= detection_count),
                       required_n = threshold, retained = keep,
                       reason = ifelse(keep, "included", "below_prevalence_cut"))
  filter
}
