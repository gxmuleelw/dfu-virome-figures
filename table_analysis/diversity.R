diversity_values <- function(p) {
  plogp <- p
  positive <- p > 0
  plogp[positive] <- p[positive] * log(p[positive])
  rbind(Shannon = -colSums(plogp), inverse_Simpson = 1 / colSums(p^2))
}
alpha_values <- function(rel, counts) {
  stopifnot(identical(dimnames(rel), dimnames(counts)),
            all(is.finite(rel)), all(rel >= 0),
            all(abs(colSums(rel) - 1) < 1e-8), all(colSums(counts) > 0))
  values <- diversity_values(rel)
  data.frame(sample = colnames(rel),
             Shannon = values["Shannon", ],
             inverse_Simpson = values["inverse_Simpson", ],
             observed_mapped_votu = colSums(counts > 0), row.names = NULL)
}
