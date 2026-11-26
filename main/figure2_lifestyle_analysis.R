# Figure 2: percentages are calculated within each patient before the group comparison.
# Input: breadth-filtered RPKM (vOTUs x patients), a named category per vOTU,
# and patient metadata (sample, severity, Age, Sex). No study data are included.

source(local({
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
  current <- if (length(files)) tail(files,1L)[[1L]] else
    sub("^--file=", "", grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
  file.path(dirname(normalizePath(current)), "hc3.R")
}), local=TRUE)

safe_lifestyle_fit<-function(z){
 z<-z[complete.cases(z[,c('value','severity','Age','Sex')]),];positive<-tapply(z$value>0,z$severity,sum)
 status<-if(length(unique(z$value))<=1L)'descriptive_constant'else if(length(positive)!=2L||any(is.na(positive)|positive<5L))'descriptive_sparse'else'estimated'
 fit<-data.frame(estimate=NA_real_,se=NA_real_,lower=NA_real_,upper=NA_real_,p=NA_real_,n=nrow(z),df=NA_real_)
 if(status=='estimated'){attempt<-tryCatch(fit_hc3(value~severity+Age+Sex,z,'severityWagner_4'),error=identity);if(inherits(attempt,'error'))status<-'skipped_model_error'else fit<-attempt}
 cbind(data.frame(status),fit)
}

# One scenario at a time: pass excluded vOTU IDs for a catalogue exclusion.
# Without exclusions, pass all directly lytic-labelled IDs for the broader-label
# sensitivity and/or the original cohort's E_coli_z values in metadata.
lifestyle_analysis <- function(rpkm, category, meta, exclude=character(), direct_lytic=NULL) {
  classes <- c("lytic", "temperate", "unclassified")
  ids <- rownames(rpkm); patients <- colnames(rpkm)
  stopifnot(is.matrix(rpkm), !is.null(ids), !is.null(patients),
            !anyDuplicated(ids), !anyDuplicated(patients), !anyDuplicated(meta$sample),
            !anyDuplicated(names(category)), all(ids %in% names(category)),
            all(patients %in% meta$sample), all(category[ids] %in% classes),
            all(is.finite(rpkm)), all(rpkm >= 0), all(exclude %in% ids),
            all(direct_lytic %in% ids))
  kept <- setdiff(ids, exclude)
  m <- rpkm[kept,,drop=FALSE]
  totals <- colSums(m); valid <- totals > 0
  stopifnot(any(valid))
  rel <- sweep(m[,valid,drop=FALSE], 2, totals[valid], "/")
  d <- meta[match(colnames(rel), meta$sample),,drop=FALSE]
  d$severity <- factor(d$severity, levels=c("Wagner_1_3","Wagner_4"))
  d$Sex <- factor(d$Sex, levels=c("Female","Male"))
  pct <- function(members) 100 * colSums(rel[intersect(kept,members),,drop=FALSE])
  patient_pct <- vapply(classes,function(cl) pct(ids[category[ids]==cl]),
                        numeric(ncol(rel)))
  dim(patient_pct) <- c(ncol(rel),length(classes)); dimnames(patient_pct) <- list(colnames(rel),classes)
  stopifnot(all(abs(rowSums(patient_pct)-100)<1e-8))
  fit_value <- function(y) { d$value <- y; safe_lifestyle_fit(d) }
  fits <- do.call(rbind,lapply(classes,function(cl) fit_value(patient_pct[,cl])))
  fits$category <- classes
  fits$q <- p.adjust(fits$p,method="BH",n=length(classes))
  extra <- list()
  if (!is.null(direct_lytic)) {
    stopifnot(length(exclude)==0L)
    extra$broader_label <- fit_value(pct(direct_lytic))
  }
  if ("E_coli_z" %in% names(d)) {
    d$value <- patient_pct[,"lytic"]
    extra$ecoli_adjusted <- fit_hc3(value~severity+Age+Sex+E_coli_z,d,"severityWagner_4")
  }
  # Extra fits are separate contrasts with no BH q in the figure.
  list(patient_percent=data.frame(sample=d$sample,patient_pct),models=fits,
       extra=extra,zero_denominator=patients[!valid])
}

# Panel d uses the saved best-MAG table. Coverage is in percent. The denominator
# includes every member of its category, even with a missing genus match.
same_genus_coverage <- function(category, genus_match, coverage_pct) {
  stopifnot(length(category)==length(genus_match),length(category)==length(coverage_pct),
            is.logical(genus_match),!anyNA(category))
  hit <- !is.na(genus_match) & genus_match & !is.na(coverage_pct) & coverage_pct>50
  do.call(rbind,lapply(unique(category),function(cl) {
    ii <- category==cl; k <- sum(hit[ii]); n <- sum(ii)
    data.frame(category=cl,count=k,denominator=n,percent=100*k/n)
  }))
}

# Panel e: reproduce the five displayed contrasts from already annotated inputs.
# catalog columns: vOTU, category, host_genus, host_genes, viral_genes, direct_lytic
# (the last is the nonexclusive exact label, not the three-way category).
figure2_sensitivity <- function(rpkm,catalog,meta) {
  stopifnot(!anyDuplicated(catalog$vOTU),setequal(rownames(rpkm),catalog$vOTU),
            is.logical(catalog$direct_lytic),!anyNA(catalog$direct_lytic),
            !anyNA(catalog$host_genes),!anyNA(catalog$viral_genes),
            "E_coli_z" %in% names(meta))
  classes <- setNames(catalog$category,catalog$vOTU)
  ordinary <- lifestyle_analysis(rpkm,classes,meta,direct_lytic=catalog$vOTU[catalog$direct_lytic])
  unadjusted <- meta[,setdiff(names(meta),"E_coli_z"),drop=FALSE]
  without_host <- lifestyle_analysis(rpkm,classes,unadjusted,
                    exclude=catalog$vOTU[which(catalog$host_genus=="Escherichia")])
  without_markers <- lifestyle_analysis(rpkm,classes,unadjusted,
                    exclude=catalog$vOTU[catalog$host_genes>catalog$viral_genes])
  extract <- function(result) subset(result$models,category=="lytic")
  rows <- list(original=extract(ordinary), remove_host=extract(without_host),
               broader_label=ordinary$extra$broader_label,
               ecoli_adjusted=ordinary$extra$ecoli_adjusted,
               remove_host_gt_viral=extract(without_markers))
  do.call(rbind,lapply(names(rows),function(s) {
    r <- rows[[s]]
    data.frame(scenario=s,r[,c("estimate","lower","upper","p","n","df")],
               q=if ("q" %in% names(r)) r$q else NA_real_,row.names=NULL)
  }))
}
