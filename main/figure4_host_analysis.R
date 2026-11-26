# Figure 4: compute ranks for panels a-b; model detection and positive abundance
# (panels c-e) from filtered counts/RPKM and prior MAG features.
source(local({
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
  current <- if (length(files)) tail(files,1L)[[1L]] else
    sub("^--file=", "", grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
  file.path(dirname(normalizePath(current)), "hc3.R")
}), local=TRUE)

# The input matrix contains ALL MAGs (rows) and the full cohort (columns).
# Select one E. coli and two previously specified Collinsella MAGs in that order.
prior_mag_score <- function(mag_counts, prior_mags) {
  stopifnot(is.matrix(mag_counts), length(prior_mags)==3L, !anyDuplicated(prior_mags),
            !is.null(colnames(mag_counts)), !anyDuplicated(colnames(mag_counts)),
            !anyDuplicated(rownames(mag_counts)), all(prior_mags %in% rownames(mag_counts)),
            all(is.finite(mag_counts)), all(mag_counts>=0), all(colSums(mag_counts)>0))
  pct <- 100*t(sweep(mag_counts[prior_mags,,drop=FALSE],2,colSums(mag_counts),"/"))
  stopifnot(any(pct>0))
  pseudo <- min(pct[pct>0])/2
  log_pct <- log10(pct+pseudo)
  stopifnot(all(apply(log_pct,2,sd)>0))
  z <- scale(log_pct)
  data.frame(sample=rownames(pct), E_coli_z=z[,1], Collinsella_1_z=z[,2],
             Collinsella_2_z=z[,3], bacterial_score=z[,1]-(z[,2]+z[,3])/2,
             row.names=NULL)
}

# One genus per call: include zeros, omit unavailable compositions, rank across
# both Wagner groups. The fitted line is descriptive; it supplies no test.
paired_ranks <- function(paired) {
  d <- paired[which(paired$composition_available),,drop=FALSE]
  stopifnot(nrow(d)>=2L,!anyDuplicated(d$sample),all(is.finite(d$bacterial_relative)),
            all(is.finite(d$viral_relative)),length(unique(d$bacterial_relative))>1L)
  d$bacterial_rank <- rank(d$bacterial_relative,ties.method="average")
  d$viral_rank <- rank(d$viral_relative,ties.method="average")
  fit <- lm(viral_rank~bacterial_rank,d)
  list(ranks=d,descriptive_line=coef(fit),x_span=range(d$bacterial_rank))
}

# Standalone: Rscript figure4_host_analysis.R EXTERNAL_WORKDIR
# External tables: rpkm.tsv and counts.tsv (features x samples), metadata.tsv,
# bacteria.tsv (prior_mag_score output), hosts.tsv (host assignments).
run_host_models <- function(out) {
  previous <- options(digits=17,stringsAsFactors=FALSE)
  on.exit(options(previous),add=TRUE)
  stopifnot(dir.exists(out))
  dir.create(file.path(out,"models"),showWarnings=FALSE)
  read_tsv<-function(f)read.delim(file.path(out,f),check.names=FALSE,quote='',comment.char='',fileEncoding='UTF-8')
  write_tsv<-function(d,f)write.table(d,file.path(out,f),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
  getmat<-function(f){d<-read_tsv(f);m<-as.matrix(d[,-1]);rownames(m)<-d[[1]];storage.mode(m)<-'double';m}
  mat<-getmat('rpkm.tsv');counts<-getmat('counts.tsv');meta<-read_tsv('metadata.tsv');bac<-read_tsv('bacteria.tsv');hosts<-read_tsv('hosts.tsv')
  stopifnot(identical(dimnames(mat),dimnames(counts)), !anyDuplicated(rownames(mat)), !anyDuplicated(colnames(mat)),
   !anyDuplicated(meta$sample), !anyDuplicated(bac$sample), !anyDuplicated(hosts$vOTU), all(colnames(mat) %in% meta$sample),
   all(colnames(mat) %in% bac$sample), all(rownames(mat) %in% hosts$vOTU), all(is.finite(mat)), all(mat>=0), all(is.finite(counts)), all(counts>=0))
  meta<-meta[match(colnames(mat),meta$sample),];bac<-bac[match(meta$sample,bac$sample),]
  meta$severity<-factor(meta$severity,levels=c('Wagner_1_3','Wagner_4'));meta$Sex<-factor(meta$Sex,levels=c('Female','Male'))
  total<-colSums(mat);keep<-total>0;stopifnot(any(keep))
  rel<-sweep(mat[,keep,drop=FALSE],2,total[keep],'/')
  ids<-sort(hosts$vOTU[hosts$host_group=='Escherichia']);stopifnot(length(ids)>0)
  diagnostics<-list()
  family_n<-8L
  result_row <- function(variant,host,test) data.frame(variant,host,test,status='not_estimable',reason='',
    method='',formula='',term='',effect_scale='',estimate=NA_real_,se=NA_real_,lower=NA_real_,
    upper=NA_real_,p=NA_real_,q=NA_real_,n=0L,n_ref=0L,n_case=0L,positive_ref=0L,positive_case=0L,
    df=NA_real_,family_n=family_n)
  continuous_fit <- function(d,variant,host,test,formula,term) {
    r <- result_row(variant,host,test)
    r$method <- 'OLS HC3; residual-df t inference'
    r$formula <- paste(deparse(formula),collapse=' '); r$term <- term
    r$effect_scale <- if(term=='bacterial_score') 'log2 viral fraction per 1 original bacterial score unit' else 'log2 viral fraction; Wagner4 minus Wagner1-3'
    d <- d[complete.cases(d[,all.vars(formula),drop=FALSE]),]
    r$n <- nrow(d); r$n_ref <- sum(d$severity=='Wagner_1_3'); r$n_case <- sum(d$severity=='Wagner_4')
    r$positive_ref <- r$n_ref; r$positive_case <- r$n_case
    if(min(r$n_ref,r$n_case)<5L) {r$reason <- 'fewer_than_5_positive_complete_in_a_severity_group';return(r)}
    if(r$n<20L) {r$reason <- 'fewer_than_20_positive_complete';return(r)}
    fit <- tryCatch(lm(formula,d,na.action=na.fail),error=identity)
    if(inherits(fit,'error')) {r$reason <- paste0('model_error:',conditionMessage(fit));return(r)}
    if(fit$rank<ncol(model.matrix(fit))) {r$reason <- 'rank_deficient_design';return(r)}
    if(df.residual(fit)<10L) {r$reason <- 'fewer_than_10_residual_df';return(r)}
    ans <- tryCatch(hc3_result(fit,term),error=identity)
    if(inherits(ans,'error') || !all(is.finite(unlist(ans))) || ans$se<=0) {
      r$reason <- 'nonfinite_or_nonpositive_HC3_inference';return(r)
    }
    r[names(ans)] <- ans; r$status <- 'estimated'; r$reason <- ''
    key <- paste(variant,host,test,sep='__')
    saveRDS(list(fit=fit,target_term=term,covariance=sandwich::vcovHC(fit,type='HC3'),data=d,
      standardization='Bacterial variables frozen on original cohort; viral response not standardized'),
      file.path(out,'models',paste0(key,'.rds')))
    diagnostics[[key]] <<- data.frame(variant,host,test,sample=d$sample,
      fitted=fitted(fit),residual=residuals(fit),leverage=hatvalues(fit),cooks_distance=cooks.distance(fit))
    r
  }
  sets<-list(Escherichia=ids,Collinsella=sort(hosts$vOTU[hosts$host_group=='Collinsella']))

  results<-list();metrics<-list();members<-list()
  for(s in names(sets)) {
   vid<-sets[[s]];variant<-'score90';host<-s
   members[[s]]<-data.frame(variant,host,vOTU=vid,confidence=hosts$best_score[match(vid,hosts$vOTU)],n_candidates=hosts$n_candidates[match(vid,hosts$vOTU)])
   d<-meta[,c('sample','severity','Age','Sex')];d$host_z<-if(host=='Escherichia') bac$E_coli_z else (bac$Collinsella_1_z+bac$Collinsella_2_z)/2
   d$bacterial_score<-bac$bacterial_score;d$detected<-colSums(counts[vid,,drop=FALSE])>0
   d$viral_relative<-NA_real_;d$viral_relative[keep]<-colSums(rel[vid,,drop=FALSE]);d$positive_abundance_available<-is.finite(d$viral_relative)&d$viral_relative>0
   d$viral_log2_positive<-NA_real_;d$viral_log2_positive[d$positive_abundance_available]<-log2(d$viral_relative[d$positive_abundance_available])
   r<-result_row(variant,host,'detection_severity');r$method<-'Fisher exact two-sided; no covariate adjustment';r$formula<-'detected by severity';r$term<-'severityWagner_4';r$effect_scale<-'odds ratio; Wagner4 versus Wagner1-3'
   r$n<-nrow(d);r$n_ref<-sum(d$severity=='Wagner_1_3');r$n_case<-sum(d$severity=='Wagner_4');r$positive_ref<-sum(d$detected&d$severity=='Wagner_1_3');r$positive_case<-sum(d$detected&d$severity=='Wagner_4')
   if(length(unique(d$detected))<2L) r$reason<-'constant_detection' else {
    contingency<-matrix(c(r$positive_case,r$n_case-r$positive_case,r$positive_ref,r$n_ref-r$positive_ref),nrow=2,byrow=TRUE)
    f<-fisher.test(contingency,alternative='two.sided',conf.level=.95);r$estimate<-unname(f$estimate);r$lower<-f$conf.int[1];r$upper<-f$conf.int[2];r$p<-f$p.value;r$status<-'estimated'
   }
   results[[length(results)+1]]<-r
   pos<-d[d$positive_abundance_available,]
   results[[length(results)+1]]<-continuous_fit(pos,variant,host,'positive_abundance_severity',viral_log2_positive~severity+Age+Sex,'severityWagner_4')
   results[[length(results)+1]]<-continuous_fit(pos,variant,host,'positive_abundance_score',viral_log2_positive~bacterial_score+severity+Age+Sex,'bacterial_score')
   results[[length(results)+1]]<-continuous_fit(pos,variant,host,'host_adjusted_severity',viral_log2_positive~host_z+severity+Age+Sex,'severityWagner_4')
   d$variant<-variant;d$host<-host;metrics[[s]]<-d
  }
  res<-do.call(rbind,results);valid<-is.finite(res$p);res$q[valid]<-p.adjust(res$p[valid],method='BH',n=family_n)
  res$raw_p_candidate<-is.finite(res$p)&res$p<.05;res$BH_candidate<-is.finite(res$q)&res$q<.05
  write_tsv(res,'primary_tests.tsv');write_tsv(do.call(rbind,members),'host_set_membership.tsv');write_tsv(do.call(rbind,metrics),'sample_metrics.tsv');write_tsv(do.call(rbind,diagnostics),'model_diagnostics.tsv')
}
if (sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  stopifnot(length(args)==1L)
  run_host_models(args[[1L]])
}
