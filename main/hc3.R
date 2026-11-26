# OLS coefficient inference shared by the lifestyle and focal host analyses.
hc3_result<-function(fit,term){
 stopifnot(fit$rank==ncol(model.matrix(fit)),df.residual(fit)>0,term %in% names(coef(fit)))
 se<-sqrt(diag(sandwich::vcovHC(fit,type='HC3')))[term];estimate<-unname(coef(fit)[term]);df<-df.residual(fit)
 data.frame(estimate=estimate,se=unname(se),lower=estimate-qt(.975,df)*se,upper=estimate+qt(.975,df)*se,p=2*pt(-abs(estimate/se),df),n=nobs(fit),df=df,row.names=NULL)
}
fit_hc3<-function(formula,data,term)hc3_result(lm(formula,data=data,na.action=na.omit),term)
