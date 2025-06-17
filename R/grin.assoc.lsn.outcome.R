
#' Associate Lesions with Clinical Outcomes
#'
#' @description
#' Performs statistical association analysis between binary gene-lesion events and clinical outcomes of interest, including binary outcomes (e.g., Minimal Residual Disease (MRD)) and time-to-event outcomes (e.g., Event-Free Survival (EFS), and Overall Survival (OS)). The function supports both univariate and covariate-adjusted logistic regression and Cox proportional hazards models.
#'
#' @param lsn.mtx A binary lesion matrix where each row represents a unique gene-lesion pair (e.g., \code{ENSG00000148400_gain}). Each column represents a patient. Values are denoted as \code{1} if the patient harbors the specified lesion, and \code{0} otherwise. This matrix is typically produced using the \code{\link{prep.binary.lsn.mtx}} function.
#' @param clin.data A clinical data \code{data.frame}, where the first column \code{"ID"} represent patient identifiers matching those in \code{lsn.mtx}.
#' @param annotation.data A gene annotation \code{data.frame}, either provided by the user or retrieved using the \code{\link{get.ensembl.annotation}} function. Must include the columns: \code{"gene"} (Ensembl ID), \code{"chrom"} (chromosome), \code{"loc.start"} (gene start position), and \code{"loc.end"} (gene end position).
#' @param clinvars A character vector of clinical outcome variables to analyze. Binary variables (e.g., MRD) should be coded as \code{0} (negative) and \code{1} (positive). Survival outcomes (e.g., EFS, OS) must be precomputed using the \code{\link[survival]{Surv}} function and added as new columns to \code{clin.data}.
#' @param covariate Optional. A character vector specifying covariates to include as model adjustments (e.g., risk group, age, gender, etc...).
#'
#' @details
#' For each gene-lesion pair in the binary lesion matrix, the function can performs:
#' \itemize{
#'   \item \strong{Logistic regression} for binary outcomes (e.g., MRD), producing odds ratios (OR), 95_confidence intervals (CI), p-values, and FDR-adjusted q-values.
#'   \item \strong{Cox proportional hazards models} for survival outcomes (e.g., EFS, OS), producing hazard ratios (HR), 95\% CI, p-values, and FDR-adjusted q-values.
#' }
#' Models can optionally be adjusted for covariates such as clinical or demographic factors. Summary counts of patients with and without lesions, stratified by outcome status, are also included in the output.
#'
#' @return
#' A results \code{data.frame} containing gene annotation and association statistics for each gene-lesion pair across the specified clinical outcomes. The output includes:
#' \itemize{
#'   \item Odds ratio (OR), lower and upper 95CI, p-value, and q-value (FDR) for logistic regression models.
#'   \item Hazard ratio (HR), lower and upper 95CI, p-value, and q-value for Cox proportional hazards models.
#'   \item Patient counts for those with and without lesions, and corresponding outcome event statuses.
#' }
#'
#' @export
#'
#' @importFrom survival coxph is.Surv
#' @importFrom stats glm p.adjust
#' @importFrom stringr str_split_fixed
#' @importFrom utils head
#'
#' @references
#' Andersen, P. K., & Gill, R. D. (1982). Cox's regression model for counting processes: A large sample study.
#'
#' Therneau, T. M., & Grambsch, P. M. (2000). \emph{Modeling Survival Data: Extending the Cox Model}.
#'
#' Dobson, A. J. (1990). \emph{An Introduction to Generalized Linear Models}.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link{prep.binary.lsn.mtx}}, \code{\link[survival]{coxph}}, \code{\link[stats]{glm}}
#'
#' @examples
#' data(lesion_data)
#' data(hg38_gene_annotation)
#' data(clin_data)
#'
#' # Step 1: Prepare gene-lesion overlap
#' gene.lsn <- prep.gene.lsn.data(lesion_data, hg38_gene_annotation)
#' gene.lsn.overlap <- find.gene.lsn.overlaps(gene.lsn)
#'
#' # Step 2: Create a binary lesion matrix (minimum 5 patients per lesion)
#' lsn.binary.mtx <- prep.binary.lsn.mtx(gene.lsn.overlap, min.ngrp = 5)
#'
#' # Step 3: Create survival objects and add to clinical data
#' library(survival)
#' clin_data$EFS <- Surv(clin_data$efs.time, clin_data$efs.censor)
#' clin_data$OS <- Surv(clin_data$os.time, clin_data$os.censor)
#'
#' # Step 4: Specify outcomes of interest
#' clinvars <- c("MRD.binary", "EFS", "OS")
#'
#' # Step 5: Run association analysis
#' assc.outcomes <- grin.assoc.lsn.outcome(lsn.binary.mtx,
#'                                         clin_data,
#'                                         hg38_gene_annotation,
#'                                         clinvars)
#'
#' # Optional: Adjust for covariates using the 'covariate' argument

grin.assoc.lsn.outcome=function(lsn.mtx,            # output of prep.binary.lsn.mtx function, each lesion affecting a gene is represented in a separate row and patient will be coded as 1 if he/she is affected with the lesion and 0 otherwise
                                clin.data,          # clinical data table. First column "ID" has patient IDs and each row is a patient clinical data
                                annotation.data,    # gene annotation data file
                                clinvars,           # clinical variables (survival variables such as EFS and OS should be first coded as survival objects using Surv function and added as new columns to the clinical data file, binary variables such as MRD should be coded as 0, 1)
                                covariate=NULL)     # covariates that the model will adjust for if any
{
  lsn.mtx=t(lsn.mtx)
  lsn.df=as.data.frame(lsn.mtx)

  clin.data=clin.data
  clin.ids=clin.data$ID

  lsn.df=lsn.df[(rownames(lsn.df) %in% clin.ids),]
  lsn.ids=rownames(lsn.df)
  clin.data=clin.data[clin.data$ID %in% lsn.ids,]
  lsn.df=lsn.df[order(rownames(lsn.df)),]
  clin.data=clin.data[order(clin.data$ID),]
  # To make sure that both lesion and clinical data have same patient IDs
  check.rownames=all(rownames(lsn.df)==clin.data$ID)
  if (check.rownames==FALSE)
    stop("Gene-lesion matrix rownames must match patient IDs in the clinical data.")

  merged.data=cbind(lsn.df, clin.data)

  res<-data.frame(matrix(NA, length(lsn.df),))

  for (i in 1:length(clinvars))
  {
    var.name=clinvars[i]
    thisvar<- merged.data[, clinvars[i]]

    # If the variable is a survival object, run COXPH models
    if (survival::is.Surv(thisvar))
    {
      message(paste0("Running COXPH models for association with ", var.name, ":",date()))
      surv.time=thisvar[,1]
      surv.censor=thisvar[,2]
      surv.data=cbind(merged.data, surv.time, surv.censor)
      surv.data=surv.data[!is.na(surv.data$surv.time),]
      surv.censor=surv.data$surv.censor
      surv.lsn.clms=surv.data[ ,grepl("ENSG", names(surv.data))]

      surv.data.count=cbind(surv.censor, surv.lsn.clms)

      surv.event= surv.data.count[surv.data.count$surv.censor==1,]
      surv.without.event= surv.data.count[surv.data.count$surv.censor==0,]

      surv.event=surv.event[,-1]
      surv.without.event=surv.without.event[,-1]

      # count the number of event free and patients with events who are affected or not affected by genomic lesions
      surv.event.with.lsn=as.data.frame(colSums(surv.event== 1))
      surv.event.without.lsn=as.data.frame(colSums(surv.event== 0))
      surv.no.event.with.lsn=as.data.frame(colSums(surv.without.event== 1))
      surv.no.event.without.lsn=as.data.frame(colSums(surv.without.event== 0))

      surv.count=cbind(surv.event.with.lsn, surv.event.without.lsn,
                       surv.no.event.with.lsn, surv.no.event.without.lsn)

      colnames(surv.count)=c(paste0(var.name,".event.with.lsn"), paste0(var.name,".event.without.lsn"),
                             paste0(var.name,".no.event.with.lsn"), paste0(var.name,".no.event.without.lsn"))

      # Run COXPH models for association between genomic lesions and the survival object without adjustment for any covariate
      if (is.null(covariate)) {
        surv.model<-lapply(surv.lsn.clms,function(x) survival::coxph(survival::Surv(surv.time, surv.censor) ~ x , data=surv.data))

        ## To extract coefficients:
        surv.Coeff<-lapply(surv.model,function(f) summary(f)$coefficients[,1])
        ## For univariate mosels, number of columns should be 1
        surv.coeff.mat <- do.call(rbind,lapply(surv.Coeff,matrix,ncol=1,byrow=TRUE))

        ## To extract hazard.ratio:
        surv.hazard<-lapply(surv.model,function(f) summary(f)$coefficients[,2])
        surv.hazard.mat = as.matrix(surv.hazard)
        surv.hazard.mat=as.numeric(as.character(surv.hazard.mat))

        ## To extract model error:
        surv.error = lapply(surv.model,function(f) summary(f)$coefficients[,3])
        surv.error.mat <- do.call(rbind,lapply(surv.error,matrix,ncol=1,byrow=TRUE))

        ## To calculate lower95% CI:
        surv.lower95 = exp(surv.coeff.mat[,1] - 1.96*surv.error.mat[,1])
        surv.lower.mat = as.matrix(surv.lower95)
        surv.lower.mat=as.numeric(as.character(surv.lower.mat))

        ## To calculate upper95% CI:
        surv.upper95 = exp(surv.coeff.mat[,1] + 1.96*surv.error.mat[,1])
        surv.upper.mat = as.matrix(surv.upper95)
        surv.upper.mat=as.numeric(as.character(surv.upper.mat))

        ## To extract p-value:
        surv.Pvalue <-lapply(surv.model,function(f) summary(f)$coefficients[,5])
        surv.Pvalue.mat = as.matrix(surv.Pvalue)
        surv.pvalue.mat=do.call(rbind,lapply(surv.Pvalue,matrix,ncol=1,byrow=TRUE))
        pvalue.surv= surv.pvalue.mat[,1]
        pvalue.surv=as.numeric(as.character(pvalue.surv))

        # Compute FDR with the Pounds & Cheng (2006) estimator of the proportion of tests with a true null (pi.hat
        pi.hat=min(1,2*mean(pvalue.surv,na.rm=T))
        q.surv=pi.hat*stats::p.adjust(pvalue.surv,method="fdr")
        q.surv=as.numeric(as.character(q.surv))

        ## To prepare final results without covariates adjustment:
        surv.results=cbind(surv.hazard.mat, surv.lower.mat, surv.upper.mat, pvalue.surv, q.surv)
        surv.lsn=as.data.frame(colnames(lsn.mtx))
        colnames(surv.lsn)="Gene_lsn"
        results.surv.final=cbind(surv.lsn, surv.results)
        colnames(results.surv.final)=c("Gene_lsn",paste0(var.name, ".HR") ,
                                       paste0(var.name, ".lower95"), paste0(var.name, ".upper95"),
                                       paste0(var.name, ".p-value"), paste0(var.name, ".q-value.adj"))
        thisres=cbind(results.surv.final, surv.count)

      } else {

        # Run COXPH models for association between genomic lesions and survival object while adjusting for some covariates
        covariates=covariate
        covariates=surv.data[,covariates]

        surv.model.adj<-lapply(surv.lsn.clms,function(x) survival::coxph(survival::Surv(surv.time, surv.censor) ~ x+covariates , data=surv.data))

        ## To extract coefficients:
        surv.Coeff.adj<-lapply(surv.model.adj,function(f) summary(f)$coefficients[,1])
        ## Number of columns depend on how many variables we have in the model
        surv.coeff.mat.adj <- do.call(rbind,lapply(surv.Coeff.adj,matrix,byrow=TRUE))
        num.grps=length(surv.coeff.mat.adj)/ncol(lsn.mtx)
        surv.coeff.mat.adj <- do.call(rbind,lapply(surv.Coeff.adj,matrix,ncol=num.grps,byrow=TRUE))

        ## To extract hazard.ratio:
        surv.hazard.adj<-lapply(surv.model.adj,function(f) summary(f)$coefficients[,2])
        surv.hazard.mat.adj <- do.call(rbind,lapply(surv.hazard.adj,matrix,ncol=num.grps,byrow=TRUE))

        ## To extract model error:
        surv.error.adj = lapply(surv.model.adj,function(f) summary(f)$coefficients[,3])
        surv.error.mat.adj <- do.call(rbind,lapply(surv.error.adj,matrix,ncol=num.grps,byrow=TRUE))

        ## To calculate lower95% CI:
        surv.lower95.adj = exp(surv.coeff.mat.adj[,1] - 1.96*surv.error.mat.adj[,1])
        surv.lower.mat.adj = as.matrix(surv.lower95.adj)
        surv.lower.mat.adj=as.numeric(as.character(surv.lower.mat.adj))

        ## To calculate upper95% CI:
        surv.upper95.adj = exp(surv.coeff.mat.adj[,1] + 1.96*surv.error.mat.adj[,1])
        surv.upper.mat.adj = as.matrix(surv.upper95.adj)
        surv.upper.mat.adj=as.numeric(as.character(surv.upper.mat.adj))

        ## To extract p-value:
        surv.Pvalue.adj <-lapply(surv.model.adj,function(f) summary(f)$coefficients[,5])
        surv.Pvalue.mat.adj = as.matrix(surv.Pvalue.adj)
        surv.pvalue.mat.adj=do.call(rbind,lapply(surv.Pvalue.adj,matrix,ncol=num.grps,byrow=TRUE))
        pvalue.surv.adj=surv.pvalue.mat.adj[,1]
        pvalue.surv.adj=as.numeric(as.character(pvalue.surv.adj))

        # Compute FDR with the Pounds & Cheng (2006) estimator of the proportion of tests with a true null (pi.hat
        pi.hat=min(1,2*mean(pvalue.surv.adj,na.rm=T))
        q.surv.adj=pi.hat*stats::p.adjust(pvalue.surv.adj,method="fdr")
        q.surv.adj=as.numeric(as.character(q.surv.adj))

        ## To prepare final results (covariates adjusted):
        surv.results.adj=cbind(surv.hazard.mat.adj[,1], surv.lower.mat.adj, surv.upper.mat.adj,
                               pvalue.surv.adj, q.surv.adj)
        surv.lsn=as.data.frame(colnames(lsn.mtx))
        colnames(surv.lsn)="Gene_lsn"
        results.surv.final.adj=cbind(surv.lsn, surv.results.adj)
        colnames(results.surv.final.adj)=c("Gene_lsn", paste0(var.name, ".HR.adj"), paste0(var.name, ".lower95.adj"),
                                           paste0(var.name, ".upper95.adj"), paste0(var.name, ".p-value.adj"),
                                           paste0(var.name, ".q-value.adj"))
        thisres=cbind(results.surv.final.adj, surv.count)
      }
    }

    # if the variable is numeric such as MRD coded as 0 for MRD negative and 1 for MRD positive patients, run logistic regression models
    else if (is.numeric(thisvar))
    {
      message(paste0("Running logistic regression models for association with ", var.name, ":",date()))
      num.variable=thisvar
      num.data=cbind(merged.data, num.variable)
      num.data=num.data[!is.na(num.data$num.variable),]

      num.lsn.clms=num.data[ ,grepl("ENSG", names(num.data))]
      num.variable=num.data$num.variable

      num.data.count=cbind(num.variable, num.lsn.clms)

      num.event= num.data.count[num.data.count$num.variable==1,]
      num.without.event= num.data.count[num.data.count$num.variable==0,]

      num.event=num.event[,-1]
      num.without.event=num.without.event[,-1]

      # count the number of event free and patients with events who are affected or not affected by genomic lesions
      num.event.with.lsn=as.data.frame(colSums(num.event== 1))
      num.event.without.lsn=as.data.frame(colSums(num.event== 0))
      num.no.event.with.lsn=as.data.frame(colSums(num.without.event== 1))
      num.no.event.without.lsn=as.data.frame(colSums(num.without.event== 0))

      num.count=cbind(num.event.with.lsn, num.event.without.lsn,
                      num.no.event.with.lsn, num.no.event.without.lsn)

      colnames(num.count)=c(paste0(var.name,".event.with.lsn"), paste0(var.name,".event.without.lsn"),
                            paste0(var.name,".no.event.with.lsn"), paste0(var.name,".no.event.without.lsn"))


      # Run logistic regression models for association between genomic lesions and binary variables without adjustment for any covariate
      if (is.null(covariate)) {
        log.fit<-lapply(num.lsn.clms,function(x) stats::glm(num.variable  ~ x ,data = num.data, family = "binomial"))
        # To extract model coefficients
        coeff = lapply(log.fit,function(f) summary(f)$coefficients[,1])

        ## for univariate models, number of columns "ncol" should be 2 (column for intercept and column for model coefficient):
        coeff.mat <- do.call(rbind,lapply(coeff,matrix,ncol=2,byrow=TRUE))
        mod_coeff = lapply(coeff, utils::head, 2)

        ## To extract model error:
        error = lapply(log.fit,function(f) summary(f)$coefficients[,2])
        error.mat <- do.call(rbind,lapply(error,matrix,ncol=2,byrow=TRUE))

        ## To calculate odds ratio and 95% CI:
        ## N.B) Here, we are using the standard errors not the robust SE
        odds.ratio= exp(coeff.mat[,2])
        odds.mat=as.matrix(odds.ratio)
        odds.mat=as.numeric(as.character(odds.mat))

        lower.confint = exp(coeff.mat[,2] - 1.96*error.mat[,2])
        lower.mat=as.matrix(lower.confint)
        lower.mat=as.numeric(as.character(lower.mat))

        higher.confint = exp(coeff.mat[,2] + 1.96*error.mat[,2])
        higher.mat=as.matrix(higher.confint)
        higher.mat=as.numeric(as.character(higher.mat))

        ## To extract P-value:
        Pvalue<-lapply(log.fit,function(f) summary(f)$coefficients[,4])
        pvalue.mat=do.call(rbind,lapply(Pvalue,matrix,ncol=2,byrow=TRUE))
        pvalue.num= pvalue.mat[,2]
        pvalue.num=as.numeric(as.character(pvalue.num))

        # Compute FDR with the Pounds & Cheng (2006) estimator of the proportion of tests with a true null (pi.hat
        pi.hat=min(1,2*mean(pvalue.num,na.rm=T))
        q.num=pi.hat*stats::p.adjust(pvalue.num,method="fdr")
        q.num=as.numeric(as.character(q.num))

        # prepare final results without adjustment for any covariate
        results.logistic=cbind(odds.mat, lower.mat, higher.mat, pvalue.num, q.num)
        logistic.lsn=as.data.frame(colnames(lsn.mtx))
        colnames(logistic.lsn)="Gene_lsn"
        results.logistic.final=cbind(logistic.lsn, results.logistic)
        colnames(results.logistic.final)=c("Gene_lsn", paste0(var.name, ".odds.ratio"),
                                           paste0(var.name, ".lower95"), paste0(var.name, ".upper95"),
                                           paste0(var.name, ".p-value"), paste0(var.name, ".q-value"))
        num.results.count.added=cbind(results.logistic.final, num.count)
        thisres=num.results.count.added

      } else {

        # Run logistic regression models for association between genomic lesions and binary variable with covariate adjustment
        covariates=covariate
        covariates=num.data[,covariates]
        log.fit.adj<-lapply(num.lsn.clms,function(x) stats::glm(data = num.data, num.variable  ~ x + covariates , family = "binomial"))
        ## To extract coefficients
        coeff.adj = lapply(log.fit.adj,function(f) summary(f)$coefficients[,1])
        coeff.mat.adj <- do.call(rbind,lapply(coeff.adj,matrix,byrow=TRUE))
        num.grps=round(length(coeff.mat.adj)/ncol(lsn.mtx))
        coeff.mat.adj <- do.call(rbind,lapply(coeff.adj,matrix,ncol=num.grps,byrow=TRUE))

        ## To extract model error:
        error.adj = lapply(log.fit.adj,function(f) summary(f)$coefficients[,2])
        error.mat.adj <- do.call(rbind,lapply(error.adj,matrix,ncol=num.grps,byrow=TRUE))

        ## To calculate odds ratio and 95% CI:
        ## N.B) Here, we are using the standard errors not the robust SE
        odds.ratio.adj= exp(coeff.mat.adj[,2])
        odds.mat.adj=as.matrix(odds.ratio.adj)
        odds.mat.adj=as.numeric(as.character(odds.mat.adj))

        lower.confint.adj = exp(coeff.mat.adj[,2] - 1.96*error.mat.adj[,2])
        lower.mat.adj=as.matrix(lower.confint.adj)
        lower.mat.adj=as.numeric(as.character(lower.mat.adj))

        higher.confint.adj = exp(coeff.mat.adj[,2] + 1.96*error.mat.adj[,2])
        higher.mat.adj=as.matrix(higher.confint.adj)
        higher.mat.adj=as.numeric(as.character(higher.mat.adj))

        ## To extract P-value:
        Pvalue.adj<-lapply(log.fit.adj,function(f) summary(f)$coefficients[,4])
        pvalue.mat.adj=do.call(rbind,lapply(Pvalue.adj,matrix,ncol=num.grps,byrow=TRUE))
        pvalue.num.adj= pvalue.mat.adj[,2]
        pvalue.num.adj=as.numeric(as.character(pvalue.num.adj))

        # Compute FDR with the Pounds & Cheng (2006) estimator of the proportion of tests with a true null (pi.hat
        pi.hat=min(1,2*mean(pvalue.num.adj,na.rm=T))
        q.num.adj=pi.hat*stats::p.adjust(pvalue.num.adj,method="fdr")
        q.num.adj=as.numeric(as.character(q.num.adj))

        # prepare final results while adjsuting for covariates
        results.logistic.adj=cbind(odds.mat.adj, lower.mat.adj, higher.mat.adj, pvalue.num.adj, q.num.adj)
        logistic.lsn=as.data.frame(colnames(lsn.mtx))
        colnames(logistic.lsn)="Gene_lsn"
        results.logistic.final.adj=cbind(logistic.lsn, results.logistic.adj)
        colnames(results.logistic.final.adj)=c("Gene_lsn", paste0(var.name, ".odds.ratio.adj"),
                                               paste0(var.name, ".lower95.adj"), paste0(var.name, ".upper95.adj"),
                                               paste0(var.name, ".p-value.adj"), paste0(var.name, ".q-value.adj"))
        num.results.count.added.adj=cbind(results.logistic.final.adj, num.count)
        thisres=num.results.count.added.adj

      }
    }
    res<-cbind.data.frame(res, thisres)
    res$matrix.NA..length.lsn.df....<-NULL

    gene.list=stringr::str_split_fixed(res$Gene_lsn, "_", 2)
    colnames(gene.list)=c("gene","lsn")
    annotation.merged=merge(annotation.data,gene.list,by="gene", all.y=TRUE)
    res.final=cbind(annotation.merged, res)
    res.final$lsn<-NULL
  }

  return(res.final)
}
