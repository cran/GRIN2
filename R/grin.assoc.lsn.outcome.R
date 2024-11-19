
#' Associate Lesions with Clinical Outcomes
#'
#' @description
#' The function run association analysis between the binary lesion matrix (output of prep.binary.lsn.mtx function) and clinical outcomes of interest such as Minimal Residual Disease (MRD), Event-free Survival (EFS) and Overall Survival (OS), etc...
#'
#' @param lsn.mtx Binary lesion matrix in which each type of lesions affecting certain gene is represented in a separate row for example ENSG00000148400_gain. If the gene is affected by this specific type of lesion, patient entry will be coded as 1 or 0 otherwisw. This matrix is the output of the prep.binary.lsn.mtx function.
#' @param clin.data Clinical data table in which the first column "ID" should has the patient ID.
#' @param annotation.data Gene annotation data either provided by the user or retrieved from ensembl BioMart database using get.ensembl.annotation function included in the GRIN2.0 library. Data.frame should has four columns: "gene" which is the ensembl ID of annotated genes, "chrom" which is the chromosome on which the gene is located, "loc.start" which is the gene start position, and "loc.end" the gene end position.
#' @param clinvars Clinical outcome variables of interest (survival variables such as EFS and OS should be first coded as survival objects using Surv function and added as new columns to the clinical data file, binary variables such as MRD should be coded as 0, 1).
#' @param covariate Covariates that the model will adjust for if any.
#'
#' @details
#' The function run association analysis between the binary lesion matrix in which each type of lesions affecting certain gene is represented in a separate row (output of prep.binary.lsn.mtx function) and clinical outcomes.Function will run logistic regression models for association between each gene-lesion pair with numeric variables such as MRD that should be coded as 0 if the patient is MRD-negative and 1 if the patient is MRD positive. Function will also run COX-Proportional hazard models for association between lesions and survival objects such as Event-free survival (EFS) and oveall survival (OS). EFS and OS should be first coded as survival objects using Surv function and added as new columns to the clinical data file. If specified, the models can be also adjusted for one or a group of covariates such as risk group assignment, gender, age, etc...
#'
#' @return
#' Function returns a results table that has gene annotation data, and multiple columns showing results of the logistic regression model for association with binary variables such as MRD that include odds.ratio, lower and upper 95 confidence interval (CI), model p and FDR adjusted q values, in addition to the number of patients with/without lesion who experienced or did not experience the event. Results table will also include results of COXPH models for association between lesions with survival variables such as EFS, OS that include COXPH hazard ratio, lower and upper 95 CI, model p and FDR adjusted q values, in addition to the number of patients with/without the lesion who experienced or did not experience the event.
#'
#' @export
#'
#' @importFrom survival coxph is.Surv
#' @importFrom stats glm p.adjust
#' @importFrom stringr str_split_fixed
#' @importFrom utils head
#'
#' @references
#' Andersen, P. and Gill, R. (1982). Cox's regression model for counting processes, a large sample study.
#'
#' Therneau, T., Grambsch, P. (2000) Modeling Survival Data: Extending the Cox Model.
#'
#' Dobson, A. J. (1990) An Introduction to Generalized Linear Models.
#'
#' @author {Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [stats::glm()], [survival::coxph()]
#'
#' @examples
#' data(lesion.data)
#' data(hg19.gene.annotation)
#' data(clin.data)
#'
#' # prepare lesion data and find gene lesion overlaps:
#' gene.lsn=prep.gene.lsn.data(lesion.data, hg19.gene.annotation)
#' gene.lsn.overlap= find.gene.lsn.overlaps(gene.lsn)
#'
#' # Prepare a binary lesion matrix for genes affected by a certain type of lesion in at least
#' # 5 subjects using prep.binary.lsn.mtx function:
#' lsn.binary.mtx=prep.binary.lsn.mtx(gene.lsn.overlap, min.ngrp=5)
#'
#' # Prepare EFS and OS survival objects and add two new columns to the clinical data file:
#' library(survival)
#' clin.data$EFS <- Surv(clin.data$efs.time, clin.data$efs.censor)
#' clin.data$OS <- Surv(clin.data$os.time, clin.data$os.censor)
#' # define clinical outcome variables to be included in the analysis:
#' clinvars=c("MRD.binary", "EFS", "OS")
#'
#' # Run association analysis between lesions in the binary lesion matrix and clinical variables
#' # in the clinvars object:
#' assc.outcomes=grin.assoc.lsn.outcome(lsn.binary.mtx,
#'                                      clin.data,
#'                                      hg19.gene.annotation,
#'                                      clinvars)
#'
#' # to adjust the models for one or a group of covariates, user can specify one or a group
#' # of covariates using the 'covariate' argument.

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
