
#' Get Ensembl Gene and Regulatory Features Annotation Data
#'
#' @description
#' Function directly retrieve gene and regulatory features annotation data from Ensembl BioMart database based on the specified genome assembly.
#'
#' @param genome.assembly User can specify one of four genome assemblies that include "Human_GRCh38", "Human_GRCh37", "Mouse_HGCm39" and "Mouse_HGCm38".
#'
#' @details
#' Based on the genome assembly specified by the user, the function will directly retrieve gene and regulatory features annotation data from ensembl BioMart database. Annotation data include enesembl ID, the chromosome on which the gene is located, gene start and gene end position, gene name, gene description, biotype, chromosome strand and chromosome band. Gene classes (biotypes) include protein coding genes, long noncoding RNAs (lncRNAs), microRNAs (miRNAs), small nuclear RNAs (snRNA), small nucleolar RNAs (snoRNA), immunoglobulins (IGs), T-cell receptors (TCRs) and pseudogens. Regulatory features data retrieved from Ensembl regulatory build are categorized in 6 classes that include promoters, promoter flanking regions, predicted enhancers, CTCF binding sites, transcription factor (TF) binding sites and the open chromatin regions.Ensembl first imports publicly available data from different large epigenomic consortia such as ENCODE, Roadmap Epigenomics and Blueprint. All high-throughput sequencing data sets are then uniformly processed using the Ensembl Regulation Sequence Analysis (ERSA) pipeline to generate signal tracks for enriched regions also referred to as annotated features or peaks. Segmentation data provide information about promoter, promoter flanking regions, enhancers and CTCF binding sites. If TF binding probability is >0 in areas outside previously mentioned regions, it takes a label of TF binding site. If any open chromatin region did not overlap with the above features, it takes a label of unannotated open chromatin. users will also have the chance to use a list of experimentally verified enhancers/transcription start sites (TSS) using the CAGE (Cap Analysis of Gene Expression) experiment on a multitude of different primary cells and tissues from the Functional Annotation of the Mouse/Mammalian Genome (FANTOM5) project.
#'
#' @return
#' A list of three components:
#' \item{gene.annotation}{A 9 columns data.frame with gene annotation data that include enesembl ID, chromosome, gene start and gene end position, gene name, gene description, biotype, chromosome strand, and chromosome band.}
#' \item{reg.annotation.predicted}{A 5 columns data.frame with regulatory features annotation data directly retreived from Ensembl regulatory build that include enesembl ID, chromosome, description(promoter, enhancer, etc..), feature start and end positions.}
#' \item{reg.annotation.validated}{A 5 columns data.frame with regulatory features annotation data for experimentally verified features retreived from FANTOM5 project that include feature ID, chromosome, description(enhancer, transcription start site (TSS)), feature start and end positions.}
#'
#' @export
#' @importFrom biomaRt useEnsembl getBM
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' Zerbino, Daniel R., et al. (2015). The ensembl regulatory build.
#'
#' Kinsella, Rhoda J., et al. (2011). Ensembl BioMarts: a hub for data retrieval across taxonomic space.
#'
#' @author {Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [biomaRt::useEnsembl()], [biomaRt::getBM()]
#'
#' @examples
#' \donttest{
#' # Retrieve annotation data for human hg19 genome assembly:
#' hg19.ann=get.ensembl.annotation("Human_GRCh37")
#' # gene annotation data:
#' hg19.gene.annotation=hg19.ann$gene.annotation
#' # regulatory features annotation data retrieved from Ensembl regulatory build:
#' hg19.reg.annotation=hg19.ann$reg.annotation.predicted
#' # regulatory features annotation data retrieved from FANTOM5 project:
#' hg19.fantom.annotation=hg19.ann$reg.annotation.validated
#'
#' # "Human_GRCh38" can be used instead of "Human_GRCh37" to retrieve gene and regulatory features
#' # annotation data for human hg38 genome assembly.
#' }

get.ensembl.annotation=function(genome.assembly)  # one of four genome assemblies that include "Human_GRCh38", "Human_GRCh37", "Mouse_HGCm39" and "Mouse_HGCm38" can be specified

{
  # retrieve gene annotation data for human_GRCh38 (hg38) genome assembly from ensembl biomaRt version 110
  if (genome.assembly=="Human_GRCh38")
  {
    message(paste0("Retrieving gene annotation data from Ensembl BioMart: ",date()))
    ensembl_GRCh38 = biomaRt::useEnsembl(biomart="genes",
                                         dataset="hsapiens_gene_ensembl",  # specify dataset for homosapiens
                                         version = "110")                  # specifying version is critical to get stable search results. If we did not specify version, query will extract data from the most updated version

    chromosomes = c(1:22, "X", "Y")                               # specify chromosomes of interest
    # specify data elements to retreive for each gene
    hg38_gene_annotation <- biomaRt::getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                                        'start_position','end_position',
                                                        'description','external_gene_name',
                                                        'gene_biotype', 'strand', 'band'),
                                           filters = 'chromosome_name',  values = chromosomes,
                                           mart = ensembl_GRCh38)         # attributes specify data to be retrieved from biomaRt database

    gene=hg38_gene_annotation[,1]
    chrom=hg38_gene_annotation[,2]
    loc.start=hg38_gene_annotation[,3]
    loc.end=hg38_gene_annotation[,4]
    description=hg38_gene_annotation[,5]
    gene.name=hg38_gene_annotation[,6]
    biotype=hg38_gene_annotation[,7]
    chrom.strand=hg38_gene_annotation[,8]
    chrom.band=hg38_gene_annotation[,9]

    gene.data=cbind.data.frame(gene=gene,
                               chrom=chrom,
                               loc.start=loc.start,
                               loc.end=loc.end,
                               description=description,
                               gene.name=gene.name,
                               biotype=biotype,
                               chrom.strand=chrom.strand,
                               chrom.band=chrom.band)

    # To retrieve data for computationally predicted regulatory features mapped to GRCh38 from ensembl regulatory build
    message(paste0("Retrieving computationally predicted regulatory features from Ensembl regulatory build: ",date()))

    regulatory.hg38 = biomaRt::useEnsembl(biomart="regulation",
                                          dataset="hsapiens_regulatory_feature", # regulatory features includes (promoters, promoter flanking regions, enhancers, CTCF binding sites, TF binding sites and open chromatin regions)
                                          version = '110')                       # specify version for stable results

    # specify classes of regulatory features of interest
    regulatory.features=c("Promoter", "Promoter Flanking Region", "Enhancer",
                          "CTCF Binding Site", "TF binding site", "Open chromatin")

    # specify data elements to retreive for each regulatory feature
    hg38.regulatory= biomaRt::getBM(attributes=c("regulatory_stable_id", "chromosome_name",
                                                 "chromosome_start","chromosome_end",
                                                 "feature_type_description"),
                                    filters =c("regulatory_feature_type_name", "chromosome_name"),
                                    values =list(regulatory.features,chromosomes) ,      # We are excluding regulatory features not mapped to chr1:22, X, Y
                                    mart =regulatory.hg38)

    gene=hg38.regulatory[,1]
    chrom=hg38.regulatory[,2]
    loc.start=hg38.regulatory[,3]
    loc.end=hg38.regulatory[,4]
    description=hg38.regulatory[,5]

    regulation.data=cbind.data.frame(gene=gene,
                                     chrom=chrom,
                                     loc.start=loc.start,
                                     loc.end=loc.end,
                                     description=description)

    namevector <- c("gene.name", "biotype", "chrom.strand", "chrom.band")
    regulation.data[ , namevector] <- NA
    regulation.data.final=rbind.data.frame(regulation.data, gene.data)

    # To retrieve data for experimentally validated regulatory features mapped to GRCh38 (FANTOM project)
    message(paste0("Retrieving experimentally validated regulatory features (FANTOM project): ",date()))

    reg.hg38.val = biomaRt::useEnsembl(biomart="regulation",
                                       dataset="hsapiens_external_feature", # regulatory features includes (promoters, promoter flanking regions, enhancers, CTCF binding sites, TF binding sites and open chromatin regions)
                                       version = '110')                     # specify version for stable results

    # specify data elements to retreive for each regulatory feature
    hg38.reg.val= biomaRt::getBM(attributes=c("display_label", "chromosome_name",
                                              "chromosome_start","chromosome_end",
                                              "feature_type_description", "feature_type_class"),
                                 filters =c("external_feature_set_name", "chromosome_name"),
                                 values =list("FANTOM predictions",chromosomes) ,      # We are excluding regulatory features not mapped to chr1:22, X, Y
                                 mart =reg.hg38.val)

    gene=hg38.reg.val[,1]
    chrom=hg38.reg.val[,2]
    loc.start=hg38.reg.val[,3]
    loc.end=hg38.reg.val[,4]
    description=hg38.reg.val[,5]

    regulation.data.val=cbind.data.frame(gene=gene,
                                         chrom=chrom,
                                         loc.start=loc.start,
                                         loc.end=loc.end,
                                         description=description)

    namevector <- c("gene.name", "biotype", "chrom.strand", "chrom.band")
    regulation.data.val[ , namevector] <- NA

    regulation.data.val.final=rbind.data.frame(regulation.data.val, gene.data)

    # return gene and regulatory features annotation data in one list
    res=list(gene.annotation=gene.data,                          # gene.data represent annotated genes
             reg.annotation.predicted=regulation.data.final,     # predicted regulatory features (ensembl regulatory build)
             reg.annotation.validated=regulation.data.val.final) # experimentally validated regulatory features (FANTOM project)

    return(res)
  }

  # retrieve gene annotation data for human_GRCh37 (hg19) genome assembly from ensembl biomaRt version 75
  if (genome.assembly=="Human_GRCh37")
  {
    message(paste0("Retrieving gene annotation data from Ensembl BioMart: ",date()))
    ensemblGRCh37 <- biomaRt::useEnsembl(biomart = 'ensembl',
                                         dataset = 'hsapiens_gene_ensembl',
                                         version = '75')

    chromosomes = c(1:22, "X", "Y")    # specify chromosomes of interst

    # specify data elements to retreive for each gene
    hg19_gene_annotation <- biomaRt::getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                                        'start_position','end_position',
                                                        'description','external_gene_id',
                                                        'gene_biotype', 'strand', 'band'),
                                           filters = 'chromosome_name',  values = chromosomes,
                                           mart = ensemblGRCh37)

    gene=hg19_gene_annotation[,1]
    chrom=hg19_gene_annotation[,2]
    loc.start=hg19_gene_annotation[,3]
    loc.end=hg19_gene_annotation[,4]
    description=hg19_gene_annotation[,5]
    gene.name=hg19_gene_annotation[,6]
    biotype=hg19_gene_annotation[,7]
    chrom.strand=hg19_gene_annotation[,8]
    chrom.band=hg19_gene_annotation[,9]

    gene.data=cbind.data.frame(gene=gene,
                               chrom=chrom,
                               loc.start=loc.start,
                               loc.end=loc.end,
                               description=description,
                               gene.name=gene.name,
                               biotype=biotype,
                               chrom.strand=chrom.strand,
                               chrom.band=chrom.band)

    # To retrieve data for regulatory features mapped to GRCh37
    message(paste0("Retrieving computationally predicted regulatory features from Ensembl regulatory build: ",date()))

    regulatory.hg19 = biomaRt::useEnsembl(biomart="regulation",
                                          dataset="hsapiens_regulatory_feature",
                                          version = 'GRCh37')

    # specify data elements to retreive for each regulatory feature
    regulatory.features=c("Promoter", "Promoter Flanking Region", "Enhancer",
                          "CTCF Binding Site", "TF binding site", "Open chromatin")

    hg19.regulatory= biomaRt::getBM(attributes=c("regulatory_stable_id", "chromosome_name",
                                                 "chromosome_start","chromosome_end",
                                                 "feature_type_description"),
                                    filters =c("regulatory_feature_type_name", "chromosome_name"),
                                    values =list(regulatory.features,chromosomes) ,      # We are excluding regulatory features not mapped to chr1:22, X, Y
                                    mart =regulatory.hg19)

    gene=hg19.regulatory[,1]
    chrom=hg19.regulatory[,2]
    loc.start=hg19.regulatory[,3]
    loc.end=hg19.regulatory[,4]
    description=hg19.regulatory[,5]

    regulation.data=cbind.data.frame(gene=gene,
                                     chrom=chrom,
                                     loc.start=loc.start,
                                     loc.end=loc.end,
                                     description=description)

    namevector <- c("gene.name", "biotype", "chrom.strand", "chrom.band")
    regulation.data[ , namevector] <- NA
    regulation.data.final=rbind.data.frame(regulation.data, gene.data)

    # To retrieve data for experimentally validated regulatory features mapped to GRCh37 (FANTOM project)
    message(paste0("Retrieving experimentally validated regulatory features (FANTOM project): ",date()))

    reg.hg19.val = biomaRt::useEnsembl(biomart="regulation",
                                       dataset="hsapiens_external_feature", # regulatory features includes (promoters, promoter flanking regions, enhancers, CTCF binding sites, TF binding sites and open chromatin regions)
                                       version = 'GRCh37')                  # specify version for stable results

    hg19.reg.val= biomaRt::getBM(attributes=c("display_label", "chromosome_name",
                                              "chromosome_start","chromosome_end",
                                              "feature_type_description", "feature_type_class"),
                                 filters =c("external_feature_set_name", "chromosome_name"),
                                 values =list("FANTOM predictions",chromosomes) ,      # We are excluding regulatory features not mapped to chr1:22, X, Y
                                 mart =reg.hg19.val)

    gene=hg19.reg.val[,1]
    chrom=hg19.reg.val[,2]
    loc.start=hg19.reg.val[,3]
    loc.end=hg19.reg.val[,4]
    description=hg19.reg.val[,5]

    regulation.data.val=cbind.data.frame(gene=gene,
                                         chrom=chrom,
                                         loc.start=loc.start,
                                         loc.end=loc.end,
                                         description=description)

    namevector <- c("gene.name", "biotype", "chrom.strand", "chrom.band")
    regulation.data.val[ , namevector] <- NA
    regulation.data.val.final=rbind.data.frame(regulation.data.val, gene.data)

    # return gene and regulatory features annotation data in one list
    res=list(gene.annotation=gene.data,                          # gene.data represent annotated genes
             reg.annotation.predicted=regulation.data.final,     # predicted regulatory features (ensembl regulatory build)
             reg.annotation.validated=regulation.data.val.final) # experimentally validated regulatory features (FANTOM5 project)

    return(res)
  }

  # retrieve gene annotation data for Mouse_HGCm39 genome assembly from ensembl biomaRt version 104
  if (genome.assembly=="Mouse_HGCm39")
  {
    message(paste0("Retrieving gene annotation data from Ensembl BioMart: ",date()))
    ensembl.HGCm39 = biomaRt::useEnsembl(biomart="genes",
                                         dataset="mmusculus_gene_ensembl",
                                         version = "104")

    chromosomes = c(1:19, "X", "Y")
    HGCm39.gene_annotation <- biomaRt::getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                                          'start_position','end_position',
                                                          'description','external_gene_name',
                                                          'gene_biotype', 'strand', 'band'),
                                             filters = 'chromosome_name',  values = chromosomes,
                                             mart = ensembl.HGCm39)

    gene=HGCm39.gene_annotation[,1]
    chrom=HGCm39.gene_annotation[,2]
    loc.start=HGCm39.gene_annotation[,3]
    loc.end=HGCm39.gene_annotation[,4]
    description=HGCm39.gene_annotation[,5]
    gene.name=HGCm39.gene_annotation[,6]
    biotype=HGCm39.gene_annotation[,7]
    chrom.strand=HGCm39.gene_annotation[,8]
    chrom.band=HGCm39.gene_annotation[,9]

    gene.data=cbind.data.frame(gene=gene,
                               chrom=chrom,
                               loc.start=loc.start,
                               loc.end=loc.end,
                               description=description,
                               gene.name=gene.name,
                               biotype=biotype,
                               chrom.strand=chrom.strand,
                               chrom.band=chrom.band)

    # To retrieve data for mouse_regulatory features mapped to mm39
    message(paste0("Retrieving computationally predicted regulatory features from Ensembl regulatory build: ",date()))

    regulatory.mm39 = biomaRt::useEnsembl(biomart="regulation",
                                          dataset="mmusculus_regulatory_feature",
                                          version = '104')

    regulatory.features=c("Promoter", "Promoter Flanking Region", "Enhancer",
                          "CTCF Binding Site", "TF binding site", "Open chromatin")
    mm39.regulatory= biomaRt::getBM(attributes=c("regulatory_stable_id", "chromosome_name",
                                                 "chromosome_start","chromosome_end",
                                                 "feature_type_description"),
                                    filters =c("regulatory_feature_type_name", "chromosome_name"),
                                    values =list(regulatory.features,chromosomes) ,      # We are excluding regulatory features not mapped to chr1:22, X, Y
                                    mart =regulatory.mm39)

    gene=mm39.regulatory[,1]
    chrom=mm39.regulatory[,2]
    loc.start=mm39.regulatory[,3]
    loc.end=mm39.regulatory[,4]
    description=mm39.regulatory[,5]

    regulation.data=cbind.data.frame(gene=gene,
                                     chrom=chrom,
                                     loc.start=loc.start,
                                     loc.end=loc.end,
                                     description=description)

    namevector <- c("gene.name", "biotype", "chrom.strand", "chrom.band")
    regulation.data[ , namevector] <- NA
    regulation.data.final=rbind.data.frame(regulation.data, gene.data)

    # return gene and regulatory features annotation data in one list
    res=list(gene.annotation=gene.data,
             reg.annotation=regulation.data.final)
    return(res)
  }

  # retrieve gene annotation data for Mouse_HGCm38 from ensembl biomaRt version 102
  if (genome.assembly=="Mouse_HGCm38")
  {
    message(paste0("Retrieving gene annotation data from Ensembl BioMart: ",date()))
    ensemblGRCm38 <- biomaRt::useEnsembl(biomart = 'genes',
                                         dataset = 'mmusculus_gene_ensembl',
                                         version = '102')

    chromosomes = c(1:19, "X", "Y") # chromosomes of interest
    Mouse_mm10.gene_annotation <- biomaRt::getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                                              'start_position','end_position',
                                                              'description','external_gene_name',
                                                              'gene_biotype', 'strand', 'band'),
                                                 filters = 'chromosome_name',  values = chromosomes,
                                                 mart = ensemblGRCm38)

    gene=Mouse_mm10.gene_annotation[,1]
    chrom=Mouse_mm10.gene_annotation[,2]
    loc.start=Mouse_mm10.gene_annotation[,3]
    loc.end=Mouse_mm10.gene_annotation[,4]
    description=Mouse_mm10.gene_annotation[,5]
    gene.name=Mouse_mm10.gene_annotation[,6]
    biotype=Mouse_mm10.gene_annotation[,7]
    chrom.strand=Mouse_mm10.gene_annotation[,8]
    chrom.band=Mouse_mm10.gene_annotation[,9]

    gene.data=cbind.data.frame(gene=gene,
                               chrom=chrom,
                               loc.start=loc.start,
                               loc.end=loc.end,
                               description=description,
                               gene.name=gene.name,
                               biotype=biotype,
                               chrom.strand=chrom.strand,
                               chrom.band=chrom.band)

    # To retrieve data for mouse_regulatory features mapped to mm10
    message(paste0("Retrieving computationally predicted regulatory features from Ensembl regulatory build: ",date()))

    regulatory.mm10 = biomaRt::useEnsembl(biomart="regulation",
                                          dataset="mmusculus_regulatory_feature",
                                          version = '102')
    regulatory.features=c("Promoter", "Promoter Flanking Region", "Enhancer",
                          "CTCF Binding Site", "TF binding site", "Open chromatin")

    mm10.regulatory= biomaRt::getBM(attributes=c("regulatory_stable_id", "chromosome_name",
                                                 "chromosome_start","chromosome_end",
                                                 "feature_type_description"),
                                    filters =c("regulatory_feature_type_name", "chromosome_name"),
                                    values =list(regulatory.features,chromosomes),      # We are excluding regulatory features not mapped to chr1:22, X, Y
                                    mart =regulatory.mm10)

    gene=mm10.regulatory[,1]
    chrom=mm10.regulatory[,2]
    loc.start=mm10.regulatory[,3]
    loc.end=mm10.regulatory[,4]
    description=mm10.regulatory[,5]

    regulation.data=cbind.data.frame(gene=gene,
                                     chrom=chrom,
                                     loc.start=loc.start,
                                     loc.end=loc.end,
                                     description=description)

    namevector <- c("gene.name", "biotype", "chrom.strand", "chrom.band")
    regulation.data[ , namevector] <- NA
    regulation.data.final=rbind.data.frame(regulation.data, gene.data)

    # return gene and regulatory features annotation data in one list
    res=list(gene.annotation=gene.data,
             reg.annotation=regulation.data.final)

    return(res)
  }
}
