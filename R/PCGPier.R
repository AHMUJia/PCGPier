#' @name PCGPier
#' @title Calculate the gene-pair based prognostic signature PCGPier for prostate cancer samples
#' @description This prognosis prediction study established a gene-pairs based RFS prediction classifier for prostate cancer patients, which is universal applicable in multiple cohorts with different gene expression detection methods.
#' @param expr A numerical expression matrix or data frame with row for gene symbol name and column for sample ID. Note: In principle, the expression profile does not need any normalization, but since the amount of gene expression is affected by the gene length, the original count or the normalized count data may not be suitable for this analysis. It is recommended to provide FPKM or TPM value and log transformation is not necessary.
#' @param res.path A string value to indicate the output path for storing the Estimated IRGPI.txt file.
#'
#' @return A data frame stored sample name and estimated IRGPI.
#' @export
#' @author Jialin Meng, Xiaofan Lu, Tong Jin
#' @examples
#' library(PCGPier)
#' load(system.file("extdata", "demo.RData", package = "PCGPier", mustWork = TRUE)) # load example data
#' riskscore  <- PCGPier(expr = demo)

PCGPier <- function(expr = NULL,res.path = getwd())
{

  # genes in Riskscore1
  PCGP1 <- c("CDCA3",
             "HOPX",
             "CLP1",
             "DUSP12",
             "CELSR3",
             "ASTE1",
             "MGAT3",
             "FSCN2",
             "TPX2",
             "CD86",
             "AURKA",
             "RAD54B",
             "MUM1",
             "AURKA",
             "CAMK4",
             "CTAGE5",
             "FGF1",
             "DNAJC9",
             "CACNA2D3",
             "CAMK1G",
             "MOGAT2",
             "FAM49B",
             "ALPK1",
             "ADCY1",
             "CCNJ",
             "HPSE2")

  # genes in Riskscore2
  PCGP2 <- c("VPS53",
             "HPSE2",
             "GGH",
             "TTC17",
             "MORC2",
             "PIGB",
             "MLYCD",
             "MAPKBP1",
             "UBE3B",
             "KPTN",
             "SMG6",
             "RNF32",
             "NUP214",
             "USPL1",
             "KIAA1644",
             "PDE2A",
             "KCNIP1",
             "MMS19",
             "CCNE2",
             "CDCA3",
             "TREM1",
             "TPX2",
             "HELLS",
             "BRCA1",
             "XYLT1",
             "ZFP37")

  # coefficient for each PCGP
  coeff <- c(-0.792135779,
             -0.403762448,
             -0.375114609,
             -0.363169154,
             -0.360724876,
             -0.337901672,
             -0.322580743,
             -0.253583755,
             -0.143328214,
             -0.123739185,
             -0.116529692,
             -0.114518586,
             -0.097742537,
             -0.069928613,
             -0.036942305,
             -0.024182837,
             -0.01575333,
             -0.004845988,
             0.02990844,
             0.134100101,
             0.158840217,
             0.199349467,
             0.225710522,
             0.276381397,
             0.291181288,
             0.386132311)

  # initial check
  if(max(expr) >= 25){
    message("--please make sure a properly normalized expression data has been provided (e.g., FPKM or TPM); count data is not suitable because it does not consider gene length.")
  }

  if(!all(is.element(unique(c(PCGP1,PCGP2)), rownames(expr)))) {
    missgene <- setdiff(unique(c(PCGP1,PCGP2)), rownames(expr))
    stop(paste0(length(missgene)," genes cannot be mapped in your data, please check for the following missing feature(s):\n",
                paste(missgene, collapse = "\n")))
  }


  # extract expression
  PCGP1.expr <- expr[PCGP1,]
  PCGP2.expr <- expr[PCGP2,]

  # calculate difference and convert to binary matrix
  PCGP.diff <- PCGP1.expr - PCGP2.expr
  PCGP.binary <- ifelse(PCGP.diff < 0, 1, 0)

  # estimate Riskscore
  Riskscore <- apply(t(PCGP.binary),1,function(x) {x %*% coeff})

  # output
  outTab <- data.frame(SampleID = names(Riskscore),
                       Riskscore = as.numeric(Riskscore),
                       Group = ifelse(Riskscore>=-2.820,"HRisk","LRisk"),
                       stringsAsFactors = F)
  write.table(outTab,
              file = file.path(res.path,"Estimated Riskscore.txt"),
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)

  return(outTab)
}

