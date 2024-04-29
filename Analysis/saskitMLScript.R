library(igraph)
library(caret)
library(randomForestSRC)
library(survival)
library(limma)
library(edgeR)
library(GSVA)
data("c2BroadSets")
library(survcomp)
library(preprocessCore)

# Classes -----------------------------------------------------------------

# SASKit container to hold all data and outputs.

saskit <- setClass("saskit", slots = list(demographic = "data.frame",
                                          clinical = "data.frame",
                                          transcriptomics = "matrix",
                                          proteomics = "matrix",
                                          transcriptomic_features = "list",
                                          proteomic_features = "list",
                                          saskit_model = "list",
                                          settings = "list"
)
)

# Function to allow the saskit object to be subsetted by subjects using [1:n, ].

saskit_subset <- setMethod("[", signature(x = "saskit", i = "ANY"),
                           function(x, i, ..., drop) {
                             if (nrow(x@demographic) > 0) {
                               new_demographic <- x@demographic[i, , drop = FALSE]
                             } else {
                               new_demographic <- x@demographic
                             }
                             if (nrow(x@clinical) > 0) {
                               new_clinical <- x@clinical[i, , drop = FALSE]
                             } else {
                               new_clinical <- x@clinical
                             }
                             if (length(x@transcriptomics) > 0) {
                               new_transcriptomics <- x@transcriptomics[, i, drop = FALSE]
                             } else {
                               new_transcriptomics <- x@transcriptomics
                             }
                             if (length(x@proteomics) > 0) {
                               new_proteomics <- x@proteomics[, i, drop = FALSE]
                             } else {
                               new_proteomics <- x@proteomics
                             }
                             new_object <- new("saskit", demographic = new_demographic, clinical = new_clinical, transcriptomics = new_transcriptomics, proteomics = new_proteomics, transcriptomic_features = x@transcriptomic_features, proteomic_features = x@proteomic_features, saskit_model = x@saskit_model, settings = x@settings)
                             return(new_object)
                           }
)

# Calculate GSVA ----------------------------------------------------------

# Function to calculate the GSVA pathway activation scores. Takes a SASkit
# object x and whether "proteomics" or "transcriptomics" are to be used to
# calculate the features.

calculateGSVA <- function(x, rankFC = F, gsvaN = 5, pathways = NULL, RSF = F, omics) {
  subjects <- x@demographic$SubjectID
  conditions <- x@demographic$condition

  # Get KEGG pathways.
  canonicalC2BroadSets <- readRDS("C2BroadSetsEnsembl.R")

  if (omics == "transcriptomics") {
    normalized_transcriptomics <- cpm(x@transcriptomics)

    gsvaPar <- gsvaParam(normalized_transcriptomics, canonicalC2BroadSets, minSize = 5)
    gsvaResult <- gsva(gsvaPar, verbose=FALSE, )
  }


  if (omics == "proteomics") {
    gsvaPar <- gsvaParam(x@proteomics, canonicalC2BroadSets, minSize = 5)
    gsvaResult <- gsva(gsvaPar, verbose=FALSE)
  }



  # Perform GSVA differential analysis.

  if (length(table(x@demographic$sex)) > 1) {
    sample_info <- data.frame(subject = x@demographic$SubjectID,
                              condition = x@demographic$condition,
                              age = x@demographic$age,
                              sex = x@demographic$sex
    )
    gsvaDesign <- model.matrix(~ 0 + condition + age + sex, data = x@demographic)
  } else {
    sample_info <- data.frame(subject = x@demographic$SubjectID,
                              condition = x@demographic$condition,
                              age = x@demographic$age
    )
    gsvaDesign <- model.matrix(~ 0 + condition + age, data = x@demographic)
  }


  fit <- lmFit(gsvaResult, gsvaDesign)
  fit2 <- eBayes(fit)

  gsvaDiff <- topTable(fit2, coef = 1, number = Inf)

  if (rankFC == F) {
    gsvaDiff <- gsvaDiff[order(gsvaDiff[ , 4], decreasing = F), ]
  } else if (rankFC == T) {
    gsvaDiff <- gsvaDiff[order(abs(gsvaDiff[ , 1]), decreasing = T), ]
  }

  # Identify the N most important pathways and return their scores for all subjects.
  if (is.null(pathways)) {
    gsvaTopN <- rownames(gsvaDiff)[1:gsvaN]
  } else {
    gsvaTopN <- pathways
  }

  gsvaFeatures <- gsvaResult[which(rownames(gsvaResult) %in% gsvaTopN), ]
  if (omics == "transcriptomics") {
    rownames(gsvaFeatures) <- paste("transcriptomics", rownames(gsvaFeatures), sep = "_")
  }
  if (omics == "proteomics") {
    rownames(gsvaFeatures) <- paste("proteomics", rownames(gsvaFeatures), sep = "_")
  }

  gsvaAll <- list(topGSVA = gsvaDiff[1:gsvaN, ], pathwaysUsed = gsvaTopN,  gsvaFeatures = gsvaFeatures,
                  gsvaAll = gsvaResult)
}

# Caluclate ExprEssence ---------------------------------------------------

# Function to calculate the ExprEssence average linkscores. Takes a SASKit object
# x and whether "proteomics" or "transcriptomics" are to be used to
# calculate the features.

calculateExprEssence <- function(x, net, labs = c("case", "control"), cutoff = 50, exprN = 5, omics, subs = NULL) {
  if (omics == "transcriptomics") {
    expressionValues <- x@transcriptomics
    expressionValues <- normalize.quantiles(cpm(expressionValues, log = TRUE))
    rownames(expressionValues) <- rownames(x@transcriptomics)
    colnames(expressionValues) <- colnames(x@transcriptomics)
  }

  if (omics == "proteomics") {
    expressionValues <- x@proteomics
    expressionValues <- normalize.quantiles(expressionValues)
    rownames(expressionValues) <- rownames(x@proteomics)
    colnames(expressionValues) <- colnames(x@proteomics)
  }

  # Trim the network of any genes without expression values.

  # Function to check if nodes are in the expressionValues row names.
  checkNode <- function(node) {
    ifelse(node %in% rownames(expressionValues), node, NA)
  }

  # Apply the checkNode function to each column of the network.
  net[, 1] <- sapply(net[, 1], checkNode)
  net[, 2] <- sapply(net[, 2], checkNode)

  # Convert to character and recreate the data frame.
  net <- data.frame(Gene1 = as.character(net[, 1]), Gene2 = as.character(net[, 2]), stringsAsFactors = FALSE)

  # Remove rows with NA values
  net <- na.omit(net)

  # Annotate the trimmed network with the required expression values.

  # First we need to take an average of expression values for each gene for each condition.

  # Function to compute average expression values for a given condition label.
  computeAverageExpression <- function(conditionLabel) {
    expression <- rowMeans(expressionValues[, x@demographic$condition == conditionLabel])
    aggregatedExpression <- aggregate(expression, by = list(Gene = names(expression)), FUN = mean, na.rm = TRUE)
    # Convert the aggregated data frame to a named numeric vector.
    setNames(as.numeric(aggregatedExpression$x), aggregatedExpression$Gene)
  }

  # Calculate expression averages for case and control.
  caseExpr <- computeAverageExpression(labs[1])
  controlExpr <- computeAverageExpression(labs[2])

  # Calculate expression changes for Gene1 and Gene2 in the network.
  net$ExprChangeGene1 <- caseExpr[net[, 1]] - controlExpr[net[, 1]]
  net$ExprChangeGene2 <- caseExpr[net[, 2]] - controlExpr[net[, 2]]

  # Sum the expression changes for Gene1 and Gene2 to compute the Linkscore.
  net$Linkscore <- net$ExprChangeGene1 + net$ExprChangeGene2

  # Create an igraph object for further manipulation and visualization.
  ig <- igraph::simplify(graph_from_data_frame(net, directed = FALSE), edge.attr.comb = "first")

  # Obtain a graph containing only the top-N (default 50) edges by linkscore.
  linkAttr <- E(ig)$Linkscore

  # If the number of edges exceeds the cutoff, only keep the top-N edges.
  if (length(linkAttr) > cutoff) {
    # Order edges by Linkscore and find indices of edges to delete
    deleteIx <- order(linkAttr, decreasing = TRUE)[(cutoff + 1):length(linkAttr)]
    ig <- delete_edges(ig, E(ig)[deleteIx])
  }

  # Remove isolated nodes after edge deletion
  ig_exprEssence <- igraph::delete.vertices(ig, which(igraph::degree(ig) == 0))

  # Calculate the average linkscores of the top-N (default 5) components.
  subnets <- components(ig_exprEssence)

  # Define a function to extract subnetworks and calculate their mean linkscores
  getSubnetwork <- function(member, ig_graph) {
    listOfNodes <- names(subnets$membership)[which(subnets$membership == member)]
    induced.subgraph(ig_graph, listOfNodes)
  }

  ig_exprEssence <- lapply(1:subnets$no, getSubnetwork, ig_exprEssence)

  getAverage <- function(ig_graph) {
    linkscores <- get.edge.attribute(ig_graph, "Linkscore")
    mean(linkscores)
  }

  linkscores <- sapply(ig_exprEssence, getAverage)


  # Determine indices to keep based on linkscores, limiting to exprN or the available number of components

  if (subnets$no == 0) {
    print("Could not calculate ExprEssence")
  } else {
    if (subnets$no < exprN) {
      if (!is.null(subs)) {
        exprN <- length(subs)
      } else {
        print(paste("Only ", subnets$no, " out of ", exprN, "ExprEssence features generated."))
        exprN <- subnets$no
      }
    }
  }
  keepIx <- sort(linkscores, index.return = TRUE, decreasing = TRUE)$ix[1:exprN]

  # Subset ig_exprEssence using the indices from keepIx
  ig_exprEssence <- ig_exprEssence[keepIx]

  # Get the average linkscores for each of the top-N networks for each subject, comparing each subject to the average control value.
  getLinkscores <- function(x) {
    controlExprOrdered <- controlExpr[order(names(controlExpr))]
    xOrdered <- x[order(rownames(x)), ]
    ExprChangeMatrix <- sweep(xOrdered, 1, controlExprOrdered)

    ExprChangeGene1 <- ExprChangeMatrix[match(net$Gene1, rownames(ExprChangeMatrix)), ]
    ExprChangeGene2 <- ExprChangeMatrix[match(net$Gene2, rownames(ExprChangeMatrix)), ]

    Linkscore <- ExprChangeGene1 + ExprChangeGene2
    rownames(Linkscore) <- paste(net$Gene1, net$Gene2, sep = ".")

    return(Linkscore)
  }

  class(expressionValues) <- "numeric"
  linkscoreResults <- getLinkscores(expressionValues)

  identifyEdges <- function(a) {
    paste(igraph::as_edgelist(ig_exprEssence[[a]])[ , 1], paste(igraph::as_edgelist(ig_exprEssence[[a]])[ , 2]), sep = "--")
  }

  if (is.null(subs)) {
    importantNets <- sapply(1:length(ig_exprEssence), identifyEdges)
  } else {
    importantNets <- subs
  }

  nodePairs <- paste(net[ , 1], net[ , 2], sep = "--")

  annotateLinkscores <- function(target, lst) {
    for (i in seq_along(lst)) {
      if (target %in% lst[[i]]) {
        if (omics == "transcriptomics") {
          if (exprN == 1) {
            return(paste("transcriptomics_exprEssence", "1", sep = "_"))
          } else {
            return(paste("transcriptomics_exprEssence", i, sep = "_"))
          }
        }
        if (omics == "proteomics") {
          return(paste("proteomics_exprEssence", i, sep = "_"))
        }
      }
    }
    return(NULL)  # Return NULL if the string is not found in any vector
  }

  rownames(linkscoreResults) <- nodePairs
  linkscoreResults <- as.data.frame(linkscoreResults)

  linkscoreResults$Subnetworks <- as.character(sapply(nodePairs, function(x) annotateLinkscores(x, importantNets) ))

  # Calculate the average linkscore for each ExprEssence subnetwork and return it for each patient.
  getAverageLinks <- function(x) {
    aggregate(x ~ linkscoreResults$Subnetworks, linkscoreResults[ , 1:(ncol(linkscoreResults)-1)], mean)[ , 2]
  }
  averageLinks <- sapply(linkscoreResults[ , 1:(ncol(linkscoreResults)-1)], getAverageLinks)

  if(exprN == 1) {
    if (length(subs) == 1) {
      names(averageLinks) <- c(paste(omics, "_exprEssence_1"))
    } else {
      rownames(averageLinks) <- c("NULL", paste(omics, "_exprEssence_1"))
    }
  } else {
    rownames(averageLinks) <- aggregate(linkscoreResults[ , 1] ~ linkscoreResults$Subnetworks, linkscoreResults[ , 1:(ncol(linkscoreResults)-1)], mean)[ , 1]

  }
  if (length(subs) != 1) {
    exprEssenceResult <- averageLinks[!rownames(averageLinks) == "NULL", ]
  }
  if (length(subs) == 1) {
    exprEssenceResult <- averageLinks[2 , ]
  }

  if (exprN == 1) {
    importantNets <- list(expressence = importantNets)
    if (omics == "transcriptomics") {
      exprEssenceResult <- t(data.frame(transcriptomics_exprEssence_1 = exprEssenceResult))
    }
    if (omics == "proteomics") {
      exprEssenceResult <- t(data.frame(proteomics_exprEssence_1 = exprEssenceResult))
    }
    #names(exprEssenceResult) <- "transcriptomics_exprEssence_1"
  }
  paste(exprEssenceResult)
  exprEssenceAll <- list(topExprEssence = importantNets, exprEssenceFeatures = exprEssenceResult)
}

# Caluclate Omics ---------------------------------------------------------

# Function to determine all required omics features. Takes a SASKit object x and
# whether "proteomics" or "transcriptomics" are to be used to
# calculate the features.

calculate_omics_features <- function(x, platform) { # Set platform to "transcriptomics" or "proteomics"
  # Determine what platform is being used.
  if (platform == "transcriptomics") {
    gsva_options <- x@settings$transcriptomic_options$gsva_options
    expressence_options <- x@settings$transcriptomic_options$expressence_options
  }
  if (platform == "proteomics") {
    gsva_options <- x@settings$proteomic_options$gsva_options
    expressence_options <- x@settings$proteomic_options$expressence_options
  }

  # Calculate GSVA features, if requested.
  if (gsva_options$run_gsva == TRUE) {
    gsva_results <- calculateGSVA(x, gsvaN = gsva_options$gsvaN, rankFC = gsva_options$rankFC, pathways = gsva_options$pathways, omics = platform)
  }
  # Calculate ExprEssence features, if requested.
  if (expressence_options$run_expressence == TRUE) {
    expressence_results <- calculateExprEssence(x, labs = expressence_options$labs, net = net, exprN = expressence_options$exprN, cutoff = expressence_options$cutoff, subs = expressence_options$subs, omics = platform)
  }
  omics_results <- list()
  if (gsva_options$run_gsva == TRUE) {
    omics_results <- c(omics_results, list(gsva_results = gsva_results))
  }
  if (expressence_options$run_expressence == TRUE) {
    omics_results <- c(omics_results, list(expressence_results = expressence_results))
  }
  return(omics_results)
}

# Perform RSF -------------------------------------------------------------

# Function to run a random survival forest analysis using either the SASKit features,
# or the expression values of all differentially expressed genes.

run_rsf <- function(x, age = TRUE, sex = TRUE) {
  saskit_options <- x@settings$saskit_options
  rsf_options <- x@settings$rsf_options

  transcriptomics_gsva_options <- x@settings$transcriptomic_options$gsva_options
  transcriptomics_expressence_options <- x@settings$transcriptomic_options$expressence_options

  proteomics_gsva_options <- x@settings$proteomic_options$gsva_options
  proteomics_expressence_options <- x@settings$proteomic_options$expressence_options

  featureTable <- data.frame(
    Status = as.numeric(x@demographic$Status),
    Time = as.numeric(x@demographic$Time)
  )

  if (rsf_options$age == TRUE) {
    Age <- x@demographic$age
    featureTable <- cbind(featureTable, Age)
  }

  if (rsf_options$sex == TRUE) {
    Sex <- as.factor(x@demographic$sex)
    featureTable <- cbind(featureTable, Sex)
  }

  if(saskit_options$transcriptomics == TRUE) {
    if (transcriptomics_gsva_options$run_gsva == TRUE) {
      featureTable <- cbind(featureTable, t(x@transcriptomic_features$gsva$gsvaFeatures))
    }
    if (transcriptomics_expressence_options$run_expressence == TRUE) {
      featureTable <- cbind(featureTable, t(x@transcriptomic_features$expressence$exprEssenceFeatures))
    }
  }
  if (saskit_options$proteomics == TRUE) {
    if (proteomics_gsva_options$run_gsva == TRUE) {
      featureTable <- cbind(featureTable, t(x@proteomic_features$gsva$gsvaFeatures))
    }
    if (proteomics_expressence_options$run_expressence == TRUE) {
      featureTable <- cbind(featureTable, t(x@proteomic_features$expressence$exprEssenceFeatures))
    }
  }
  if (saskit_options$clinical == TRUE) {
    featureTable <- cbind(featureTable, x@clinical)
  }

  featureTable <- featureTable[x@demographic$condition == x@settings$transcriptomic_options$expressence_options$labs[1], ]
  featureTable <- featureTable[complete.cases(featureTable), ]

  saskitModel <- rfsrc(Surv(Time, Status) ~ ., data = featureTable,
                       importance = TRUE,
                       refit = FALSE,
                       block.size = 1
  )
}

# Validate RSF ------------------------------------------------------------

# Function to perform 10-fold cross-validation on the rsf model produced by
# run_rsf().

validate_rsf <- function(x, rsf_model, features = TRUE) {
  saskit_options <- x@settings$saskit_options
  rsf_options <- x@settings$rsf_options

  transcriptomics_dega_options <- x@settings$transcriptomic_options$dega_options
  transcriptomics_gsva_options <- x@settings$transcriptomic_options$gsva_options
  transcriptomics_expressence_options <- x@settings$transcriptomic_options$expressence_options

  proteomics_dega_options <- x@settings$proteomic_options$dega_options
  proteomics_gsva_options <- x@settings$proteomic_options$gsva_options
  proteomics_expressence_options <- x@settings$proteomic_options$expressence_options

  featureTable <- data.frame(
    Status = as.numeric(x@demographic$Status),
    Time = as.numeric(x@demographic$Time)
  )

  if (rsf_options$age == TRUE) {
    Age <- x@demographic$age
    featureTable <- cbind(featureTable, Age)
  }

  if (rsf_options$sex == TRUE) {
    Sex <- as.factor(x@demographic$sex)
    featureTable <- cbind(featureTable, Sex)
  }

  if (features == TRUE) {
    if(saskit_options$transcriptomics == TRUE) {
      if (transcriptomics_gsva_options$run_gsva == TRUE) {
        featureTable <- cbind(featureTable, t(x@transcriptomic_features$gsva$gsvaFeatures))
      }
      if (transcriptomics_expressence_options$run_expressence == TRUE) {
        featureTable <- cbind(featureTable, t(x@transcriptomic_features$expressence$exprEssenceFeatures))
      }
    }
    if (saskit_options$proteomics == TRUE) {
      if (proteomics_gsva_options$run_gsva == TRUE) {
        featureTable <- cbind(featureTable, t(x@proteomic_features$gsva$gsvaFeatures))
      }
      if (proteomics_expressence_options$run_expressence == TRUE) {
        featureTable <- cbind(featureTable, t(x@proteomic_features$expressence$exprEssenceFeatures))
      }
      if (proteomics_dega_options$run_dega == TRUE) {
        featureTable <- cbind(featureTable, t(x@proteomics[rownames(x@proteomics) %in% rownames(x@proteomic_features$dega_results), ]))
      }
    }
    if (saskit_options$clinical == TRUE) {
      featureTable <- cbind(featureTable, x@clinical)
    }

    featureTable <- featureTable[x@demographic$condition == x@settings$transcriptomic_options$expressence_options$labs[1], ]
    featureTable <- featureTable[complete.cases(featureTable), ]


    predictions <- predict.rfsrc(rsf_model, featureTable)

    # Calculate c-index
    cindex <- get.cindex(time = featureTable$Time, censoring = featureTable$Status, predicted = predictions$predicted)

    ibs_featureTable <- data.frame(time = featureTable$Time, event = featureTable$Status, score = predictions$predicted)

    ibs <- sbrier.score2proba(ibs_featureTable, ibs_featureTable, method = c("cox"))

    prediction_accuracy <- list(cindex = cindex, ibs = ibs$bsc.integrated)
  }
}

# Run pipeline ------------------------------------------------------------

# Function to run the overall pipeline according to the specifications of the
# dataset.

runSaskit <- function(x, transcriptomics = TRUE,
                      proteomics = TRUE,
                      clinical = TRUE,
                      survival = TRUE, # If set to FALSE, a classification will be performed.
                      comparisonLabels = c("case", "control"),
                      gsva = 5, # These set the number of features for both gsva and expressence to be generated.
                      expressence = 5,
                      k = 10 # Fold crossvalidation
) {

  case_subjects <- x@demographic$SubjectID[x@demographic$condition == comparisonLabels[1]]
  control_subjects <- x@demographic$SubjectID[x@demographic$condition == comparisonLabels[2]]
  folds <- createFolds(as.factor(x@demographic$Status[x@demographic$condition == comparisonLabels[1]]), k = k)

  prediction_accuracy <- list()

  for (fold in 1:length(folds)) {

    fold_subjects <- c(case_subjects[-folds[[fold]]], control_subjects)
    test_subjects <- c(case_subjects[folds[[fold]]], control_subjects)

    saskit_fold <- x[which(x@demographic$SubjectID %in% fold_subjects), ]

    saskit_test <- x[which(x@demographic$SubjectID %in% test_subjects), ]


    # Transcriptomics ----------------------------------------------------------


    if (isTRUE(x@settings$saskit_options$transcriptomics)) {
      platform <- "transcriptomics"
      # Calculate the transcriptomic features for the current cross-validation fold.
      print(paste("Calculating training transcriptomics features for fold", fold, sep = " "))
      transcriptomic_features_fold <- calculate_omics_features(saskit_fold, platform = platform)
      saskit_fold@transcriptomic_features <- transcriptomic_features_fold

      print(paste("Calculating testing transcriptomics features for fold", fold, sep = " "))
      saskit_test@settings$transcriptomic_options$gsva_options$pathways <- saskit_fold@transcriptomic_features$gsva_results$pathwaysUsed
      saskit_test@settings$transcriptomic_options$expressence_options$subs <- saskit_fold@transcriptomic_features$expressence_results$topExprEssence
      transcriptomic_features_test <- calculate_omics_features(saskit_test, platform = platform)
      saskit_test@transcriptomic_features <- transcriptomic_features_test
    }

    # Proteomics --------------------------------------------------------------


    if (isTRUE(x@settings$saskit_options$proteomics)) {
      # Calculate the proteomic features for the current cross-validation fold.
      print(paste("Calculating training proteomics features for fold", fold, sep = " "))
      proteomic_features_fold <- calculate_omics_features(saskit_fold, platform = "proteomics")
      saskit_fold@proteomic_features <- proteomic_features_fold

      print(paste("Calculating testing proteomic features for fold", fold, sep = " "))
      saskit_test@settings$proteomic_options$gsva_options$pathways <- saskit_fold@proteomic_features$gsva_results$pathwaysUsed
      saskit_test@settings$proteomic_options$expressence_options$subs <- saskit_fold@proteomic_features$expressence_results$topExprEssence
      proteomic_features_test <- calculate_omics_features(saskit_test, platform = "proteomics")
      saskit_test@proteomic_features <- proteomic_features_test

    }

    # Random Survival Forest --------------------------------------------------

    print(paste("Running RSF analysis for fold ", fold, sep = ""))

    saskit_model_fold <- run_rsf(saskit_fold)

    print(paste("Running RSF analysis for test ", fold, sep = ""))

    prediction_accuracy_fold <- validate_rsf(saskit_test, rsf_model = saskit_model_fold, features = TRUE)
    prediction_accuracy <- c(prediction_accuracy, prediction_accuracy_fold)

  }

  # Transctipomics ----------------------------------------------------------


  if (isTRUE(x@settings$saskit_options$transcriptomics)) {
    print("Calculating main transcriptomic features.")
    # Calculate the transcriptomic features for the current cross-validation fold.
    transcriptomic_features <- calculate_omics_features(x, platform = "transcriptomics")
    x@transcriptomic_features <- transcriptomic_features
  }

  # Proteomics --------------------------------------------------------------

  if (isTRUE(x@settings$saskit_options$proteomics)) {
    print("Calculating main proteomic features.")
    # Calculate the proteomic features for the current cross-validation fold.
    proteomic_features <- calculate_omics_features(x, platform = "proteomics")
    x@proteomic_features <- proteomic_features
  }

  print("Running RSF analysis")

  saskitModel <- run_rsf(x, sex = x@settings$rsf_options$sex)

  result_list <- list(model = saskitModel, accuracy = prediction_accuracy, saskit = x)

  return(result_list)
}



# Analysis Begins ---------------------------------------------------------

# Load in Data ------------------------------------------------------------

# The STRING network for exprEssence.
net <- readRDS("STRING.R")

# The SASKit objects of each cohort.
prediabetesSASKit <- readRDS("data/diabetesSASKit.rds")
ibdSASKit <- readRDS("data/ibdSASKit.rds")
arthritisSASKit <- readRDS("data/arthritisSASKit.rds")

cancersSASKit <- lapply(list.files("data/", pattern = "\\.rds$", full.names = TRUE), readRDS)
names(cancersSASKit) <- gsub("SASKit.rds", "", list.files("data/", pattern = "\\.rds$", full.names = FALSE))

# Run analysis ------------------------------------------------------------

set.seed(123)
print("Running HMP-T2D")
prediabetes_results <- runSaskit(prediabetesSASKit, comparisonLabels = c("Prediabetic", "Control"))

print("Running HMP-IBD")
ibd_results <- runSaskit(ibdSASKit, comparisonLabels = c("IBD", "nonIBD"))

print("Running RAMAP")
arthritis_results <- runSaskit(arthritisSASKit, comparisonLabels = c("case", "control"))

# These will take a while to run!
#print("Running TCGA")
#cancer_results <- lapply(cancersSASKit, runSaskit, comparisonLabels = c("Cancer", "Normal"))
