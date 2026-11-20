#
##
#### EVOLUTIONARY MODEL FITTING 
#### USING SHAPE DATA + CONSENSUS AND POSTERIOR TREES
##
#

#### Script by Rafael A. Rivero Vega & Jacob S. Berv
#### Based off original and contributed to by Julien Clavel, with
#### further details in the manuscript methods and acknowledgements.




#### Load the required libraries ####

if(!'mvMORPH'%in%installed.packages()
  ) install.packages("mvMORPH"); library(mvMORPH)

if(!'geomorph'%in%installed.packages()
  ) install.packages("geomorph"); library(geomorph)

if(!'future.apply'%in%installed.packages()
  ) install.packages("future.apply"); library(future.apply)

if(!'HDInterval'%in%installed.packages()
  ) install.packages("HDInterval"); library(HDInterval)

if(!'vioplot'%in%installed.packages()
  ) install.packages("vioplot"); library(vioplot)

if(!'geiger'%in%installed.packages()
  ) install.packages("geiger"); library(geiger)

#rm(list = ls()) # if you need to clean up



#####

#### Specify path to working directory ####

#wd.lf <- readClipboard() # if you copied and want to paste the wd
#wd.lf <- setwd(choose.dir()) # if you want to popup select the wd 
#wd.lf <- "~/path/here" # if you want to manually enter the wd

setwd(wd.lf)



#####

#### Set up the data for the analyses ####



### Identify the clade being worked with
### NOTE: this reads all the correct data, tunes the model parameters, 
### and adjusts the way outputs are saved/backed up further down


set.group <- "coelacanths"
#set.group <- "lungfishes"
#set.group <- "holosteans"



### Import Data & set the model fitting bounds
### (lower/upper, in that order) 
### NOTE: bounds are set as lower and upper in that order and they were
### determined by running the analyses and seeing if the bootstraps were
### repeatedly reaching the edge of the search space. They were then 
### adjusted until this was not the case. NA = no bounds; NULL = default
### bounds (should be called as 'list(NULL)' or else it breaks).


if(exists("set.group") == FALSE || is.null(set.group) == TRUE){
  
  print(
  "You haven't let the script know what clade you're working with yet!"
  )
  
} else if (set.group == "coelacanths") {
  
  # Maximum Clade Credibility tree (pruned to match shape data)
  tree.mcc.shapeanalyses.coelacanths <- 
    readRDS("tree_mcc_coelacanth_shapeanalyses.RDS")
  tree.mcc.shapeanalyses.reduced <- 
    tree.mcc.shapeanalyses.coelacanths
  
  # Posterior distribution of trees outputted from Bayesian analysis
  trees.posterior.shapeanalyses.coelacanths <- 
    readRDS("trees_posterior_coelacanth_shapeanalyses.RDS")
  trees.posterior.shapeanalyses.reduced <- 
    trees.posterior.shapeanalyses.coelacanths
  
  # Read shape data
  shapes.coelacanths <- readRDS("data_shape_coelacanth.RDS")
  shapes <- shapes.coelacanths
  
  # Upper and lower bounds for the mvgls analyses
  bounds <- list()  # 1 is always EB, 2 is always BM, 3 is always OU
  bounds[[1]] <- list(-50,10)
  bounds[[2]] <- list(NULL,NULL)
  bounds[[3]] <- list(-25,265) # finnicky up, careful changing it
  
} else if (set.group == "lungfishes") {
  
  # Maximum Clade Credibility tree (pruned to match shape data)
  tree.mcc.shapeanalyses.lungfishes <- 
    readRDS("tree_mcc_lungfish_shapeanalyses.RDS")
  tree.mcc.shapeanalyses.reduced <- 
    tree.mcc.shapeanalyses.lungfishes
  
  # Posterior distribution of trees outputted from Bayesian analysis
  trees.posterior.shapeanalyses.lungfishes <- 
    readRDS("trees_posterior_lungfish_shapeanalyses.RDS")
  trees.posterior.shapeanalyses.reduced <- 
    trees.posterior.shapeanalyses.lungfishes
  
  # Read shape data
  shapes.lungfishes <- readRDS("data_shape_lungfish.RDS")
  shapes <- shapes.lungfishes
  
  # Upper and lower bounds for the mvgls analyses
  bounds <- list()
  bounds[[1]] <- list(-50,10)
  bounds[[2]] <- list(NULL,NULL)
  bounds[[3]] <- list(-25,250) # very finnicky up, careful changing it
  
} else if (set.group == "holosteans") {
  
  # Maximum Clade Credibility tree (pruned to match shape data)
  tree.mcc.shapeanalyses.holosteans <- 
    readRDS("tree_mcc_holostean_shapeanalyses.RDS")
  tree.mcc.shapeanalyses.reduced <- 
    tree.mcc.shapeanalyses.holosteans
  
  # Posterior distribution of trees outputted from Bayesian analysis
  trees.posterior.shapeanalyses.holosteans <- 
    readRDS("trees_posterior_holostean_shapeanalyses.RDS")
  trees.posterior.shapeanalyses.reduced <- 
    trees.posterior.shapeanalyses.holosteans
  
  # Read shape data
  shapes.holosteans <- readRDS("data_shape_holostean.RDS")
  shapes <- shapes.holosteans
  
  # Upper and lower bounds for the mvgls analyses
  bounds <- list()
  bounds[[1]] <- list(-10,10)
  bounds[[2]] <- list(NULL,NULL)
  bounds[[3]] <- list(NULL,NULL)
  
}



#####

#### Transform and analyze the 2D shape data ####



### Generalized Procrustes Analysis (GPA)


## GPA transform the shape data

gpa <- gpagen(shapes)


## Store just the Procrustes-aligned coordinates
## NOTE: multiply by -1 to display them "right side" up visually for 
## supplement

procrustes.coordinates <- gpa$coords * -1



### Principal Components Analyses (PCA)


## Standard PCA

pca <- gm.prcomp(procrustes.coordinates)


## Phylogenetically-aligned Components Analysis (PaCA) for supplement

pca.phylo <- gm.prcomp(
  procrustes.coordinates, phy = tree.mcc.shapeanalyses.reduced
  )



### Extract PC values to be used in analyses


## Determine whether you want axes summarizing a % of the variance or 
## all axes
## NOTE: don't use fractional values (e.g., 0.95)! Write the full number
## you want. 

p.variance <- 100 # 100% of the variance (all axes)
#p.variance <- 95 # 95% of the variance (LOWER THAN 100 NOT USED HERE)


## Extract those axes from the PCA

if (p.variance == 100) {
  
  pca.shapevalues.matrix <- pca$x
  
} else if ( 0 < p.variance && p.variance < 100) {
  
  # Calculate the cumulative proportion of variation
  # NOTE: this has to be done like this because $importance is
  # deprecated. It was adapted from here: https://stackoverflow.com
  # /questions/23866765/getting-cumulative-proportion-in-pca
  pca.cumprop <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
  
  # Identify the axes required to summarize >95 of the variance
  pca.axes <- min(which(pca.cumprop > (p.variance/100)))
  
  # Store the PC scores (eigenvalues) to use as inputs for the rate 
  # shift analyses
  pca.shapevalues.matrix <- pca$x[,1:pca.axes]
  
} else {
  
  paste("You didn't specify a % variance as a number between 1 - 100!")
  
}



#####

#### CHECKPOINT: Export transformed procrustes and PC data ####



### Store the outputs according to the group you're testing and 
### standardize


if (set.group == "coelacanths") {
  
  # Back up the pca matrix according to the axis number selected
  pca.shapevalues.matrix.coelacanths <- pca.shapevalues.matrix
  saveRDS(pca.shapevalues.matrix
    , 'pca_coelacanth_shapevalues_matrix_100.RDS'
    )
  
} else if (set.group == "lungfishes") {
  
  # Back up the pca matrix according to the axis number selected
  pca.shapevalues.matrix.lungfishes <- pca.shapevalues.matrix
  saveRDS(pca.shapevalues.matrix
    , 'pca_lungfish_shapevalues_matrix_100.RDS'
    )
  
} else if (set.group == "holosteans") {

  # Back up the pca matrix according to the axis number selected
  pca.shapevalues.matrix.holosteans <- pca.shapevalues.matrix
  saveRDS(pca.shapevalues.matrix
    , 'pca_holostean_shapevalues_matrix_100.RDS'
    )
  
} else {
  
  print(
    "Something went wrong. Have you set the clade you're working with?"
    )
  
}



### Do you want to read in the checkpoint files?
### NOTE: use this if you don't want to rerun analyses

#load.checkpoint <- TRUE
load.checkpoint <- FALSE


if(
  exists("load.checkpoint") == FALSE || is.null(load.checkpoint) == TRUE
  ){
  
  print(
  "You haven't let the script know if you want to load from this 
  checkpoint!"
  )
  
} else if (load.checkpoint == FALSE) {
  
  print("You haven't opted to load any saved files.")
  
} else if (load.checkpoint == TRUE) {
  
  if (set.group == "coelacanths") {
    
    pca.shapevalues.matrix <- 
      readRDS('pca_coelacanth_shapevalues_matrix_100.RDS')
    
  } else if (set.group == "lungfishes") {
    
    pca.shapevalues.matrix <- 
      readRDS('pca_lungfish_shapevalues_matrix_100.RDS')
    
  } else if (set.group == "holosteans") {
    
    pca.shapevalues.matrix <- 
      readRDS('pca_holostean_shapevalues_matrix_100.RDS')
    
  }
  
}



#####

#### Prepare MCC and posterior trees for analyses ####



### Burn-in and sample trees from the posterior


## Remove 'burnin.p' % of the posterior trees as burn-in to match how 
## the MCC tree was calculated

# Set the burn-in %
burnin.p <- 0.1 # 10% to align with how the MCC was calculated
#burnin.p <- n # set another % in decimal value

# Calculate 'burnin.p' % of the sample
burnin <- 
  round(burnin.p * length(trees.posterior.shapeanalyses.reduced))

# Remove 'burnin.p' % of the sample as burn-in
trees.modeltests <- 
  trees.posterior.shapeanalyses.reduced[
    burnin:length(trees.posterior.shapeanalyses.reduced)
    ]


## Sample 'nsample' trees from the post burn-in posterior trees

# Set seed for reproducibility
set.seed(67)

# Set the total trees to be sampled
nsample <- 100 # manual
#nsample <- length(trees.posterior.shapeanalyses.reduced)#automatic, all

# Sample 'nsample' trees from the post burn-in posterior trees
trees.modeltests.sample <- sample(trees.modeltests, nsample) # random
#trees.modeltests.sample <- trees.modeltests[100:199] # manual

# Set the class to 'multiPhylo' so nothing breaks
class(trees.modeltests.sample) <- "multiPhylo"



### Transform the posterior and MCC trees to unit height
### NOTE: can be done either separately (i.e., using 'max.node.height' 
### from the sampled posterior and the 'max.node.height' from the MCC 
### to rescale each subset of trees individually) or all together (i.e.,
### using the 'max.node.height' from the sampled posterior + the 
### 'max.node.height' from the MCC to rescale all trees). You should
### use the latter if you want every single tree used from now on to be
### scaled by the same factor.


## Select the scaling type

scaling <- "together"
#scaling <- "separate"


if (scaling == "together") {
  
  # Extract the max node heights from the posterior trees
  max.node.height <- 
    lapply(trees.modeltests.sample, function(x) max(nodeHeights(x)))
  
  # Find the max node height across the *entire dataset*
  max.node.height <- max(
    c(unlist(max.node.height)
      , max(nodeHeights(tree.mcc.shapeanalyses.reduced))
      )
    )
  
  # Create an empty list object to receive the rescaled posterior trees
  trees.modeltests.sample.unitheight <- list()
  
  # Loop through the sampled trees to rescale them
  for (i in 1:length(trees.modeltests.sample)){
   trees.modeltests.sample.unitheight[[i]] <- 
    trees.modeltests.sample[[i]]
   trees.modeltests.sample.unitheight[[i]]$edge.length <-
   trees.modeltests.sample.unitheight[[i]]$edge.length / max.node.height
  }
  
  # Set the class to 'multiPhylo' so nothing breaks
  class(trees.modeltests.sample.unitheight) <- "multiPhylo"
  
  # Make a copy of the MCC tree
  tree.modeltests.unitheight <- tree.mcc.shapeanalyses.reduced
  
  # Rescale the MCC tree *using the largest max node height out of all*
  # *the trees*
  tree.modeltests.unitheight$edge.length <-
    tree.modeltests.unitheight$edge.length / max.node.height
  
} else if (scaling == "separate") {
  
  # Extract the max node heights from the posterior trees
  max.node.height <- 
    lapply(trees.modeltests.sample, function(x) max(nodeHeights(x)))
  
  # Find the max node height across the *posterior sample*
  max.node.height <- max(unlist(max.node.height))
  
  # Create an empty list object to receive the rescaled posterior trees
  trees.modeltests.sample.unitheight <- list()
  
  # Loop through the sampled trees to rescale them
  for (i in 1:length(trees.modeltests.sample)){
  trees.modeltests.sample.unitheight[[i]] <- 
    trees.modeltests.sample[[i]]
  trees.modeltests.sample.unitheight[[i]]$edge.length <-
   trees.modeltests.sample.unitheight[[i]]$edge.length / max.node.height
  }
  
  # Set the class to 'multiPhylo' so nothing breaks
  class(trees.modeltests.sample.unitheight) <- "multiPhylo"
  
  # Make a copy of the MCC tree
  tree.modeltests.unitheight <- tree.mcc.shapeanalyses.reduced
  
  # Rescale the MCC tree *using it's own max node height as the scaling*
  # *factor*
  tree.modeltests.unitheight$edge.length <-
    tree.modeltests.unitheight$edge.length / max(
                                nodeHeights(tree.modeltests.unitheight)
                                )
  
} else {
  
  paste("You haven't set the scaling correctly, rerun and try again!")
  
}



#####

#### Fit mvgls models on the MCC tree ####
#### NOTE: uses the MCC tree and PC data as inputs


### Fit EB, BM, and OU models
### NOTE: uses the PC scores representing the selected % variance and a 
### pruned, unit height scaled MCC tree

system.time(
fit.EB.pca <- mvgls(
  pca.shapevalues.matrix ~ 1
  , tree = tree.modeltests.unitheight
  , scale.height = F
  , model = "EB"
  , method = "PL-LOOCV" #default
  #, method = "LL"
  #, method = "H&L"
  #, method = "Mahalanobis"
  #, REML = FALSE #use ML to estimate parameter
  , target = "unitVariance"
  #, target = "Variance"
  #, target = "null" #only with penalty = "RidgeAlt"
  , error = T
  , penalty ="RidgeArch" #default
  #, penalty = "RidgeAlt"
  , upper = bounds[[1]][[2]]
  , lower = bounds[[1]][[1]]
  , ncores = 20
)
)

system.time(
fit.BM.pca <- mvgls(
  pca.shapevalues.matrix ~ 1
  , tree = tree.modeltests.unitheight
  , scale.height = F
  , model = "BM"
  , method = "PL-LOOCV" #default
  #, method = "LL"
  #, method = "H&L"
  #, method = "Mahalanobis"
  #, REML = FALSE #use ML to estimate parameter
  , target = "unitVariance"
  #, target = "Variance"
  #, target = "null" #only with penalty = "RidgeAlt"
  , error = T
  , penalty ="RidgeArch" #default
  #, penalty = "RidgeAlt"
  , upper = bounds[[2]][[2]]
  , lower = bounds[[2]][[1]]
  , ncores = 20
)
)

system.time(
fit.OU.pca <- mvgls(
  pca.shapevalues.matrix ~ 1
  , tree = tree.modeltests.unitheight
  , scale.height = F
  , model = "OU"
  , method = "PL-LOOCV" #default
  #, method = "LL"
  #, method = "H&L"
  #, method = "Mahalanobis"
  #, REML = FALSE #use ML to estimate parameter
  , target = "unitVariance"
  #, target = "Variance"
  #, target = "null" #only with penalty = "RidgeAlt"
  , error = T
  , penalty ="RidgeArch" #default
  #, penalty = "RidgeAlt"
  , upper = bounds[[3]][[2]]
  , lower = bounds[[3]][[1]]
  , ncores = 20
)
)



### Estimate model weights for the MCC tree


## Join the results of all the model tests

fits.modeltests.mcc.pca <- NULL
fits.modeltests.mcc.pca[[1]] <- fit.EB.pca # 1 is always EB
fits.modeltests.mcc.pca[[2]] <- fit.BM.pca # 2 is always BM
fits.modeltests.mcc.pca[[3]] <- fit.OU.pca # 3 is always OU


## Calculate the Generalized Information Criterion (GIC) scores for
## each model

fits.modeltests.mcc.pca.gic <- 
  c(EB = GIC(fits.modeltests.mcc.pca[[1]])$GIC
    , BM = GIC(fits.modeltests.mcc.pca[[2]])$GIC
    , OU = GIC(fits.modeltests.mcc.pca[[3]])$GIC
    )


## Calculate the GIC weights (GICw) for each model test using the
## 'aicw' function
## NOTE: although we're using a different information criterion 
## (GIC vs AIC), the formula to calculate the weights doesn't change,
## so we use the standard 'aicw' function from the 'geiger' package to
## calculate the GICw.

fits.modeltests.mcc.pca.gic.weights <- 
  geiger::aicw(fits.modeltests.mcc.pca.gic)



#####

#### CHECKPOINT: Export the MCC model test results ####



### Store the outputs according to the group you're testing and
### standardize


if (set.group == "coelacanths") {
  
  # Back up the MCC model fits
  fits.modeltests.mcc.pca.coelacanths <- fits.modeltests.mcc.pca

  saveRDS(fits.modeltests.mcc.pca
    , 'fittree_coelacanth_pcamodeltests_100.RDS'
    )
  
  # Back up the MCC model fit weights
  fits.modeltests.mcc.pca.gic.weights.coelacanths <- 
    fits.modeltests.mcc.pca.gic.weights

  saveRDS(fits.modeltests.mcc.pca.gic.weights
    , 'fittree_coelacanth_pcamodeltests_gicw_100.RDS'
    )
  
  
} else if (set.group == "lungfishes") {
  
  # Back up the MCC model fits
  fits.modeltests.mcc.pca.lungfishes <- fits.modeltests.mcc.pca

  saveRDS(fits.modeltests.mcc.pca
    , 'fittree_lungfish_pcamodeltests_100.RDS'
    )
  
  # Back up the MCC model fit weights
  fits.modeltests.mcc.pca.gic.weights.lungfishes <- 
    fits.modeltests.mcc.pca.gic.weights

  saveRDS(fits.modeltests.mcc.pca.gic.weights
    , 'fittree_lungfish_pcamodeltests_gicw_100.RDS'
    )
  
  
} else if (set.group == "holosteans") {
  
  # Back up the MCC model fits
  fits.modeltests.mcc.pca.holosteans <- fits.modeltests.mcc.pca

  saveRDS(fits.modeltests.mcc.pca
    , 'fittree_holostean_pcamodeltests_100.RDS'
    )
  
  # Back up the MCC model fit weights
  fits.modeltests.mcc.pca.gic.weights.holosteans <- 
    fits.modeltests.mcc.pca.gic.weights

  saveRDS(fits.modeltests.mcc.pca.gic.weights
    , 'fittree_holostean_pcamodeltests_gicw_100.RDS'
    )
  
  
} else {
  
  print("Something went wrong. Have you set your clade correctly?")
  
}



### Do you want to read in the checkpoint files?
### NOTE: use this if you don't want to rerun analyses

#load.checkpoint <- TRUE
load.checkpoint <- FALSE


if(
  exists("load.checkpoint") == FALSE || is.null(load.checkpoint) == TRUE
  ){
  
  print("You haven't set if you want to load from this checkpoint!")
  
} else if (load.checkpoint == FALSE) {
  
  print("You haven't opted to load any saved files.")
  
} else if (load.checkpoint == TRUE) {
  
  if (set.group == "coelacanths") {
    
    fits.modeltests.mcc.pca <- 
      readRDS('fittree_coelacanth_pcamodeltests_100.RDS')
      
    fits.modeltests.mcc.pca.gic.weights <- 
      readRDS('fittree_coelacanth_pcamodeltests_gicw_100.RDS')
    
  } else if (set.group == "lungfishes") {
    
    fits.modeltests.mcc.pca <- 
      readRDS('fittree_lungfish_pcamodeltests_100.RDS')

    fits.modeltests.mcc.pca.gic.weights <- 
      readRDS('fittree_lungfish_pcamodeltests_gicw_100.RDS')
    
  } else if (set.group == "holosteans") {
    
    fits.modeltests.mcc.pca <- 
      readRDS('fittree_holostean_pcamodeltests_100.RDS')

    fits.modeltests.mcc.pca.gic.weights <- 
      readRDS('fittree_holostean_pcamodeltests_gicw_100.RDS')
    
  }
  
}



#####

#### Fit mvgls models on a sample from the posterior distribution ####
#### NOTE: uses a random sample from the posterior distribution of trees
#### and PC data as inputs



### parallelize via the 'future' package


## Total number of cores to use in parallel

#n.cores <- 4 # laptops
#n.cores <- 8 # powerful laptops/desktops
#n.cores <- 16 # workstations/servers
n.cores <- availableCores() - 
  ceiling(availableCores() * 0.05) # automatic using 95% available cores
#n.cores <- N # manual, edit here


## Set up the parallelization session using 'n.cores' number of cores

plan(multisession, workers = n.cores)



### Fit EB, BM, and OU models
### NOTE: uses the PC scores representing the selected % variance and 
### pruned, unit height scaled posterior trees. Additionally, some tests
### may return errors, so 'tryCatch' has been included to insert NAs in
### those places (helps to filter them)


system.time(
fit.EB.pca.posterior <- 
  future_lapply(1:length(trees.modeltests.sample.unitheight)
  , function(x){
  tryCatch(mvgls(
    pca.shapevalues.matrix ~ 1
    , tree = trees.modeltests.sample.unitheight[[x]]
    , scale.height = F
    , model = "EB"
    , method = "PL-LOOCV" #default
    #, method = "LL"
    #, method = "H&L"
    #, method = "Mahalanobis"
    #, REML = FALSE #use ML to estimate parameter
    , target = "unitVariance"
    #, target = "Variance"
    #, target = "null" #only with penalty = "RidgeAlt"
    , error = T
    , penalty ="RidgeArch" #default
    #, penalty = "RidgeAlt"
    , upper = bounds[[1]][[2]]
    , lower = bounds[[1]][[1]]
    #, ncores = 20
    )
    ,error = function(e) return(NA), warning = function(w) return(NA)
    )
}, future.seed=TRUE)
)

system.time(
fit.BM.pca.posterior <- 
  future_lapply(1:length(trees.modeltests.sample.unitheight)
  , function(x){
  tryCatch(mvgls(
    pca.shapevalues.matrix ~ 1
    , tree = trees.modeltests.sample.unitheight[[x]]
    , scale.height = F
    , model = "BM"
    , method = "PL-LOOCV" #default
    #, method = "LL"
    #, method = "H&L"
    #, method = "Mahalanobis"
    #, REML = FALSE #use ML to estimate parameter
    , target = "unitVariance"
    #, target = "Variance"
    #, target = "null" #only with penalty = "RidgeAlt"
    , error = T
    , penalty ="RidgeArch" #default
    #, penalty = "RidgeAlt"
    , upper = bounds[[2]][[2]]
    , lower = bounds[[2]][[1]]
    #, ncores = 20
    )
    ,error = function(e) return(NA), warning = function(w) return(NA)
    )
}, future.seed=TRUE)
)

system.time(
fit.OU.pca.posterior <- 
  future_lapply(1:length(trees.modeltests.sample.unitheight)
  , function(x){
  tryCatch(mvgls(
    pca.shapevalues.matrix ~ 1
    , tree = trees.modeltests.sample.unitheight[[x]]
    , scale.height = F
    , model = "OU"
    , method = "PL-LOOCV" #default
    #, method = "LL"
    #, method = "H&L"
    #, method = "Mahalanobis"
    #, REML = FALSE #use ML to estimate parameter
    , target = "unitVariance"
    #, target = "Variance"
    #, target = "null" #only with penalty = "RidgeAlt"
    , error = T
    , penalty ="RidgeArch" #default
    #, penalty = "RidgeAlt"
    , upper = bounds[[3]][[2]]
    , lower = bounds[[3]][[1]]
    #, ncores = 20
    )
    ,error = function(e) return(NA), warning = function(w) return(NA)
    )
}, future.seed=TRUE)
)



### Estimate model weights for the posterior tree


## Create functions to compute GICs on a list of model fits

# 'tryGIC' to calculate GICs with 'tryCatch' outputting NAs in case of
# errors
tryGIC <- function(fit){
  result <- 
    tryCatch(GIC(fit)$GIC, warning=function(w) NA, error=function(e) NA)
  return(result)
}

# 'compute_GICs' to 'lapply' 'tryGIC' across all model fits while
# removing entries that return a warning or error
compute_GICs <- function(fits) {
  tmp <- lapply(fits, tryGIC)
  fits.GICs <- unlist(tmp)
  return(fits.GICs)
}


## Calculate the GIC scores for each model test
## NOTE: all GIC vectors should be the same length, with NAs where there
## was either an error with model fitting or with estimating GIC

fit.EB.pca.posterior.GICs <- compute_GICs(fit.EB.pca.posterior)
fit.BM.pca.posterior.GICs <- compute_GICs(fit.BM.pca.posterior)
fit.OU.pca.posterior.GICs <- compute_GICs(fit.OU.pca.posterior)


## Calculate GICw for each model test using the 'aicw' function

# Make a new data frame with the GIC scores across all trees
fit.all.pca.posterior.GICs <- 
  data.frame(
    EB = fit.EB.pca.posterior.GICs
    , BM = fit.BM.pca.posterior.GICs
    , OU = fit.OU.pca.posterior.GICs
    )

# Filter out any rows with NAs
fit.all.pca.posterior.GICs <- 
  fit.all.pca.posterior.GICs[
    complete.cases(fit.all.pca.posterior.GICs),
    ]

# Estimate GICw across all trees
fit.all.pca.posterior.GICs.weights <- 
  lapply(1:length(fit.all.pca.posterior.GICs[,1]), function(x){
    geiger::aicw(
      setNames(
        as.numeric(fit.all.pca.posterior.GICs[x,]),c("EB", "BM", "OU"))
      )
  }
)

# Replace GIC scores in order from 'fit.all.pca.posterior.GICs' with
# their respective GICw calculated in the previous step
weight.data <- fit.all.pca.posterior.GICs

for (i in 1:length(fit.all.pca.posterior.GICs.weights)){
  for (p in 1:3){
    weight.data[i,p] <- fit.all.pca.posterior.GICs.weights[[i]]$w[p]
  }
}

fit.all.pca.posterior.GICs.weights <- weight.data



#####

#### CHECKPOINT: Export the posterior model test results ####



### Store the outputs according to the group you're testing and 
### standardize


if (set.group == "coelacanths") {
  
  # Back up the GICs for all the posterior fits
  fit.all.pca.posterior.GICs.coelacanths <- fit.all.pca.posterior.GICs
  
  saveRDS(fit.all.pca.posterior.GICs
    , 'fittrees_coelacanth_pcamodeltests_100.RDS'
    )
  
  # Back up the GIC weights for all the posterior fits
  fit.all.pca.posterior.GICs.weights.coelacanths <- 
    fit.all.pca.posterior.GICs.weights
  
  saveRDS(fit.all.pca.posterior.GICs.weights
    , 'fittrees_coelacanth_pcamodeltests_gicw_100.RDS'
    )
  
} else if (set.group == "lungfishes") {
  
  # Back up the GICs for all the posterior fits
  fit.all.pca.posterior.GICs.lungfishes <- fit.all.pca.posterior.GICs
  saveRDS(fit.all.pca.posterior.GICs
    , 'fittrees_lungfish_pcamodeltests_100.RDS'
    )
  
  # Back up the GIC weights for all the posterior fits
  fit.all.pca.posterior.GICs.weights.lungfishes <- 
    fit.all.pca.posterior.GICs.weights
  
  saveRDS(fit.all.pca.posterior.GICs.weights
    , 'fittrees_lungfish_pcamodeltests_gicw_100.RDS'
    )
  
} else if (set.group == "holosteans") {
  
  # Back up the GICs for all the posterior fits
  fit.all.pca.posterior.GICs.holosteans <- fit.all.pca.posterior.GICs
  
  saveRDS(fit.all.pca.posterior.GICs
    , 'fittrees_holostean_pcamodeltests_100.RDS'
    )
  
  # Back up the GIC weights for all the posterior fits
  fit.all.pca.posterior.GICs.weights.holosteans <- 
    fit.all.pca.posterior.GICs.weights
  
  saveRDS(fit.all.pca.posterior.GICs.weights
    , 'fittrees_holostean_pcamodeltests_gicw_100.RDS'
    )
  
} else {
  
  print(
    "Something went wrong. Have you set the clade you're working with?"
    )
  
}



### Do you want to read in the checkpoint files?
### NOTE: use this if you don't want to rerun analyses


#load.checkpoint <- TRUE
load.checkpoint <- FALSE


if(
  exists("load.checkpoint") == FALSE || is.null(load.checkpoint) == TRUE
  ){
  
  print("You haven't set if you want to load from this checkpoint!")
  
} else if (load.checkpoint == FALSE) {
  
  print("You haven't opted to load any saved files.")
  
} else if (load.checkpoint == TRUE) {
  
  if (set.group == "coelacanths") {
    
    fit.all.pca.posterior.GICs <- 
      readRDS('fittrees_coelacanth_pcamodeltests_100.RDS')
    fit.all.pca.posterior.GICs.weights <- 
      readRDS('fittrees_coelacanth_pcamodeltests_gicw_100.RDS')
    
  } else if (set.group == "lungfishes") {
    
    fit.all.pca.posterior.GICs <- 
      readRDS('fittrees_lungfish_pcamodeltests_100.RDS')
    fit.all.pca.posterior.GICs.weights <- 
      readRDS('fittrees_lungfish_pcamodeltests_gicw_100.RDS')
    
  } else if (set.group == "holosteans") {
    
    fit.all.pca.posterior.GICs <- 
      readRDS('fittrees_holostean_pcamodeltests_100.RDS')
    fit.all.pca.posterior.GICs.weights <- 
      readRDS('fittrees_holostean_pcamodeltests_gicw_100.RDS')
    
  }
  
}



#####

#### Plot the MCC tree and posterior trees model test results ####



### Set up to automatically generate the plot titles
### NOTE: this uses the group you have as your 'set.group' and inserts 
### it into the title from that object

plot.title <- paste(
  "Model Fit across Posterior Trees\n("
  , set.group
  , sep = ""
  )



### Visualize distribution of GIC scores and weights


## Select the display style
## NOTE: you can select whether you want to display the figures side by
## side (i.e., changing par to 1,2) or individually (i.e., the default 
## par of 1,1)

display <- "joined"
#display <- "individual"


if (display == "joined") {
  
  par(mfrow=c(1,2))
  
} else if (display == "individual") {
  
  par(mfrow=c(1,1))
  
} else {
  
  paste("You haven't selected the display, rerun and try again!")
  
}


## Generate a plot of GIC scores (re-ordered to BM, OU, EB/ACDC)

vioplot(
  list(
    BM = fit.all.pca.posterior.GICs$BM
    , OU = fit.all.pca.posterior.GICs$OU
    , ACDC = fit.all.pca.posterior.GICs$EB
  )
  , ylab = "GIC scores"
  , main = paste(plot.title, "GIC scores)")
)

points(
  fits.modeltests.mcc.pca.gic.weights$fit[c(2,3,1)] # BM, OU, EB
  , col = 'red'
  , pch = 19
  , cex = 2
)


## Generate a plot of GIC weights (re-ordered to BM, OU, EB/ACDC)

vioplot(
  list(
    BM = fit.all.pca.posterior.GICs.weights$BM
    , OU = fit.all.pca.posterior.GICs.weights$OU
    , ACDC = fit.all.pca.posterior.GICs.weights$EB
  )
  , ylab = "GIC weights"
  , main = paste(plot.title, "GIC weights)")
  , ylim = c(0,1)
  )

points(
  fits.modeltests.mcc.pca.gic.weights$w[c(2,3,1)] #BM, OU, EB
  , col = 'red'
  , pch = 19
  , cex = 2
)


## Reset the display back to the default

par(mfrow=c(1,1))


#####
