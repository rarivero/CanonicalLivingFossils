#
##
#### BOOTSTRAPPED EVOLUTIONARY MODEL FITTING 
#### USING SHAPE DATA + CONSENSUS TREES
##
#

#### Script by Rafael A. Rivero Vega & Jacob S. Berv
#### Based off original and contributed to by Julien Clavel, with
#### further details in the manuscript methods and acknowledgements.




#### Load the required libraries ####

if(!'mvMORPH'%in%installed.packages()
  ) install.packages("mvMORPH"); library(mvMORPH)

if(!'future.apply'%in%installed.packages()
  ) install.packages("future.apply"); library(future.apply)

if(!'HDInterval'%in%installed.packages()
  ) install.packages("HDInterval"); library(HDInterval)

#rm(list = ls()) # if you need to clean up



#####

#### Specify path to working directory ####

#wd.lf <- readClipboard() # if you copied and want to paste the wd
#wd.lf <- setwd(choose.dir()) # if you want to popup select the wd 
#wd.lf <- "~/path/here" # if you want to manually enter the wd

setwd(wd.lf)



#####

#### Import functions ####


### Specify full path(s) and/or links to the function(s) below
### NOTE: future_mvSIM.R is an updated version of mvSIM from the mvMORPH 
### package that works with futures for parallelization. More info on
### all of this in the relevant section below.

source("./future_mvSIM.R")
#source(paste(script1))



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
    tree <- readRDS('tree_mcc_coelacanth_shapeanalyses.RDS')
    
    # Matrix representing 100% of the PCs
    pca.shapevalues.matrix <- 
      readRDS('pca_coelacanth_shapevalues_matrix_100.RDS')
    
    # Upper and lower bounds for the mvgls analyses
    bounds <- list()  # 1 is always EB, 2 is always BM, 3 is always OU
    bounds[[1]] <- list(-50,10)
    bounds[[2]] <- list(NULL,NULL)
    bounds[[3]] <- list(-25,265) # finnicky up, careful changing it
    
  } else if (set.group == "lungfishes") {
      
    # Maximum Clade Credibility tree (pruned to match shape data)
    tree <- readRDS('tree_mcc_lungfish_shapeanalyses.RDS')
    
    # Matrix representing 100% of the PCs
    pca.shapevalues.matrix <- 
      readRDS('pca_lungfish_shapevalues_matrix_100.RDS')
    
    # Upper and lower bounds for the mvgls analyses
    bounds <- list()
    bounds[[1]] <- list(-50,10)
    bounds[[2]] <- list(NULL,NULL)
    bounds[[3]] <- list(-25,250) # very finnicky up, careful changing it
      
  } else if (set.group == "holosteans") {
    
    # Maximum Clade Credibility tree (pruned to match shape data)
    tree <- readRDS('tree_mcc_holostean_shapeanalyses.RDS')
    
    # Matrix representing 100% of the PCs
    pca.shapevalues.matrix <- 
      readRDS('pca_holostean_shapevalues_matrix_100.RDS')
    
    # Upper and lower bounds for the mvgls analyses
    bounds <- list()
    bounds[[1]] <- list(-10,10)
    bounds[[2]] <- list(NULL,NULL)
    bounds[[3]] <- list(NULL,NULL)
      
  }



#####

#### Prepare MCC tree for analyses ####



### Transform the MCC tree to unit height

tree.unitheight <- tree
tree.unitheight$edge.length <- 
  tree.unitheight$edge.length/max(nodeHeights(tree))



#####

#### Fit mvgls models on the MCC tree ####
#### NOTE: uses the MCC tree and PC data as inputs



### Early Burst (EB) model fit

if(exists("set.group") == FALSE || is.null(set.group) == TRUE){
  
  print(
  "You haven't let the script know what clade you're working with yet!"
  )
  
} else if (set.group == "holosteans") {
  
  system.time(
  fit.EB.pca <- mvgls(
    pca.shapevalues.matrix ~ 1
    , tree = tree.unitheight
    , scale.height = F
    , model = "EB"
    , method = "PL-LOOCV" #default
    #, method = "LL"
    #, method = "H&L"
    #, method = "Mahalanobis"
    , REML = FALSE #FALSE = use ML to estimate parameter
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
  
} else {
  
  system.time(
  fit.EB.pca <- mvgls(
    pca.shapevalues.matrix ~ 1
    , tree = tree.unitheight
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

}



### Brownian Motion (BM) model fit

system.time(
fit.BM.pca <- mvgls(
  pca.shapevalues.matrix ~ 1
  , tree = tree.unitheight
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



### Ornstein-Uhlenbeck (OU) model fit

system.time(
fit.OU.pca <- mvgls(
  pca.shapevalues.matrix ~ 1
  , tree = tree.unitheight
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



### Store the outputs according to the group you're testing and 
### standardize


if (set.group == "coelacanths") {
  
  fit.pca.coelacanths <- NULL
  fit.pca.coelacanths[[1]] <- fit.EB.pca # 1 is always EB
  fit.pca.coelacanths[[2]] <- fit.BM.pca # 2 is always BM
  fit.pca.coelacanths[[3]] <- fit.OU.pca # 3 is always OU
  fits.pca <- fit.pca.coelacanths
  
} else if (set.group == "lungfishes") {
  
  fit.pca.lungfishes <- NULL
  fit.pca.lungfishes[[1]] <- fit.EB.pca
  fit.pca.lungfishes[[2]] <- fit.BM.pca
  fit.pca.lungfishes[[3]] <- fit.OU.pca
  fits.pca <- fit.pca.lungfishes
  
} else if (set.group == "holosteans") {
  
  fit.pca.holosteans <- NULL
  fit.pca.holosteans[[1]] <- fit.EB.pca
  fit.pca.holosteans[[2]] <- fit.BM.pca
  fit.pca.holosteans[[3]] <- fit.OU.pca
  fits.pca <- fit.pca.holosteans
  
} else {
  
  print("Something went wrong. Have you set your clade correctly?")
  
}



#####

#### Simulate traits using mvSIM ####
#### NOTE: uses the MCC tree and PC data as inputs



### Set seed for reproducibility

set.seed(46)



### Set the number of simulations to perform

n.sim <- 100



### Adjust 'future' parallelization options
### NOTE: Using mclapply in the parallel package doesn't work on all 
### systems, so it is recommended to use a solution like the
### future.apply package instead https://stackoverflow.com/questions/
### 62080053/mc-cores-1-is-not-supported-on-windows
### or doParallel + foreach
### https://cran.r-project.org/web/packages/doParallel/vignettes/
### gettingstartedParallel.pdf
### https://privefl.github.io/blog/a-guide-to-parallelism-in-r/.
### In these analyses, we've optimized everything to function using the
### future package.


## Pass the original environment so the future_mvSIM function that was
## imported at the beginning of this script works with all of the
## same functions the original mvSIM does

environment(future_mvSIM) <- environment(mvSIM)


## Toggle parallelization on or off for the simulations

parallelize <- TRUE
#parallelize <- FALSE


## Apply the appropriate parallelization plan depending on the selection 
## above

if(exists("parallelize") == FALSE || is.null(parallelize) == TRUE){
  
  print("You haven't specified whether or not to parallelize mvSIM!")
  
} else if (parallelize == FALSE) {
  
  plan(sequential)
  
} else if (parallelize == TRUE) {
  
  n.cores <- availableCores() - 
    ceiling(
    availableCores() * 0.05) # automatic using 95% available cores
  #n.cores <- 10 # manual, edit here
  
  plan(multisession, workers = n.cores)
  
}



### Run the simulations
### NOTE: theta = starting values for simulation (x$coefficients)


## EB simulations
system.time(
sim.EB <- future_mvSIM(
  fits.pca[[1]]$corrSt$phy
  , nsim = n.sim
  , model = "EB"
  , param = list(
    sigma = fits.pca[[1]]$sigma$Pinv
    , theta = fits.pca[[1]]$coefficients
    , beta = fits.pca[[1]]$param
    )
  )

)


## BM simulations
system.time(
sim.BM <- future_mvSIM(
  fits.pca[[2]]$corrSt$phy
  , nsim = n.sim
  , model = "BM1"
  , param = list(
    sigma = fits.pca[[2]]$sigma$Pinv
    , theta = fits.pca[[2]]$coefficients
    )
  )

)


## OU simulations
## NOTE: ignore decomposition errors, it sometimes adjusts method 
## automatically
system.time(
sim.OU <- future_mvSIM(
  fits.pca[[3]]$corrSt$phy
  , nsim = n.sim
  , model = "OU1"
  , param = list(
    sigma = fits.pca[[3]]$sigma$Pinv
    , theta = fits.pca[[3]]$coefficients
    , alpha = diag(fits.pca[[3]]$param
            , length(coefficients(fits.pca[[3]])
            )
      )
    )
  )

)



#####

#### Bootstrap by fitting mvgls on the simulated datasets ####
#### NOTE: uses the MCC tree and simulated data as inputs



### Apply the appropriate parallelization plan


if(exists("parallelize") == FALSE || is.null(parallelize) == TRUE){
  
  print("You haven't specified whether or not to parallelize!")
  
} else if (parallelize == FALSE) {
  
  cat(
    "Although you selected no parallelization in the mvSIM part above,"
    )
  cat(
    "this section does require parallelization, so it has been set up!"
    )
  
  # Total number of cores to use in parallel
  #n.cores <- 4 # laptops
  #n.cores <- 8 # powerful laptops/desktops
  #n.cores <- 16 # workstations/servers
  n.cores <- availableCores() - 
    ceiling(
      availableCores() * 0.05) # automatic using 95% available cores
  #n.cores <- N # manual, edit here
  
  # Set up the parallelization session using 'n.cores' number of cores
  plan(multisession, workers = n.cores)
  
} else if (parallelize == TRUE) {
  
  print("You have already set up parallelization and are ready to go!")
  
}



### Run the bootstraps


## EB model fit

if(exists("set.group") == FALSE || is.null(set.group) == TRUE){
  
  print(
  "You haven't let the script know what clade you're working with yet!"
  )
  
} else if (set.group == "holosteans") {
  
  system.time(
    boot.EB.pca <- future_lapply(1:n.sim, function(x){
      mvgls(sim.EB[[x]] ~ 1
            , tree = tree.unitheight
            , scale.height = F
            , model = "EB"
            , method = "PL-LOOCV" #default
            #, method = "LL"
            #, method = "H&L"
            #, method = "Mahalanobis"
            , REML = FALSE #FALSE = use ML to estimate parameter
            , target = "unitVariance"
            #, target = "Variance"
            #, target = "null" #only with penalty = "RidgeAlt"
            , error = T
            , penalty ="RidgeArch" #default
            #, penalty = "RidgeAlt"
            , upper = bounds[[1]][[2]]
            , lower = bounds[[1]][[1]]
      )
    },future.seed = TRUE)
  )
  
} else {
  
  system.time(
    boot.EB.pca <- future_lapply(1:n.sim, function(x){
      mvgls(sim.EB[[x]] ~ 1
            , tree = tree.unitheight
            , scale.height = F
            , model = "EB"
            , method = "PL-LOOCV" #default
            #, method = "LL"
            #, method = "H&L"
            #, method = "Mahalanobis"
            #, REML = FALSE #FALSE = use ML to estimate parameter
            , target = "unitVariance"
            #, target = "Variance"
            #, target = "null" #only with penalty = "RidgeAlt"
            , error = T
            , penalty ="RidgeArch" #default
            #, penalty = "RidgeAlt"
            , upper = bounds[[1]][[2]]
            , lower = bounds[[1]][[1]]
      )
    },future.seed = TRUE)
  )
  
}


## BM model fit

system.time(
boot.BM.pca <- future_lapply(1:n.sim, function(x){
  mvgls(sim.BM[[x]] ~ 1
        , tree = tree.unitheight
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
        )
},future.seed = TRUE)
)


## OU model fit

system.time(
boot.OU.pca <- future_lapply(1:n.sim, function(x){
  mvgls(sim.OU[[x]] ~ 1
        , tree = tree.unitheight
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
  )
},future.seed = TRUE)
)



#####

#### CHECKPOINT: Export the bootstrap results ####



### Store the outputs according to the group you're testing and 
### standardize


if (set.group == "coelacanths") {
  
  boot.pca.coelacanths <- NULL
  boot.pca.coelacanths[[1]] <- boot.EB.pca # 1 is always EB
  boot.pca.coelacanths[[2]] <- boot.BM.pca # 2 is always BM
  boot.pca.coelacanths[[3]] <- boot.OU.pca # 3 is always OU
  boots.pca <- boot.pca.coelacanths
  saveRDS(boot.pca.coelacanths
        ,'fitboot_coelacanth_pcamodeltests_100.RDS'
        )
  
} else if (set.group == "lungfishes") {
  
  boot.pca.lungfishes <- NULL
  boot.pca.lungfishes[[1]] <- boot.EB.pca
  boot.pca.lungfishes[[2]] <- boot.BM.pca
  boot.pca.lungfishes[[3]] <- boot.OU.pca
  boots.pca <- boot.pca.lungfishes
  saveRDS(boot.pca.lungfishes,'fitboot_lungfish_pcamodeltests_100.RDS')
  
} else if (set.group == "holosteans") {
  
  boot.pca.holosteans <- NULL
  boot.pca.holosteans[[1]] <- boot.EB.pca
  boot.pca.holosteans[[2]] <- boot.BM.pca
  boot.pca.holosteans[[3]] <- boot.OU.pca
  boots.pca <- boot.pca.holosteans
  saveRDS(boot.pca.holosteans,'fitboot_holostean_pcamodeltests_100.RDS')
  
} else {
  
  print("Something went wrong. Have you set your clade correctly?")
  
}



### Compare the fit to the estimated EB, BM, OU
### NOTE: we are storing all results for completeness, but are only
### assessing the results for the best fit model (EB in all cases)


## EB

# EB beta parameter
boot.EB.pca.beta <- sapply(boots.pca[[1]], function(x) x$param)

# Standard deviation of the multivariate Brownian motion process
# NOTE: This is a measure of the overall speed or rate at which the set
# of traits is changing over time due to Brownian motion.
boot.EB.pca.sdmv <- 
  sapply(boots.pca[[1]], function(x) sqrt(sum(diag(x$sigma$Pinv))))

# Square root of the determinant 
# NOTE: This represents the volume or "size" of the data cloud, which 
# can also be interpreted as a measure of the overall rate of diffusion
# for the multivariate process.
boot.EB.pca.volmv <- 
  sapply(boots.pca[[1]], function(x) sqrt(det(x$sigma$Pinv)))


## BM

# BM doesn't have a beta nor alpha parameters
#NA

# Standard deviation of the multivariate Brownian motion process
boot.BM.pca.sdmv <- 
  sapply(boots.pca[[2]], function(x) sqrt(sum(diag(x$sigma$Pinv))))

# Square root of the determinant 
boot.BM.pca.volmv <- 
  sapply(boots.pca[[2]], function(x) sqrt(det(x$sigma$Pinv)))


## OU

# OU alpha parameter
boot.OU.pca.alpha <- sapply(boots.pca[[3]], function(x) x$param)

# Standard deviation of the multivariate Brownian motion process
boot.OU.pca.sdmv <- 
  sapply(boots.pca[[3]], function(x) sqrt(sum(diag(x$sigma$Pinv))))

# Square root of the determinant 
boot.OU.pca.volmv <- 
  sapply(boots.pca[[3]], function(x) sqrt(det(x$sigma$Pinv)))



### Store the outputs according to the group you're testing and 
### standardize


if (set.group == "coelacanths") {
  
  boot.pca.coelacanths.stats <- NULL
  boot.pca.coelacanths.stats[[1]][[1]] <- boot.EB.pca.beta # 1 is beta
  boot.pca.coelacanths.stats[[1]][[2]] <- boot.EB.pca.sdmv # 2 is sdmv
  boot.pca.coelacanths.stats[[1]][[3]] <- boot.EB.pca.volmv # 3 is volmv
  boot.pca.coelacanths.stats[[2]] <- list() # initiate with empty list
  boot.pca.coelacanths.stats[[2]][[1]] <- list() # no beta
  boot.pca.coelacanths.stats[[2]][[2]] <- boot.BM.pca.sdmv
  boot.pca.coelacanths.stats[[2]][[3]] <- boot.BM.pca.volmv
  boot.pca.coelacanths.stats[[3]] <- list() # initiate with empty list
  boot.pca.coelacanths.stats[[3]][[1]] <- NA # no beta
  boot.pca.coelacanths.stats[[3]][[2]] <- boot.OU.pca.sdmv
  boot.pca.coelacanths.stats[[3]][[3]] <- boot.OU.pca.volmv
  boot.pca.coelacanths.stats[[3]][[4]] <- boot.OU.pca.alpha # 4 is alpha
  boots.pca.stats <- boot.pca.coelacanths.stats
  saveRDS(boots.pca.stats
          , "fitboot_coelacanth_pcamodeltests_stats_100.RDS"
          )
  
} else if (set.group == "lungfishes") {
  
  boot.pca.lungfishes.stats <- NULL
  boot.pca.lungfishes.stats[[1]][[1]] <- boot.EB.pca.beta # 1 is beta
  boot.pca.lungfishes.stats[[1]][[2]] <- boot.EB.pca.sdmv # 2 is sdmv
  boot.pca.lungfishes.stats[[1]][[3]] <- boot.EB.pca.volmv # 3 is volmv
  boot.pca.lungfishes.stats[[2]] <- list() # initiate with empty list
  boot.pca.lungfishes.stats[[2]][[1]] <- list() # no beta
  boot.pca.lungfishes.stats[[2]][[2]] <- boot.BM.pca.sdmv
  boot.pca.lungfishes.stats[[2]][[3]] <- boot.BM.pca.volmv
  boot.pca.lungfishes.stats[[3]] <- list() # initiate with empty list
  boot.pca.lungfishes.stats[[3]][[1]] <- list() # no beta
  boot.pca.lungfishes.stats[[3]][[2]] <- boot.OU.pca.sdmv
  boot.pca.lungfishes.stats[[3]][[3]] <- boot.OU.pca.volmv
  boot.pca.lungfishes.stats[[3]][[4]] <- boot.OU.pca.alpha # 4 is alpha
  boots.pca.stats <- boot.pca.lungfishes.stats
  saveRDS(boots.pca.stats
          , "fitboot_lungfish_pcamodeltests_stats_100.RDS"
          )
  
} else if (set.group == "holosteans") {
  
  boot.pca.holosteans.stats <- NULL
  boot.pca.holosteans.stats[[1]][[1]] <- boot.EB.pca.beta # 1 is beta
  boot.pca.holosteans.stats[[1]][[2]] <- boot.EB.pca.sdmv # 2 is sdmv
  boot.pca.holosteans.stats[[1]][[3]] <- boot.EB.pca.volmv # 3 is volmv
  boot.pca.holosteans.stats[[2]] <- list() # initiate with empty list
  boot.pca.holosteans.stats[[2]][[1]] <- list() # no beta
  boot.pca.holosteans.stats[[2]][[2]] <- boot.BM.pca.sdmv
  boot.pca.holosteans.stats[[2]][[3]] <- boot.BM.pca.volmv
  boot.pca.holosteans.stats[[3]] <- list() # initiate with empty list
  boot.pca.holosteans.stats[[3]][[1]] <- list() # no beta
  boot.pca.holosteans.stats[[3]][[2]] <- boot.OU.pca.sdmv
  boot.pca.holosteans.stats[[3]][[3]] <- boot.OU.pca.volmv
  boot.pca.holosteans.stats[[3]][[4]] <- boot.OU.pca.alpha # 4 is alpha
  boots.pca.stats <- boot.pca.holosteans.stats
  saveRDS(boots.pca.stats
          , "fitboot_holostean_pcamodeltests_stats_100.RDS"
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
    
    boots.pca <- readRDS('fitboot_coelacanth_pcamodeltests_100.RDS')

    boots.pca.stats <- 
      readRDS("fitboot_coelacanth_pcamodeltests_stats_100.RDS")
    
  } else if (set.group == "lungfishes") {
    
    boots.pca <- readRDS('fitboot_lungfish_pcamodeltests_100.RDS')

    boots.pca.stats <- 
      readRDS("fitboot_lungfish_pcamodeltests_stats_100.RDS")
    
  } else if (set.group == "holosteans") {
    
    boots.pca <- readRDS('fitboot_holostean_pcamodeltests_100.RDS')

    boots.pca.stats <- 
      readRDS("fitboot_holostean_pcamodeltests_stats_100.RDS")
    
  }
  
}



#####

#### Estimate the 95% confidence interval for the bootstraps ####


if (set.group == "coelacanths") {
  
  print("EB beta parameter")
  print(fits.pca[[1]]$param[[1]]) # EB (1)
  
  cat("\n")
  
  print("EB beta parameter (95% confidence interval)")
  print(hdi(boots.pca.stats[[1]][[1]])) # EB (1) beta (1)
  
  cat("\n")
  
} else if (set.group == "lungfishes") {
  
  print("EB beta parameter")
  print(fits.pca[[1]]$param[[1]]) # EB (1)
  
  cat("\n")
  
  print("EB beta parameter (95% confidence interval)")
  print(hdi(boots.pca.stats[[1]][[1]])) # EB (1) beta (1)
  
  cat("\n")
  
} else if (set.group == "holosteans") {
  
  print("EB beta parameter")
  print(fits.pca[[1]]$param[[1]]) # EB (1)
  
  cat("\n")
  
  print("EB beta parameter (95% confidence interval)")
  print(hdi(boots.pca.stats[[1]][[1]])) # EB (1) beta (1)
  
  cat("\n")
  
} else {
  
  print("Something went wrong. Have you set your clade correctly?")
  
}



#####

#### Plot the bootstrap results ####



### Identify the clade being worked with
### NOTE: as long as you backed up all analyses for all groups using the 
### checkpoint above, you can use this section to reload whichever 
### dataset you want from the R environment. Otherwise, it won't affect
### anything.


set.group <- "coelacanths"
#set.group <- "lungfishes"
#set.group <- "holosteans"


if (set.group == "coelacanths") {
  
  fits.pca <- fit.pca.coelacanths
  boots.pca <- boot.pca.coelacanths
  boots.pca.stats <- boot.pca.coelacanths.stats

} else if (set.group == "lungfishes") {
  
  fits.pca <- fit.pca.lungfishes
  boots.pca <- boot.pca.lungfishes
  boots.pca.stats <- boot.pca.lungfishes.stats

} else if (set.group == "holosteans") {
  
  fits.pca <- fit.pca.holosteans
  boots.pca <- boot.pca.holosteans
  boots.pca.stats <- boot.pca.holosteans.stats

} else {
  
  print("Something went wrong. Have you set your clade correctly?")
  
}



### Set up to automatically generate the plot titles
### NOTE: this uses the group you have as your 'set.group' and inserts 
### it into the title from that object.

plot.title <- paste("bootstrap samples (",set.group,")", sep = "")



### Generate a plot of bootstrap parameter estimates
### NOTE: the histogram will be the bootstraps and the abline will be
### the empirical estimate from the MCC tree.

hist(
  boots.pca.stats[[1]][[1]]
  , xlim = c(-30,0)
  , ylim = c(0,20)
  , breaks = 20
  , main = paste("EB", plot.title)
  , xlab = c("beta")
  )

abline(
  v = fits.pca[[1]]$param
  , lwd = 2
  , col = "red"
  )


#####
