#
##
#### PROCESSING OF RAW DATA 
#### FOR EVOLUTIONARY RATE, SATURATION, AND SHAPE ANALYSES
##
#

#### Script by Rafael A. Rivero Vega
#### Morphological rate visualization based off original in 
#### Close et al 2016 (https://doi.org/10.1111/pala.12245), with
#### further details in the manuscript methods and acknowledgements.



#### Load the required libraries ####


if(!'phytools'%in%installed.packages()
  ) install.packages("phytools"); library(phytools)

if(!'ape'%in%installed.packages()
  ) install.packages("ape"); library(ape)

if(!'phylotate'%in%installed.packages()
  ) install.packages("phylotate"); library(phylotate)

if(!'coda'%in%installed.packages()
  ) install.packages("coda"); library(coda)

if(!'geomorph'%in%installed.packages()
  ) install.packages("geomorph"); library(geomorph)

if(!'StereoMorph'%in%installed.packages()
  ) install.packages("StereoMorph"); library(StereoMorph)

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
### NOTE: read.annotated.nexus is from a deprecated package 
### called OutbreakTools, originally published in Jombart et al 2014
### (https://doi.org/10.1016/j.epidem.2014.04.003)

source("./readannotatednexus.R")
#source(paste(script1))



#####

#### Import the raw data to process the trees ####



### Identify the clade being worked with
### NOTE: this setting allows for the import of all the raw data, 
### preps objects, and adjusts the way outputs are saved and backed
### up further down


set.group <- "coelacanths"
#set.group <- "lungfishes"
#set.group <- "holosteans"



### Import Raw Data
### NOTE: this will use the read.annotated.nexus function to import the
### mcc tree, the posterior trees file, and a list of outgroups to prune


if(exists("set.group") == FALSE || is.null(set.group) == TRUE){
  
  print(
  "You haven't let the script know what clade you're working with yet!"
  )
  
} else if (set.group == "coelacanths") {
  
  # Maximum Clade Credibility tree (TreeAnnotator output)
  tree <- 
    read.annotated.nexus(
      "./beast_coelacanths_morphology_Torino-etal_tree-mcc.tre"
      )
  
  # Posterior distribution of trees (BEAST output)
  # NOTE: this may take a bit depending on the filesize and your 
  # computer
  trees <- 
    read.annotated.nexus(
      "./beast_coelacanths_morphology_Torino-etal.trees"
      )
  
  # List of outgroup taxa to prune from the tree
  # NOTE: should be named *exactly* the same as the tips on the trees
  # above
  tip.drop <- 
    readLines("tree_taxa_coelacanth_outgroups.txt")
  
  
} else if (set.group == "lungfishes") {
  
  # Maximum Clade Credibility tree (TreeAnnotator output)
  tree <- 
    read.annotated.nexus(
      "./beast_lungfishes_morphology_Lloyd-etal_tree-mcc.tre"
      )
  
  # Posterior distribution of trees (BEAST output)
  # NOTE: this may take a bit depending on the filesize and your 
  # computer
  trees <- 
    read.annotated.nexus(
      "./beast_lungfishes_morphology_Lloyd-etal.trees"
      )
  
  # List of outgroup taxa to prune from the tree
  # NOTE: should be named *exactly* the same as the tips on the trees
  # above
  tip.drop <- 
    readLines("tree_taxa_lungfish_outgroups.txt")
  
  
} else if (set.group == "holosteans") {
  
  # Maximum Clade Credibility tree (TreeAnnotator output)
  tree <- 
    read.annotated.nexus(
    "./beast_neopterygii_morphology_LopezArbarello-Sferco_tree-mcc.tre"
    )
  
  # Posterior distribution of trees (BEAST output)
  # NOTE: this may take a bit depending on the filesize and your 
  # computer
  trees <- 
    read.annotated.nexus(
      "./beast_neopterygii_morphology_LopezArbarello-Sferco.trees"
      )
  
  # List of outgroup taxa to prune from the tree
  # NOTE: should be named *exactly* the same as the tips on the trees
  # above
  tip.drop <- 
    readLines("tree_taxa_holostean_outgroups.txt")
  
  
}



#####

#### Prune outgroups from all trees and standardize ####



### Drop tips on the MCC tree while preserving the imported annotations


## Preserve edge lengths with their ascending numbers to use as an
## index further down
## NOTE: it is crucial that this is done before dropping tips!

edge.lengths <- tree$edge.length


## Bind the edge lengths together with their numerical values
## NOTE: it is crucial that you preserve the numbers in ascending order
## before joining edge lengths here! 

edge.lengths <- cbind(1:length(edge.lengths),edge.lengths)



### Create a final tree object with the adjusted annotations


## Exclude tips using outgroup data

# Drop outgroup tips that may skew rates
tree.drop <- drop.tip(tree,tip.drop)

# Create a post-drop edge length data frame
edge.lengths.drop <- tree.drop$edge.length


## Unite the adjusted tree and annotations 

# Isolate the edge values that are in the first index but not the 
# second (aka the edges that have been dropped)
edge.lengths.drop <- 
  edge.lengths[!(edge.lengths[,2] %in% edge.lengths.drop),]

# Adjust the index of dropped edges so it goes in descending order
# NOTE: this order is crucial because of the way list items are 
# edited (see loop below)
edge.lengths.drop <- 
  edge.lengths.drop[order(edge.lengths.drop[,1], decreasing = TRUE),]

# Extract the full list of annotations from the tree file
# NOTE: even though we've dropped tips, we need to drop the associated
# annotations too, which is what this whole section is supposed to do
annotations <- tree.drop$annotations

# Loop through each edge number in descending order and remove it from
# the annotations object
# NOTE: because the index is in descending order, this loop will move
# from the *bottom* of the list to the top, deleting the relevant
# entries as it goes. It *needs* to be done this way or else the
# indexing doesn't work because the list numbers are adjusted to fill
# the gap every time any part of the list is removed.
for (i in edge.lengths.drop[,1]){
  annotations[[i]] <- NULL
}

# Replace the full list of annotations with the updated list we just
# created as the annotations object
tree.drop$annotations <- annotations

# Extract the new root time (post tip drops) and add it as an extra
# annotation for later use
tree.drop$root.time <- max(diag(vcv(tree.drop)))



### Drop tips on the posterior trees while preserving the imported
### annotations
### NOTE: this is essentially repeating what was already done for the
### MCC tree above but across all posterior trees


## Prepare for the process

# Create an empty object to hold the list of trees with tips dropped
trees.drop <- NULL

## Process posterior trees

# Loop through all posterior trees and apply the same method to drop
# tips as for the MCC
for (i in 1:length(trees)){
  
  # Preserve edge lengths with their ascending numbers to use as an
  # index further down
  # NOTE: it is crucial that this is done before dropping tips!
  edge.lengths <- trees[[i]]$edge.length

  # Bind the edge lengths together with their numerical values
  # NOTE: it is crucial that you preserve the numbers in ascending 
  # order before joining edge lengths here! 
  edge.lengths <- cbind(1:length(edge.lengths),edge.lengths)
  
  # Drop outgroup tips that may skew rates
  trees.drop[[i]] <- drop.tip(trees[[i]],tip.drop)
  
  # Create a post-drop edge length data frame
  edge.lengths.drop <- trees.drop[[i]]$edge.length
  
  # Isolate the edge values that are in the first index but not the 
  # second (aka the edges that have been dropped)
  edge.lengths.drop <- 
    edge.lengths[!(edge.lengths[,2] %in% edge.lengths.drop),]
  
  # Adjust the index of dropped edges so it goes in descending order
  # of edge number
  # NOTE: this order is crucial because of the way list items are
  # edited (see loop below)
  edge.lengths.drop <- 
    edge.lengths.drop[order(edge.lengths.drop[,1], decreasing = TRUE),]
  
  # Extract the full list of annotations from the tree file
  # NOTE: even though we've dropped tips, we need to drop the 
  # associated annotations too, which is what this whole section is
  # supposed to do
  annotations <- trees.drop[[i]]$annotations
  
  # Loop through each edge number in descending order and remove it 
  # from the annotations object
  # NOTE: because the index is in descending order, this loop will move
  # from the *bottom* of the list to the top, deleting the relevant
  # entries as it goes. It *needs* to be done this way or else the
  # indexing doesn't work because the list numbers are adjusted to fill
  # the gap every time any part of the list is removed.
  for (j in edge.lengths.drop[,1]){
    annotations[[j]] <- NULL
  }
  
  # Replace the full list of annotations with the updated list we just
  # created as the annotations object
  trees.drop[[i]]$annotations <- annotations
  
  # Extract the new root time (post tip drops) and add it as an extra
  # annotation for later use
  trees.drop[[i]]$root.time <- max(diag(vcv(trees.drop[[i]])))
}



### Standardize tip labels
### NOTE: this only has to be done if your tree$tip.label contains
### extra info beyond Genus_species (e.g., in BEAST where you can
### include ages appended to tip labels). In our case, that means
### this script is written to only apply this to lungfishes and
### holosteans, but not to coelacanths


if(exists("set.group") == FALSE || is.null(set.group) == TRUE){
  
  print(
  "You haven't let the script know what clade you're working with yet!"
  )
  
} else if (set.group == "coelacanths") {
  
  # Save the MCC tree to its own object
  tree.mcc.coelacanth.morph <- tree.drop
  
  # Save the posterior trees to their own object
  trees.posterior.coelacanth.morph <- trees.drop

  # Ensure that the posterior tree object is of class multiPhylo so it
  # can be read and processed without issue
  class(trees.posterior.coelacanth.morph) <- "multiPhylo"
  
} else if (set.group == "lungfishes") {
  
  # Standardize tips to 'Genus_species' in the MCC tree

  #"F_" fossil identifier
  tree.drop$tip.label <- 
    gsub("F_","",tree.drop$tip.label)
  
  #Ages in the name
  tree.drop$tip.label <- 
    gsub("(.*)_[0-9].*","\\1",tree.drop$tip.label)
  
  # Standardize tips to 'Genus_species' in the posterior trees

  #Loop to apply across all posterior trees
  for (i in 1:length(trees.drop)){
    
    #"F_" fossil identifier
    trees.drop[[i]]$tip.label <- 
      gsub("F_","",trees.drop[[i]]$tip.label)
    
    #Ages at the end
    trees.drop[[i]]$tip.label <- 
      gsub("(.*)_[0-9].*","\\1",trees.drop[[i]]$tip.label)
    
  }
  
  # Save the MCC tree to its own object
  tree.mcc.lungfish.morph <- tree.drop
  
  # Save the posterior trees to their own object
  trees.posterior.lungfish.morph <- trees.drop

  # Ensure this object is also of class multiPhylo
  class(trees.posterior.lungfish.morph) <- "multiPhylo"
  
} else if (set.group == "holosteans") {
  
  # Standardize tips to 'Genus_species' in the MCC tree
  
  #"F_" fossil identifier
  tree.drop$tip.label <- 
    gsub("F_","",tree.drop$tip.label)

  #Ages in the name
  tree.drop$tip.label <- 
    gsub("(.*)_[0-9].*","\\1",tree.drop$tip.label)
  
  # Standardize tips to 'Genus_species' in the posterior trees
  
  #Loop to apply across all posterior trees
  for (i in 1:length(trees.drop)){
    
    #"F_" fossil identifier
    trees.drop[[i]]$tip.label <- 
      gsub("F_","",trees.drop[[i]]$tip.label)
    
    #Ages at the end
    trees.drop[[i]]$tip.label <- 
      gsub("(.*)_[0-9].*","\\1",trees.drop[[i]]$tip.label)
  }
  
  # Save the MCC tree to its own object
  tree.mcc.holostean.morph <- tree.drop
  
  # Save the posterior trees to their own object
  trees.posterior.holostean.morph <- trees.drop
  
  # Ensure this object is also of class multiPhylo
  class(trees.posterior.holostean.morph) <- "multiPhylo"
  
}


#####


#### CHECKPOINT: Export pruned and standardized trees ####



### Store the outputs according to the group you're testing and
### standardize


if (set.group == "coelacanths") {
  
  # Back up the pruned and standardized MCC tree 
  saveRDS(tree.mcc.coelacanth.morph
          ,"tree_mcc_coelacanth_pruned.RDS"
          )
  
  # Back up the pruned and standardized posterior trees
  saveRDS(trees.posterior.coelacanth.morph
          ,"trees_posterior_coelacanth_pruned.RDS"
          )
  
} else if (set.group == "lungfishes") {
  
  # Back up the pruned and standardized MCC tree 
  saveRDS(tree.mcc.lungfish.morph
          ,"tree_mcc_lungfish_pruned.RDS"
          )
  
  # Back up the pruned and standardized posterior trees
  saveRDS(trees.posterior.lungfish.morph
          ,"trees_posterior_lungfish_pruned.RDS"
          )
  
} else if (set.group == "holosteans") {
  
  # Back up the pruned and standardized MCC tree 
  saveRDS(tree.mcc.holostean.morph
          ,"tree_mcc_holostean_pruned.RDS"
          )
  
  # Back up the pruned and standardized posterior trees
  saveRDS(trees.posterior.holostean.morph
          ,"trees_posterior_holostean_pruned.RDS"
          )
  
} else {
  
  print(
    "Something went wrong. Have you set the clade you're working with?"
    )
  
}



### Do you want to read in the checkpoint files?
### NOTE: use this if you don't want to reload the raw data

#load.checkpoint <- TRUE
load.checkpoint <- FALSE


if(
  exists("load.checkpoint") == FALSE || is.null(load.checkpoint) == TRUE
  ){
  
  print(
  "You haven't set if you want to load from this checkpoint!"
  )
  
} else if (load.checkpoint == FALSE) {
  
  print("You haven't opted to load any saved files.")
  
} else if (load.checkpoint == TRUE) {
  
  if (set.group == "coelacanths") {
    
    tree.mcc.coelacanth.morph <- 
      readRDS("tree_mcc_coelacanth_pruned.RDS")
    
    trees.posterior.coelacanth.morph <- 
      readRDS("trees_posterior_coelacanth_pruned.RDS")
    
  } else if (set.group == "lungfishes") {
    
    tree.mcc.lungfish.morph <- 
      readRDS("tree_mcc_lungfish_pruned.RDS")
    
    trees.posterior.lungfish.morph <- 
      readRDS("trees_posterior_lungfish_pruned.RDS")
    
  } else if (set.group == "holosteans") {
    
    tree.mcc.holostean.morph <- 
      readRDS("tree_mcc_holostean_pruned.RDS")
    
    trees.posterior.holostean.morph <- 
      readRDS("trees_posterior_holostean_pruned.RDS")
    
  }
  
}



#####


#### Extract data from the trees ####



### Extract data from MCC trees

if(exists("set.group") == FALSE || is.null(set.group) == TRUE){
  
  print(
  "You haven't let the script know what clade you're working with yet!"
  )
  
} else if (set.group == "coelacanths") {
  
  # Create a new object
  tree.coelacanths <- tree.mcc.coelacanth.morph
  
  # Extract rates
  coelacanth.rates.morph.mcc <- 
    as.vector(unlist(sapply(tree.coelacanths$annotations
                            , FUN = function(x) print(100*x$rate_median)
                            )))
  
  coelacanth.heights.morph.mcc <- 
    as.vector(unlist(sapply(tree.coelacanths$annotations
                            , FUN = function(x) print(x$height_median)
                            )))
  
  coelacanth.heights.morph.mcc <- 
    round(coelacanth.heights.morph.mcc, 3)
  
  # Get row numbers of MCC terminal branches
  coelacanth.branches.terminal.mcc <- 
    match(1:Ntip(tree.coelacanths), tree.coelacanths$edge[,2])
  
  # Get row numbers of MCC internal branches
  coelacanth.branches.internal.mcc <- 
    setdiff(1:nrow(tree.coelacanths$edge)
      , coelacanth.branches.terminal.mcc
      )
  
  # Get MCC node heights
  coelacanth.node.heights.mcc <- nodeHeights(tree.coelacanths)
  
  # Get MCC node depths
  coelacanth.node.depths.mcc <- 
    zapsmall(tree.coelacanths$root.time - coelacanth.node.heights.mcc)
  
  # Get MCC age midpoints
  coelacanth.midpoint.ages.mcc <- 
    (
    (coelacanth.node.depths.mcc[,1] + coelacanth.node.depths.mcc[,2])/2
    )
  
  
} else if (set.group == "lungfishes") {
  
  # Create a new object
  tree.lungfishes <- tree.mcc.lungfish.morph
  
  # Extract rates
  lungfish.rates.morph.mcc <- 
    as.vector(unlist(sapply(tree.lungfishes$annotations
                            , FUN = function(x) print(100*x$rate_median)
    )))
  
  lungfish.heights.morph.mcc <- 
    as.vector(unlist(sapply(tree.lungfishes$annotations
                            , FUN = function(x) print(x$height_median)
    )))
  
  lungfish.heights.morph.mcc <- 
    round(lungfish.heights.morph.mcc, 3)
  
  # Get row numbers of MCC terminal branches
  lungfish.branches.terminal.mcc <- 
    match(1:Ntip(tree.lungfishes), tree.lungfishes$edge[,2])
  
  # Get row numbers of MCC internal branches
  lungfish.branches.internal.mcc <- 
    setdiff(1:nrow(tree.lungfishes$edge)
      , lungfish.branches.terminal.mcc
      )
  
  # Get MCC node heights
  lungfish.node.heights.mcc <- nodeHeights(tree.lungfishes)
  
  # Get MCC node depths
  lungfish.node.depths.mcc <- 
    zapsmall(tree.lungfishes$root.time - lungfish.node.heights.mcc)
  
  # Get MCC age midpoints
  lungfish.midpoint.ages.mcc <- 
    ((lungfish.node.depths.mcc[,1] + lungfish.node.depths.mcc[,2])/2)
  
} else if (set.group == "holosteans") {
  
  # Create a new object
  tree.holosteans <- tree.mcc.holostean.morph
  
  # Extract rates
  holostean.rates.morph.mcc <- 
    as.vector(unlist(sapply(tree.holosteans$annotations
                            , FUN = function(x) print(100*x$rate_median)
    )))
  
  holostean.heights.morph.mcc <- 
    as.vector(unlist(sapply(tree.holosteans$annotations
                            , FUN = function(x) print(x$height_median)
    )))
  
  holostean.heights.morph.mcc <- 
    round(holostean.heights.morph.mcc, 3)
  
  # Get row numbers of MCC terminal branches
  holostean.branches.terminal.mcc <- 
    match(1:Ntip(tree.holosteans), tree.holosteans$edge[,2])
  
  # Get row numbers of MCC internal branches
  holostean.branches.internal.mcc <- 
    setdiff(1:nrow(tree.holosteans$edge)
      , holostean.branches.terminal.mcc
      )
  
  # Get MCC node heights
  holostean.node.heights.mcc <- nodeHeights(tree.holosteans)
  
  # Get MCC node depths
  holostean.node.depths.mcc <- 
    zapsmall(tree.holosteans$root.time - holostean.node.heights.mcc)
  
  # Get MCC age midpoints
  holostean.midpoint.ages.mcc <- 
    ((holostean.node.depths.mcc[,1] + holostean.node.depths.mcc[,2])/2)
  
}



### Extract data from a sample of the posterior trees


## Amount of trees to sample
## NOTE: this sample will be taken post burn-in, so the first 10% of
## trees are ignored in the calculations that follow

nsample <- 100


## Check Effective Sample Size (ESS) values

# Import the appropriate log file

if(exists("set.group") == FALSE || is.null(set.group) == TRUE){
  
  print(
  "You haven't let the script know what clade you're working with yet!"
  )
  
  
} else if (set.group == "coelacanths") {
  
  data <- read.table(
    file = "beast_coelacanths_morphology_Torino-etal_log.log"
    , header = T
  )
  
  
} else if (set.group == "lungfishes") {
  
  data <- read.table(
    file="beast_lungfishes_morphology_Lloyd-etal_log.log"
    , header = T
  )
  
  
} else if (set.group == "holosteans") {
  
  data <- read.table(
    file="beast_neopterygii_morphology_LopezArbarello-Sferco_log.log"
    , header = T
  )
  
  
}

# Plot the trace

plot(tail(data$likelihood
          , length(data$likelihood) - (0.1*length(data$likelihood))
          )
     , type='b'
     , cex=0.01
     )            

# Create an MCMC object

c.mcmc <- as.mcmc(tail(data$likelihood
                       , length(data$likelihood
                        ) - (0.1*length(data$likelihood))
                       )
                  )

# Estimate the ESS

esses <- list()

for(i in 1:nsample) {
  esses[i] <- effectiveSize(sample(c.mcmc, nsample))
}

# Find the mean ESS across

mean(unlist(esses)) # all should be >100 in our case


## Burn in, sample, and extract data

if(exists("set.group") == FALSE || is.null(set.group) == TRUE){
  
  print(
  "You haven't let the script know what clade you're working with yet!"
  )
  
  
} else if (set.group == "coelacanths") {
  
  # Burn in
  burnin.coelacanth <- 
    round(length(trees.posterior.coelacanth.morph)*0.1, digits = -1)
  
  # Sample
  trees.posterior.coelacanth.morph.sample <- 
    sample(
      trees.posterior.coelacanth.morph[
        burnin.coelacanth:length(trees.posterior.coelacanth.morph)
      ]
      , nsample
    )
  
  class(trees.posterior.coelacanth.morph.sample) <- "multiPhylo"
  
  # Prepare empty list objects
  coelacanth.rates.morph.posterior <- list()
  coelacanth.node.heights.posterior <- list()
  coelacanth.node.depths.posterior <- list()
  coelacanth.midpoint.ages.posterior <- list()
  
  # Loop through all posterior trees, extract individual data, and log 
  # scale
  for (i in 1:length(trees.posterior.coelacanth.morph.sample)) {
    
    trees.posterior.coelacanth.morph.sample[[i]]$root.time <- 
      max(diag(vcv(trees.posterior.coelacanth.morph.sample[[i]])))
    
    coelacanth.rates.morph.posterior[[i]] <- 
      as.vector(
        unlist(
          sapply(
            trees.posterior.coelacanth.morph.sample[[i]]$annotations
            , FUN = function(x) print(100*x$rate)
            )
          )
        )
    
    coelacanth.node.heights.posterior[[i]] <- 
      nodeHeights(trees.posterior.coelacanth.morph.sample[[i]])
    
    coelacanth.node.depths.posterior[[i]] <- 
      zapsmall(
        trees.posterior.coelacanth.morph.sample[[i]]$root.time - 
        coelacanth.node.heights.posterior[[i]]
        )
    
    coelacanth.midpoint.ages.posterior[[i]] <- 
      ((coelacanth.node.depths.posterior[[i]][,1] + 
        coelacanth.node.depths.posterior[[i]][,2])/2
      )
  
    }
  
  
} else if (set.group == "lungfishes") {
  
  # Burn in
  burnin.lungfish <- 
    round(length(trees.posterior.lungfish.morph)*0.1, digits = -1)
  
  # Sample
  trees.posterior.lungfish.morph.sample <- 
    sample(
      trees.posterior.lungfish.morph[
        burnin.lungfish:length(trees.posterior.lungfish.morph)
      ]
      , nsample
    )
  
  class(trees.posterior.lungfish.morph.sample) <- "multiPhylo"
  
  # Prepare empty list objects
  lungfish.rates.morph.posterior <- list()
  lungfish.node.heights.posterior <- list()
  lungfish.node.depths.posterior <- list()
  lungfish.midpoint.ages.posterior <- list()
  
  # Loop through all posterior trees, extract individual data, and log 
  # scale
  for (i in 1:length(trees.posterior.lungfish.morph.sample)) {
    
    trees.posterior.lungfish.morph.sample[[i]]$root.time <- 
      max(diag(vcv(trees.posterior.lungfish.morph.sample[[i]])))
    
    lungfish.rates.morph.posterior[[i]] <- 
      as.vector(
        unlist(
          sapply(
            trees.posterior.lungfish.morph.sample[[i]]$annotations
            , FUN = function(x) print(100*x$rate)
            )
          )
        )
    
    lungfish.node.heights.posterior[[i]] <- 
      nodeHeights(trees.posterior.lungfish.morph.sample[[i]])
    
    lungfish.node.depths.posterior[[i]] <- 
      zapsmall(
        trees.posterior.lungfish.morph.sample[[i]]$root.time - 
        lungfish.node.heights.posterior[[i]]
        )
    
    lungfish.midpoint.ages.posterior[[i]] <- 
      ((lungfish.node.depths.posterior[[i]][,1] + 
        lungfish.node.depths.posterior[[i]][,2])/2
      )
    
  }
  
  
} else if (set.group == "holosteans") {
  
  # Burn in
  burnin.holostean <- 
    round(length(trees.posterior.holostean.morph)*0.1, digits = -1)
  
  # Sample
  trees.posterior.holostean.morph.sample <- 
    sample(
      trees.posterior.holostean.morph[
        burnin.holostean:length(trees.posterior.holostean.morph)
      ]
      , nsample
    )
  
  class(trees.posterior.holostean.morph.sample) <- "multiPhylo"
  
  # Prepare empty list objects
  holostean.rates.morph.posterior <- list()
  holostean.node.heights.posterior <- list()
  holostean.node.depths.posterior <- list()
  holostean.midpoint.ages.posterior <- list()
  
  # Loop through all posterior trees, extract individual data, and log 
  # scale
  for (i in 1:length(trees.posterior.holostean.morph.sample)) {
    
    trees.posterior.holostean.morph.sample[[i]]$root.time <- 
      max(diag(vcv(trees.posterior.holostean.morph.sample[[i]])))
    
    holostean.rates.morph.posterior[[i]] <- 
      as.vector(
        unlist(
          sapply(
            trees.posterior.holostean.morph.sample[[i]]$annotations
            , FUN = function(x) print(100*x$rate)
            )
          )
        )
    
    holostean.node.heights.posterior[[i]] <- 
      nodeHeights(trees.posterior.holostean.morph.sample[[i]])
    
    holostean.node.depths.posterior[[i]] <- 
      zapsmall(trees.posterior.holostean.morph.sample[[i]]$root.time - 
        holostean.node.heights.posterior[[i]]
        )
    
    holostean.midpoint.ages.posterior[[i]] <- 
      ((holostean.node.depths.posterior[[i]][,1] + 
        holostean.node.depths.posterior[[i]][,2])/2
      )
    
  }
  
  
}





#####


#### Plot discrete rate data ####
#### NOTE: these plots are manually configurable, so they're not
#### automated and should be adjusted as needed



### Set plot properties


## Adjust plot settings

# Alphas

setting.alpha.branch.posterior = 0.15
setting.alpha.branch.mcc = 0.6

# Line widths

setting.lwd.branch.mcc = 2
setting.lwd.branch.posterior = 0.2
setting.lwd.loess = 6

# Colors

setting.color = "black"

# Axis ranges (manual)

#setting.range.xlim <- c(425,0)
#setting.range.ylim <- c(-7,2)



### Plot rate through time with a 'rate cloud'


## Coelacanths

# Axis ranges

setting.range.xlim <- 
  c(round(max(unlist(coelacanth.node.depths.mcc)), digits = -1), 0)

setting.range.ylim <- 
  ceiling(range(log(unlist(coelacanth.rates.morph.posterior))))

# Output info

setwd(wd.lf)

#pdf(
#  paste("./Coelacanths_MorphRates_BEAST_LineSegmentPlot-Posterior.pdf"
#  , sep = ""), height = 8, width = 12)

png(
  filename = 
  "./Coelacanths_MorphRates_BEAST_LineSegmentPlot-Posterior.png"
  , width = 1280
  , height = 720
  )

# Plot area

par(mar=c(5,4,4,5)+.1)

# Create plot layer

plot(0
     ,type = "n"
     ,xlim = setting.range.xlim
     #,ylim = setting.range.ylim
     ,ylim = c(-7,2)
     ,yaxt = "n"
     ,xaxs = "i"
     ,yaxs = "i"
     ,axes = FALSE
     ,xlab=""
     ,ylab=""
)

# Add posterior branch line segments
for (i in 1:length(coelacanth.rates.morph.posterior)) {
  segments(x0 = coelacanth.node.depths.posterior[[i]][,2]
           ,y0 = log(coelacanth.rates.morph.posterior[[i]])
           ,x1 = coelacanth.node.depths.posterior[[i]][,1]
           ,y1 = log(coelacanth.rates.morph.posterior[[i]])
           ,lwd = setting.lwd.branch.posterior
           ,col = adjustcolor(setting.color
            , alpha.f = setting.alpha.branch.posterior
            )
  )  
}

# Create new plot layer

par(new=TRUE)

plot(0
     ,type = "n"
     ,xlim = setting.range.xlim
     #,ylim = setting.range.ylim
     ,ylim = c(-7,2)
     ,xaxs = "i"
     ,yaxs = "i"
     ,axes = FALSE
     ,xlab=""
     ,ylab=""
)

# Add mcc branch line segments

segments(x0 = coelacanth.node.depths.mcc[,2]
         ,y0 = log(coelacanth.rates.morph.mcc)
         ,x1 = coelacanth.node.depths.mcc[,1]
         ,y1 = log(coelacanth.rates.morph.mcc)
         ,lwd = setting.lwd.branch.mcc
         ,col = adjustcolor(setting.color
          , alpha.f = setting.alpha.branch.mcc
          )
)

# Create new plot layer

par(new=TRUE)

plot(0
     ,type = "n"
     ,xlim = setting.range.xlim
     #,ylim = setting.range.ylim
     ,ylim = c(-7,2)
     ,xaxs = "i"
     ,yaxs = "i"
     ,axes = FALSE
     ,xlab=""
     ,ylab=""
)

# Add MCC loess line

lw1 <- 
  loess(log(coelacanth.rates.morph.mcc) ~ coelacanth.midpoint.ages.mcc
    , span = 0.35
    )

j <- order(coelacanth.midpoint.ages.mcc)

lines(coelacanth.midpoint.ages.mcc[j]
      , lw1$fitted[j]
      , col = "#009988"
      , lwd = setting.lwd.loess
      )

# Make axes visible

axis(1)
axis(2)

# Add axis text

mtext("Log Morphological Rate\n(% Change / Lineage Myr)"
  , side = 2, line = 3
  )

# Output

dev.off()


## Lungfishes

# Axis ranges

setting.range.xlim <- 
  c(round(max(unlist(lungfish.node.depths.mcc)),digits = 1), 0)

setting.range.ylim <- 
  ceiling(range(log(unlist(lungfish.rates.morph.posterior))))

# Output info

setwd(wd.lf)

#pdf(
# paste("./Lungfishes_MorphRates_BEAST_LineSegmentPlot-Posterior.pdf"
#, sep = ""), height = 8, width = 12)

png(
  filename = 
    "./Lungfishes_MorphRates_BEAST_LineSegmentPlot-Posterior.png"
    , width = 1280
    , height = 720
  )

# Plot area

par(mar=c(5,4,4,5)+.1)

# Create plot layer

plot(0
     ,type = "n"
     ,xlim = setting.range.xlim
     #,ylim = setting.range.ylim
     ,ylim = c(-7,2)
     ,yaxt = "n"
     ,xaxs = "i"
     ,yaxs = "i"
     ,axes = FALSE
     ,xlab=""
     ,ylab=""
) 

# Add posterior branch line segments

for (i in 1:length(lungfish.rates.morph.posterior)) {
  segments(x0 = lungfish.node.depths.posterior[[i]][,2]
           ,y0 = log(lungfish.rates.morph.posterior[[i]])
           ,x1 = lungfish.node.depths.posterior[[i]][,1]
           ,y1 = log(lungfish.rates.morph.posterior[[i]])
           ,lwd = setting.lwd.branch.posterior
           ,col = adjustcolor(setting.color
            , alpha.f = setting.alpha.branch.posterior
            )
  )  
}

# Create new plot layer

par(new=TRUE)

plot(0
     ,type = "n"
     ,xlim = setting.range.xlim
     #,ylim = setting.range.ylim
     ,ylim = c(-7,2)
     ,xaxs = "i"
     ,yaxs = "i"
     ,axes = FALSE
     ,xlab=""
     ,ylab=""
)

# Add mcc branch line segments

segments(x0 = lungfish.node.depths.mcc[,2]
         ,y0 = log(lungfish.rates.morph.mcc)
         ,x1 = lungfish.node.depths.mcc[,1]
         ,y1 = log(lungfish.rates.morph.mcc)
         ,lwd = setting.lwd.branch.mcc
         ,col = adjustcolor(setting.color
          , alpha.f = setting.alpha.branch.mcc
          )
)

# Create new plot layer

par(new=TRUE)

plot(0
     ,type = "n"
     ,xlim = setting.range.xlim
     #,ylim = setting.range.ylim
     ,ylim = c(-7,2)
     ,xaxs = "i"
     ,yaxs = "i"
     ,axes = FALSE
     ,xlab=""
     ,ylab=""
)

# Add MCC loess line

lw1 <- 
  loess(
    log(lungfish.rates.morph.mcc) ~ lungfish.midpoint.ages.mcc
    , span = 0.35
    )

j <- order(lungfish.midpoint.ages.mcc)

lines(lungfish.midpoint.ages.mcc[j]
      , lw1$fitted[j]
      , col = "#ED7733"
      , lwd = setting.lwd.loess
      )

# Make axes visible

axis(1)
axis(2)

# Add axes text

mtext("Log Morphological Rate\n(% Change / Lineage Myr)"
  , side = 2
  , line = 3
  )

# Output

dev.off()


## Holosteans

# Axis ranges

setting.range.xlim <- 
  c(round(max(unlist(holostean.node.depths.mcc)),digits = -1), 0)

setting.range.ylim <- 
  ceiling(range(log(unlist(holostean.rates.morph.posterior))))

# Output info

setwd(wd.lf)

#pdf(
# paste("./Holosteans_MorphRates_BEAST_LineSegmentPlot-Posterior.pdf"
# , sep = ""), height = 8, width = 12)

png(
  filename = 
    "./Lungfishes_MorphRates_BEAST_LineSegmentPlot-Posterior.png"
    , width = 1280
    , height = 720
    )

# Plot area

par(mar=c(5,4,4,5)+.1)

# Create plot layer

plot(0
     ,type = "n"
     ,xlim = setting.range.xlim
     #,ylim = setting.range.ylim
     ,ylim = c(4,-4)
     ,yaxt = "n"
     ,xaxs = "i"
     ,yaxs = "i"
     ,axes = FALSE
     ,xlab=""
     ,ylab=""
) 

# Add posterior branch line segments

for (i in 1:length(holostean.rates.morph.posterior)) {
  segments(x0 = holostean.node.depths.posterior[[i]][,2]
           ,y0 = log(holostean.rates.morph.posterior[[i]])
           ,x1 = holostean.node.depths.posterior[[i]][,1]
           ,y1 = log(holostean.rates.morph.posterior[[i]])
           ,lwd = setting.lwd.branch.posterior
           ,col = adjustcolor(setting.color
            , alpha.f = setting.alpha.branch.posterior
            )
  )  
}

# Create new plot layer

par(new=TRUE)

plot(0
     ,type = "n"
     ,xlim = setting.range.xlim
     #,ylim = setting.range.ylim
     ,ylim = c(4,-4)
     ,xaxs = "i"
     ,yaxs = "i"
     ,axes = FALSE
     ,xlab=""
     ,ylab=""
)

# Add mcc branch line segments

segments(x0 = holostean.node.depths.mcc[,2]
         ,y0 = log(holostean.rates.morph.mcc)
         ,x1 = holostean.node.depths.mcc[,1]
         ,y1 = log(holostean.rates.morph.mcc)
         ,lwd = setting.lwd.branch.mcc
         ,col = adjustcolor(setting.color
          , alpha.f = setting.alpha.branch.mcc
          )
)

# Create new plot layer

par(new=TRUE)

plot(0
     ,type = "n"
     ,xlim = setting.range.xlim
     #,ylim = setting.range.ylim
     ,ylim = c(4,-4)
     ,xaxs = "i"
     ,yaxs = "i"
     ,axes = FALSE
     ,xlab=""
     ,ylab=""
)

# Add MCC loess line

lw1 <- 
  loess(log(holostean.rates.morph.mcc) ~ holostean.midpoint.ages.mcc
    , span = 0.35
    )

j <- order(holostean.midpoint.ages.mcc)

lines(holostean.midpoint.ages.mcc[j]
      , lw1$fitted[j]
      , col = "#38BBEC"
      , lwd = setting.lwd.loess
      )

# Make axes visible

axis(1)
axis(4)

# Add axis text

mtext("Log Morphological Rate\n(% Change / Lineage Myr)"
  , side = 2, line = 3
  )

# Output

dev.off()



### Plot average rate through time for MCC and sampled posterior trees


## Coelacanths

# Axis ranges

setting.range.xlim <- 
  c(round(max(unlist(coelacanth.node.depths.posterior)),digits = -1), 0) 
setting.range.ylim <- 
  ceiling(range(log(unlist(coelacanth.rates.morph.posterior))))

# Output info

setwd(wd.lf)

#pdf(
# paste(
#   "./Figure_Coelacanths_MorphRates_BEAST_LoessLinesPlot-Posterior.pdf"
#   , sep = ""), height = 8, width = 12)

png(
  filename = 
  "./Figure_Coelacanths_MorphRates_BEAST_LoessLinesPlot-Posterior.png"
  , width = 1280
  , height = 720
  )

# Set up window
par(mar=c(5,4,4,5)+.1)

# Create plot area
plot(0
     ,type = "n"
     ,xlim = setting.range.xlim
     #,ylim = setting.range.ylim
     ,ylim = c(-7,2) #this is manual so all plots have the same y limit
     ,yaxt = "n"
     ,xaxs = "i"
     ,yaxs = "i"
     ,axes = FALSE
     ,xlab=""
     ,ylab=""
) 

# Loop through all posterior trees to plot loess lines for each
for (i in 1:length(coelacanth.rates.morph.posterior)){
  
  lw1 <- 
    loess(
      log(coelacanth.rates.morph.posterior[[i]]) ~ 
        coelacanth.midpoint.ages.posterior[[i]]
        , span = 0.35
        )
  
  j <- order(coelacanth.midpoint.ages.posterior[[i]])
  
  lines(coelacanth.midpoint.ages.posterior[[i]][j]
        , lw1$fitted[j]
        , col = adjustcolor(setting.color
          , alpha.f = setting.alpha.branch.posterior
          )
        , lwd = setting.lwd.loess
  )
  
}

# Loess line from the MCC is plotted now so it is over everything else
lw1 <- 
  loess(
    log(coelacanth.rates.morph.mcc) ~ 
      coelacanth.midpoint.ages.mcc
      , span = 0.35
      )

j <- order(coelacanth.midpoint.ages.mcc)

lines(coelacanth.midpoint.ages.mcc[j]
      , lw1$fitted[j]
      , col = "#009988"
      , lwd = setting.lwd.loess
)

# Make the axes visible
axis(1)
axis(2)

# Add text to the figure
mtext("Log Morphological Rate\n(% Change / Lineage Myr)"
  , side = 2, line = 2
  )

# Finish and output
dev.off()


## Lungfishes

# Axis ranges

setting.range.xlim <- 
  c(round(max(unlist(lungfish.node.depths.posterior)),digits = 1), 0)

setting.range.ylim <- 
  ceiling(range(log(unlist(lungfish.rates.morph.posterior))))

# Output info

setwd(wd.lf)

#pdf(
# paste(
#   "./Figure_Lungfishes_MorphRates_BEAST_LoessLinesPlot-Posterior.pdf"
#   , sep = ""), height = 8, width = 12)

png(
  filename = 
    "./Figure_Lungfishes_MorphRates_BEAST_LoessLinesPlot-Posterior.png"
    , width = 1280
    , height = 720
    )

# Set up window
par(mar=c(5,4,4,5)+.1)

# Create plot area
plot(0
     ,type = "n"
     ,xlim = setting.range.xlim
     #,ylim = setting.range.ylim
     ,ylim = c(-7,2) #this is manual so all plots have the same y limit
     ,yaxt = "n"
     ,xaxs = "i"
     ,yaxs = "i"
     ,axes = FALSE
     ,xlab=""
     ,ylab=""
) 

# Loop through all posterior trees to plot loess lines for each
for (i in 1:length(lungfish.rates.morph.posterior)){
  
  lw1 <- 
    loess(
      log(lungfish.rates.morph.posterior[[i]]) ~ 
        lungfish.midpoint.ages.posterior[[i]]
        , span = 0.35
        )
  
  j <- order(lungfish.midpoint.ages.posterior[[i]])
  
  lines(lungfish.midpoint.ages.posterior[[i]][j]
        , lw1$fitted[j]
        , col = adjustcolor(setting.color
          , alpha.f = setting.alpha.branch.posterior
          )
        , lwd = setting.lwd.loess
  )
  
}

# Loess line from the MCC is plotted now so it is over everything else
lw1 <- 
  loess(
    log(lungfish.rates.morph.mcc) ~ 
      lungfish.midpoint.ages.mcc
      , span = 0.35
      )

j <- order(lungfish.midpoint.ages.mcc)

lines(lungfish.midpoint.ages.mcc[j]
      , lw1$fitted[j]
      , col = "#ED7733"
      , lwd = setting.lwd.loess
)

# Make the axes visible
axis(1)
axis(2)

# Add text to the figure
mtext("Log Morphological Rate\n(% Change / Lineage Myr)"
  , side = 2, line = 2
  )

# Finish and output
dev.off()


## Holosteans

# Axis ranges

setting.range.xlim <- 
  c(round(max(unlist(holostean.node.depths.posterior)),digits = -1), 0) 
setting.range.ylim <- 
  ceiling(range(log(unlist(holostean.rates.morph.posterior))))

# Output info

setwd(wd.lf)
#pdf(
# paste(
#   "./Figure_Holosteans_MorphRates_BEAST_LoessLinesPlot-Posterior.pdf"
#   , sep = ""), height = 8, width = 12)

png(
  filename = 
    "./Figure_Holosteans_MorphRates_BEAST_LoessLinesPlot-Posterior.png"
    , width = 1280
    , height = 720
    )

# Set up window
par(mar=c(5,4,4,5)+.1)

# Create plot area
plot(0
     ,type = "n"
     ,xlim = setting.range.xlim
     #,ylim = setting.range.ylim
     ,ylim = c(-7,2) #this is manual so all plots have the same y limit
     ,yaxt = "n"
     ,xaxs = "i"
     ,yaxs = "i"
     ,axes = FALSE
     ,xlab=""
     ,ylab=""
) 

# Loop through all posterior trees to plot loess lines for each
for (i in 1:length(holostean.rates.morph.posterior)){
  
  lw1 <- 
    loess(
      log(holostean.rates.morph.posterior[[i]]) ~ 
        holostean.midpoint.ages.posterior[[i]]
        , span = 0.35
        )
  
  j <- order(holostean.midpoint.ages.posterior[[i]])
  
  lines(holostean.midpoint.ages.posterior[[i]][j]
        , lw1$fitted[j]
        , col = adjustcolor(setting.color
          , alpha.f = setting.alpha.branch.posterior
          )
        , lwd = setting.lwd.loess
  )
  
}

# Loess line from the MCC is plotted now so is over everything else
lw1 <- 
  loess(
    log(holostean.rates.morph.mcc) ~ 
      holostean.midpoint.ages.mcc
      , span = 0.35
      )

j <- order(holostean.midpoint.ages.mcc)

lines(holostean.midpoint.ages.mcc[j]
      , lw1$fitted[j]
      , col = "#38BBEC"
      , lwd = setting.lwd.loess
)

# Make the axes visible
axis(1)
axis(2)

# Add text to the figure
mtext("Log Morphological Rate\n(% Change / Lineage Myr)"
  , side = 2, line = 2
  )

# Finish and output
dev.off()


#####


#### Import the raw data to process shapes ####



### Import Raw Data
### NOTE: this step assumes you have digitized data in a folder called 
### 'shapes_[taxon]' within the working directory


if(exists("set.group") == FALSE || is.null(set.group) == TRUE){
  
  print(
  "You haven't let the script know what clade you're working with yet!"
  )
  
} else if (set.group == "coelacanths") {
  
  # Save the tree objects as another name to manipulate in the future
  tree.mcc.shapeanalyses <- tree.mcc.coelacanth.morph
  trees.posterior.shapeanalyses <- trees.posterior.coelacanth.morph
  
  # Read the raw shape data
  shapes <- readShapes('shapes_coelacanths')
  
  # Create points at even spacing (i.e., semilandmarks) along the curves
  # NOTE: this can vary from group to group, so make sure to test
  #curve_head = 25
  #curve_dorsal-1 = 5
  #curve_dorsal-2 = 5
  #curve_dorsal-3 = 15
  #curve_dorsal-4 = 10
  #curve_fin_epichordal = 30
  #curve_fin_hypochordal = 30
  #curve_ventral-1 = 5
  #curve_ventral-2 = 30
  semilandmarks <- c(25,5,5,15,10,30,30,5,30)
  
  # Read file listing taxa you don't have shape data for
  tip.drop.noshapedata <- 
    readLines("tree_taxa_coelacanth_noshapedata.txt")
  
  
} else if (set.group == "lungfishes") {
  
  # Save the tree objects as another name to manipulate in the future
  tree.mcc.shapeanalyses <- tree.mcc.lungfish.morph
  trees.posterior.shapeanalyses <- trees.posterior.lungfish.morph
  
  # Read the raw shape data
  shapes <- readShapes('shapes_lungfishes')
  
  # Create points at even spacing (i.e., semilandmarks) along the curves
  # NOTE: this can vary from group to group, so make sure to test
  #curve_head = 15
  #curve_dorsal = 20
  #curve_fin_dorsal = 30
  #curve_fin_anal = 20
  #curve_ventral = 30
  semilandmarks <- c(15,20,30,20,30)
  
  # Read file listing taxa you don't have shape data for
  tip.drop.noshapedata <- 
    readLines("tree_taxa_lungfish_noshapedata.txt")
  
  
} else if (set.group == "holosteans") {
  
  # Save the tree objects as another name to manipulate in the future
  tree.mcc.shapeanalyses <- tree.mcc.holostean.morph
  trees.posterior.shapeanalyses <- trees.posterior.holostean.morph
  
  # Read the raw shape data
  shapes <- readShapes('shapes_holosteans')
  
  # Create points at even spacing (i.e., semilandmarks) along the curves
  # NOTE: this can vary from group to group, so make sure to test
  #curve_head = 25
  #curve_dorsal_1 = 15
  #curve_dorsal_2 = 10
  #curve_dorsal_3 = 10
  #curve_fin_caudal = 20
  #curve_ventral_1 = 10
  #curve_ventral_2 = 10
  #curve_ventral_3 = 20
  semilandmarks <- c(25,15,10,10,20,10,10,20)
  
  # Read file listing taxa you don't have shape data for
  tip.drop.noshapedata <- 
    readLines("tree_taxa_holostean_noshapedata.txt")
  
  
}




#####


#### Process and standardize shape data ####



### Set the landmark and semilandmark data for all subsequent analyses

shapes <- readland.shapes(shapes, nCurvePts = semilandmarks)



### Remove data missing from both datasets


## Compare tree tip and landmark labels

setdiff(tree.mcc.shapeanalyses$tip.label,names(shapes$landmarks))


## Drop tips from MCC tree that are not in the shape data

tree.mcc.shapeanalyses.reduced <- 
  drop.tip(tree.mcc.shapeanalyses, tip.drop.noshapedata)


## Drop tips from the posterior trees that are not in the shape data

class(trees.posterior.shapeanalyses) <- "multiPhylo"

trees.posterior.shapeanalyses.reduced <- 
  drop.tip.multiPhylo(trees.posterior.shapeanalyses
    , tip.drop.noshapedata
    )


## Clean and standardize data so it agrees between shape and phylo
## NOTE: this section will change if naming conventions or sampling
## change, so it's important to do it right, or else the model fitting
## will have errors

if(exists("set.group") == FALSE || is.null(set.group) == TRUE) {
  
  print(
  "You haven't let the script know what clade you're working with yet!"
  )
  
  
} else if (set.group == "coelacanths") {
  
  # Remove shapedata for taxa not in the tree
  #"Changxingia" = 4
  #"Chaohuichthys" = 5
  #"Macropomoides_orientalis_rec" = 21
  shapes$landmarks <- shapes$landmarks[-c(4,5,21)]
  
  # Fix shape data taxa so they match Genus (same as tree)
  names(shapes$landmarks) <- 
    gsub("(.*)_.*_.*","\\1",names(shapes$landmarks))
  
  # Fix any remaining problematic taxa in the shape data
  names(shapes$landmarks) <- 
    gsub("(Heptanema)_.*","\\1",names(shapes$landmarks))
  
  
} else if (set.group == "lungfishes") {
  
  # Fix shape data taxa so they match Genus_species (same as tree)
  names(shapes$landmarks) <- 
    gsub("(.*_.*)_.*","\\1",names(shapes$landmarks))
  
  # Fix any taxon you want to generalize to Genus_spp in the shape data
  names(shapes$landmarks) <- 
    gsub("(Protopterus)_.*","\\1_spp",names(shapes$landmarks))

  names(shapes$landmarks) <- 
    gsub("(Soederberghia)_.*","\\1_spp",names(shapes$landmarks))

  names(shapes$landmarks) <- 
    gsub("(Sagenodus)_.*","\\1_spp",names(shapes$landmarks))
  
  # Fix any taxon you want to generalize to Genus_spp in the phylo MCC 
  # data
  tree.mcc.shapeanalyses.reduced$tip.label <- 
    gsub("(Protopterus)_.*","\\1_spp"
      , tree.mcc.shapeanalyses.reduced$tip.label
      )

  tree.mcc.shapeanalyses.reduced$tip.label <- 
    gsub("(Soederberghia)_.*","\\1_spp"
      , tree.mcc.shapeanalyses.reduced$tip.label
      )

  tree.mcc.shapeanalyses.reduced$tip.label <- 
    gsub("(Sagenodus)_.*","\\1_spp"
      , tree.mcc.shapeanalyses.reduced$tip.label
      )
  
  # Fix any taxon you want to generalize to Genus_spp in the phylo 
  # posterior data
  for (i in 1:length(trees.posterior.shapeanalyses.reduced)){

    trees.posterior.shapeanalyses.reduced[[i]]$tip.label <- 
      gsub("(Protopterus)_.*","\\1_spp"
        , trees.posterior.shapeanalyses.reduced[[i]]$tip.label
        )

    trees.posterior.shapeanalyses.reduced[[i]]$tip.label <- 
      gsub("(Soederberghia)_.*","\\1_spp"
        , trees.posterior.shapeanalyses.reduced[[i]]$tip.label
        )

    trees.posterior.shapeanalyses.reduced[[i]]$tip.label <- 
      gsub("(Sagenodus)_.*","\\1_spp"
        , trees.posterior.shapeanalyses.reduced[[i]]$tip.label
        )
  }
  
  
} else if (set.group == "holosteans") {
  
  # Remove shapedata for taxa not in the tree
  #"Plesiofuro_mingshuica" = 43
  shapes$landmarks <- shapes$landmarks[-c(43)]
  
  # Fix shape data names so they match Genus_species (same as tree)
  names(shapes$landmarks) <- 
    gsub("(.*_.*)_.*","\\1",names(shapes$landmarks))
  
  # Fix any remaining problematic taxa in the phylo MCC data
  tree.mcc.shapeanalyses.reduced$tip.label <- 
    gsub("Ophiopsiella_(.*)","Ophiopsis_\\1"
      , tree.mcc.shapeanalyses.reduced$tip.label
      )
  
  # Fix any remaining problematic taxa in the phylo posterior data
  for (i in 1:length(trees.posterior.shapeanalyses.reduced)){
    trees.posterior.shapeanalyses.reduced[[i]]$tip.label <- 
      gsub("Ophiopsiella_(.*)","Ophiopsis_\\1"
        , trees.posterior.shapeanalyses.reduced[[i]]$tip.label
        )
  }
  
  
}




#####


#### CHECKPOINT: Export pruned, standardized shape data ####



### Store the outputs according to the group you're testing and 
### standardize


if (set.group == "coelacanths") {
  
  # Back up the further pruned and standardized MCC tree 
  tree.mcc.shapeanalyses.reduced.coelacanth <- 
    tree.mcc.shapeanalyses.reduced
  
  saveRDS(tree.mcc.shapeanalyses.reduced.coelacanth
          ,"tree_mcc_coelacanth_shapeanalyses.RDS"
          )
  
  # Back up the further pruned and standardized posterior trees
  trees.posterior.shapeanalyses.reduced.coelacanth <- 
    trees.posterior.shapeanalyses.reduced
  
  saveRDS(trees.posterior.shapeanalyses.reduced.coelacanth
          ,"trees_posterior_coelacanth_shapeanalyses.RDS"
          )
  
  # Back up the processed 2d shape data
  shapes.coelacanth <- shapes
  
  saveRDS(shapes.coelacanth,"data_shape_coelacanth.RDS")
  
  
} else if (set.group == "lungfishes") {
  
  # Back up the further pruned and standardized MCC tree 
  tree.mcc.shapeanalyses.reduced.lungfish <- 
    tree.mcc.shapeanalyses.reduced
  
  saveRDS(tree.mcc.shapeanalyses.reduced.lungfish
          ,"tree_mcc_lungfish_shapeanalyses.RDS"
  )
  
  # Back up the further pruned and standardized posterior trees
  trees.posterior.shapeanalyses.reduced.lungfish <- 
    trees.posterior.shapeanalyses.reduced
  
  saveRDS(trees.posterior.shapeanalyses.reduced.lungfish
          ,"trees_posterior_lungfish_shapeanalyses.RDS"
  )
  
  # Back up the processed 2d shape data
  shapes.lungfish <- shapes
  
  saveRDS(shapes.lungfish,"data_shape_lungfish.RDS")
  
  
} else if (set.group == "holosteans") {
  
  # Back up the further pruned and standardized MCC tree 
  tree.mcc.shapeanalyses.reduced.holostean <- 
    tree.mcc.shapeanalyses.reduced
  
  saveRDS(tree.mcc.shapeanalyses.reduced.holostean
          ,"tree_mcc_holostean_shapeanalyses.RDS"
  )
  
  # Back up the further pruned and standardized posterior trees
  trees.posterior.shapeanalyses.reduced.holostean <- 
    trees.posterior.shapeanalyses.reduced
  
  saveRDS(trees.posterior.shapeanalyses.reduced.holostean
          ,"trees_posterior_holostean_shapeanalyses.RDS"
  )
  
  # Back up the processed 2d shape data
  shapes.holostean <- shapes
  
  saveRDS(shapes.holostean,"data_shape_holostean.RDS")
  
  
} else {
  
  print(
    "Something went wrong. Have you set the clade you're working with?"
    )
  
}




#####
