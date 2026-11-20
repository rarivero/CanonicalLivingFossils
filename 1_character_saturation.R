#
##
#### CHARACTER SATURATION ANALYSES 
#### USING MORPHOLOGICAL MATRICES + CONSENSUS TREES
##
#

#### Script by Rafael A. Rivero Vega
#### Based off original in Brocklehurst et al 2021
#### (https://doi.org/10.1016/j.cub.2021.04.044), with further
#### details in the manuscript methods and acknowledgements.




#### Load the required libraries ####

if(!'Claddis'%in%installed.packages()
  ) install.packages("Claddis"); library(Claddis)

if(!'adephylo'%in%installed.packages()
  ) install.packages("adephylo"); library(adephylo)

if(!'drc'%in%installed.packages()
  ) install.packages("drc"); library(drc)

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

#### Set up the data for the analyses ####



### Identify the clade being worked with
### NOTE: this setting allows for the import of all the data, preps 
### objects, and adjusts the way outputs are saved/backed up further
### down


set.group <- "coelacanths"
#set.group <- "lungfishes"
#set.group <- "holosteans"



### Import additional data
### NOTE: this will use the read.annotated.nexus function to import the
### mcc tree, morphological matrices, and a list of "affinities"


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
  
  # Morphological matrix in Nexus format (modified from the original 
  # as noted within the file itself)
  matrix.morph <-
    read_nexus_matrix(
      "./data_coelacanths_morphology-Torino-etal_matrix.nex"
    )
  
  # Taxon affinities (clade and inclusion data) in CSV format
  affinities <- 
    as.matrix(read.csv("./taxa_coelacanths_affinities.csv",header=T))
  
  
} else if (set.group == "lungfishes") {
  
  # Maximum Clade Credibility tree (TreeAnnotator output)
  tree <- 
    read.annotated.nexus(
      "./beast_lungfishes_morphology_Lloyd-etal_tree-mcc.tre"
      )
  
  # Morphological matrix in Nexus format (modified from the original 
  # as noted within the file itself)
  matrix.morph <-
    read_nexus_matrix(
      "./data_lungfishes_morphology_Lloyd-etal_matrix.nex"
    )
  
  # Taxon affinities (clade and inclusion data) in CSV format
  affinities <- 
    as.matrix(read.csv("./taxa_lungfishes_affinities.csv",header=T))
  
  
} else if (set.group == "holosteans") {
  
  # Maximum Clade Credibility tree (TreeAnnotator output)
  tree <- 
    read.annotated.nexus(
    "./beast_neopterygii_morphology_LopezArbarello-Sferco_tree-mcc.tre"
    )
  
  # Morphological matrix in Nexus format (modified from the original 
  # as noted within the file itself)
  matrix.morph <-
    read_nexus_matrix(
      "./data_neopterygii_morphology_LopezArbarello-Sferco_matrix.nex"
      )
  
  # Taxon affinities (clade and inclusion data) in CSV format
  affinities <- 
    as.matrix(read.csv("./taxa_holosteans_affinities.csv",header=T))
  
  
}



#####

#### Prepare MCC tree data for saturation analyses ####



### Extract and transform morphological rate data
### NOTE: this section is superficially similar to the beginning of
### the 0_process_raw_data.R script and does overwrite some of the 
### objects that were previously created if you ran that first, so
### ensure that you've reached the end of that script and all the 
### checkpoints before running this one. However, there are extra 
### steps in this section that mathematically change edge length data,
### so don't assume it can be skipped!


## Preserve edge lengths with their ascending numbers to use as an
## index further down
## NOTE: it is crucial that this is done before dropping tips!

edge.lengths <- tree$edge.length


## Transform branch lengths according to rate data 

# Create new object to store the transformed data 
edge.lengths.new <- NULL

# Use the rates to transform the branch lengths into units of character
# change (rate = change / time; so change = rate * time)
for(i in 1:length(edge.lengths)){
  edge.lengths.new[i] <- 
    edge.lengths[i] * tree$annotations[[i]]$rate_median
}

# Bind the new edge lengths together with their numerical values
# NOTE: it is crucial that you preserve the numbers in ascending order
# before joining with the new edge lengths here! 
edge.lengths <- 
  cbind(1:length(edge.lengths),edge.lengths,edge.lengths.new)



### Create a final tree object with the adjusted annotations


## Exclude tips using affinity data

# Drop taxa not to be included in subsequent analyses from tree
tree.drop <- drop.tip(tree,affinities[which(affinities[,3]=="0"),1])

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

# Replace the edge.lengths of the drop.tips tree with the transformed
# edge.lengths object
tree.drop$edge.length <- 
  edge.lengths[!(edge.lengths[,2] %in% edge.lengths.drop),3]

# Create a new standardized object so it works with the next steps
tree.morph <- tree.drop



### Standardize taxon labels
### NOTE: unlike in 0_process_raw_data.R, this standardization will be
### performed on the MCC tree (tree$tip.label) and the morphological
### matrix (matrix.morph), which will result in differences between
### the object used in these analyses and the rest of the scripts.
### Due to this, unlike previous times, the objects standardized here
### won't be backed up. They can, however, be recreated very easily.
### As before, the purpose of this section is to remove all extra info
### beyond Genus_species that can be found prepended or appended to
### each taxon (e.g., "F_Psarolepis_romeri_410.8" instead becomes
### "Psarolepis_romeri"). In our case, that means this script is 
### written to only apply this to lungfishes and holosteans, but not 
### to coelacanths.


if(exists("set.group") == FALSE || is.null(set.group) == TRUE){
  
  print(
  "You haven't let the script know what clade you're working with yet!"
  )
  
} else if (set.group == "coelacanths") {
  
  # Standardize tips to 'Genus_species' in the MCC tree
  
  #Save the MCC tree to its own object without any adjustments
  morph.tree.standardized <- tree.morph
  
  # Standardize tips to 'Genus_species' in the morphological matrix
  
  #Save the matrix to its own object without any adjustments
  morph.matrix.standardized <- matrix.morph
  
} else if (set.group == "lungfishes" || set.group == "holosteans") {
  
  # Standardize tips to 'Genus_species' in the MCC tree
  
  #Save the MCC tree to its own object
  morph.tree.standardized <- tree.morph
  
  #"F_" fossil identifier
  morph.tree.standardized$tip.label <- 
    gsub("F_", "", morph.tree.standardized$tip.label)
  
  #Ages in the name
  morph.tree.standardized$tip.label <- 
    gsub("(.*)_[0-9].*", "\\1", morph.tree.standardized$tip.label)
  
  # Standardize tips to 'Genus_species' in morphological matrix
  
  #Save the matrix to its own object
  morph.matrix.standardized <- matrix.morph
  
  #"F_" fossil identifier
  rownames(morph.matrix.standardized$matrix_1$matrix) <- 
    gsub("F_",""
         , rownames(morph.matrix.standardized$matrix_1$matrix)
         )
  
  #Ages in the name
  rownames(morph.matrix.standardized$matrix_1$matrix) <- 
    gsub("(.*)_[0-9].*","\\1"
         , rownames(morph.matrix.standardized$matrix_1$matrix)
         )
  
} else {
  
  print("Something went wrong. Have you set your clade correctly?")
  
}


#####

#### Prepare morphological matrices and perform calculations ####



### Calculate moprhological distances and drop excluded taxa


## Calculation using Claddis

dist <- 
  calculate_morphological_distances(morph.matrix.standardized
                                      ,distance_metric = "mord"
                                      ,distance_transformation = "none"
                                      )


## Trim excluded taxa

morph.dist <- dist[[2]][morph.tree.standardized$tip.label
                        , morph.tree.standardized$tip.label
                        ]

morph.dist <- trim_matrix(morph.dist)[[1]]


## Sometimes, Claddis doesn't replace all morph distances >1 with 1 
## automatically even though it is supposed to in unordered datasets.
## Check if that has happened in your case using the following:
## morph.dist[morph.dist > 1]
## If you see values >1, simply use the command below to adjust the
## values manually

#morph.dist[morph.dist > 1] # check
#morph.dist[morph.dist > 1] <- 1 # replace



### Calculate patristic distances


## Calculation using Claddis

pat.dist <- as.matrix(distTips(morph.tree.standardized
                               , method = "patristic")
                      )


## Add names from the morphological matrix to the object

pat.dist <- pat.dist[rownames(morph.dist), colnames(morph.dist)]



### Perform pairwise comparison


## Create an empty matrix object to fill with data

pairs <- matrix(nrow=0,ncol=4)


## Add column names

colnames(pairs)<-c("Taxon 1","Taxon 2","Morph dist","Patrist dist")


## Loop through the morph.dist object and extract data

# Identify all pairs, then pull out their calculated pairwise distances
# (using row and column information)
for(i in 1:nrow(morph.dist)){

  for(j in i:ncol(morph.dist)){

    if(i!=j){
      tax.1 <- rownames(morph.dist)[i]
      tax.2 <- colnames(morph.dist)[j]
      pair.morph <- morph.dist[tax.1,tax.2]
      pair.pat <- pat.dist[tax.1,tax.2]
      pairs <- rbind(pairs,c(tax.1,tax.2,pair.morph,pair.pat))
    }

  }

}



### Fit a regression curve


## Bring patristic and morphological distances together

dists <- 
  cbind(as.numeric(pairs[,"Morph dist"])
        , as.numeric(pairs[,"Patrist dist"])
      )


## Fit loess line

loess.line <- loess(V1~V2,as.data.frame(dists))



### Find the point of saturation


## Fit a Michalis Menten curve

mm.curve <- drm(V1~V2, data = as.data.frame(dists), fct = MM.2())


## Calculate the observed Vmax (aka point of saturation)

vmax.saturation <- coef(mm.curve)[1]




#####

#### CHECKPOINT: Export character saturation data ####



### Store the outputs according to the group you're testing and
### standardize


if (set.group == "coelacanths") {
  
  # Back up character saturation data
  dists.coelacanths <- dists

  # Export character saturation data
  saveRDS(dists
    , 'saturation_coelacanth.RDS'
    )
  
  # Back up vmax data
  vmax.saturation.coelacanths <- vmax.saturation

  # Export vmax data
  saveRDS(vmax.saturation
    , 'saturation_coelacanth_vmax.RDS'
    )
  
  
} else if (set.group == "lungfishes") {
  
  # Back up character saturation data
  dists.lungfishes <- dists

  # Export character saturation data
  saveRDS(dists
    , 'saturation_lungfish.RDS'
    )
  
  # Back up vmax data
  vmax.saturation.lungfishes <- vmax.saturation

  # Export vmax data
  saveRDS(vmax.saturation
    , 'saturation_lungfish_vmax.RDS'
    )
  
  
} else if (set.group == "holosteans") {
  
  # Back up character saturation data
  dists.holosteans <- dists

  # Export character saturation data
  saveRDS(dists
    , 'saturation_holostean.RDS'
    )
  
  # Back up vmax data
  vmax.saturation.holosteans <- vmax.saturation

  # Export vmax data
  saveRDS(vmax.saturation
    , 'saturation_holostean_vmax.RDS'
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
    
    # Import character saturation data
    dists <- readRDS('saturation_coelacanth.RDS')

    # Import Vmax data 
    vmax.saturation <- readRDS('saturation_coelacanth_vmax.RDS')
    
  } else if (set.group == "lungfishes") {
    
    # Import character saturation data
    dists <- readRDS('saturation_lungfish.RDS')

    # Import Vmax data 
    vmax.saturation <- readRDS('saturation_lungfish_vmax.RDS')
    
  } else if (set.group == "holosteans") {
    
    # Import character saturation data
    dists <- readRDS('saturation_holostean.RDS')

    # Import Vmax data 
    vmax.saturation <- readRDS('saturation_holostean_vmax.RDS')
    
  }
  
}



#####

#### Plot the results ####
#### NOTE: the plots are manually configurable, so they're not
#### automated and should be adjusted as needed



### Prepare settings


## Adjust window

par(mar = c(4,4,1,1))
#par(mar = c(1,1,1,1)) # reset if needed



### Plot for Coelacanths and export


## Ensure working directory is correct

setwd(wd.lf)


## Open the device
## NOTE: only one dev can be outputting at a time. If you want to
## output to PDF *and* PNG, it should be run twice, one at a time,
## all the way until the dev.off() command.

# Use and adjust the below if exporting to PDF
pdf(
  paste("./Coelacanths_CharacterSaturation_ScatterLoessPlot.pdf"
        , sep = ""), height = 8, width = 12)

# Use and adjust the below if exporting to PNG
#png(
#  filename = 
#    "./Coelacanths_CharacterSaturation_ScatterLoessPlot.png"
#  , width = 1280
#  , height = 720
#)


## Create the plot

# Scatterplot the character saturation
plot(dists[,2]
     , dists[,1]
     , pch = 16
     , cex = 0.3
     , col = "blue"
     , bty = "l"
     , main = "Character Saturation: Coelacanths"
     #, main = "Character Saturation: Lungfishes"
     #, main = "Character Saturation: Holosteans"
     , xlab = "Patristic Distance"
     , ylab = "Morphological Dissimilarity"
     , xlim = c(min(pat.dist),max(pat.dist))
     , ylim = c(min(morph.dist),max(morph.dist))
     )

# Overlay loess line on scatterplot
lines(sort(dists[,2]),sort(predict(loess.line)),col="black",lwd=2)


## Close the device

dev.off()



### Plot for Lungfishes and export


## Ensure working directory is correct

setwd(wd.lf)


## Open the device
## NOTE: only one dev can be outputting at a time. If you want to
## output to PDF *and* PNG, it should be run twice, one at a time,
## all the way until the dev.off() command.

# Use and adjust the below if exporting to PDF
pdf(
  paste("./Lungfishes_CharacterSaturation_ScatterLoessPlot.pdf"
        , sep = ""), height = 8, width = 12)

# Use and adjust the below if exporting to PNG
#png(
#  filename = 
#    "./Lungfishes_CharacterSaturation_ScatterLoessPlot.png"
#  , width = 1280
#  , height = 720
#)


## Create the plot

# Scatterplot the character saturation
plot(dists[,2]
     , dists[,1]
     , pch = 16
     , cex = 0.3
     , col = "blue"
     , bty = "l"
     #, main = "Character Saturation: Coelacanths"
     , main = "Character Saturation: Lungfishes"
     #, main = "Character Saturation: Holosteans"
     , xlab = "Patristic Distance"
     , ylab = "Morphological Dissimilarity"
     , xlim = c(min(pat.dist),max(pat.dist))
     , ylim = c(min(morph.dist),max(morph.dist))
)

# Overlay loess line on scatterplot
lines(sort(dists[,2]),sort(predict(loess.line)),col="black",lwd=2)


## Close the device

dev.off()



### Plot for Holosteans and export


## Ensure working directory is correct

setwd(wd.lf)


## Open the device
## NOTE: only one dev can be outputting at a time. If you want to
## output to PDF *and* PNG, it should be run twice, one at a time,
## all the way until the dev.off() command.

# Use and adjust the below if exporting to PDF
pdf(
  paste("./Holosteans_CharacterSaturation_ScatterLoessPlot.pdf"
        , sep = ""), height = 8, width = 12)

# Use and adjust the below if exporting to PNG
#png(
#  filename = 
#    "./Holosteans_CharacterSaturation_ScatterLoessPlot.png"
#  , width = 1280
#  , height = 720
#)


## Create the plot

# Scatterplot the character saturation
plot(dists[,2]
     , dists[,1]
     , pch = 16
     , cex = 0.3
     , col = "blue"
     , bty = "l"
     #, main = "Character Saturation: Coelacanths"
     #, main = "Character Saturation: Lungfishes"
     , main = "Character Saturation: Holosteans"
     , xlab = "Patristic Distance"
     , ylab = "Morphological Dissimilarity"
     , xlim = c(min(pat.dist),max(pat.dist))
     , ylim = c(min(morph.dist),max(morph.dist))
)

# Overlay loess line on scatterplot
lines(sort(dists[,2]),sort(predict(loess.line)),col="black",lwd=2)


## Close the device

dev.off()


## Reset display settings if needed

par(mar = c(1,1,1,1))



#####
