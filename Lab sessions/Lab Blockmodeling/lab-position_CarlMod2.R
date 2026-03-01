# ===================================================
# R script illustrating how to extract positional
# information from network data
# originally written by Christian Steglich 
# extended and modified by Carl Nordlund
# ===================================================

# Note:
# Some earlier versions (e.g. 3.5.2) don't work well with 'network' package


# The main examples (preselected below) is one of the all-time
# network data classics, Zachary's data set with
# interactions among the members of a karate club.
#
# Zachary W. (1977). An information flow model for conflict
#    and fission in small groups. Journal of Anthropological
#    Research, 33, 452-473. 


# load neccessary packages (install if not already there)
library(network)
library(sna)

# if needed, set working directoy to where data are:
#setwd("<..>")

# Provided data files

#file <- "baboon_t1.txt"          # Monkey interaction (5 male, 7 female); binary, symmetric
#file <- "bakerBinTable6.txt"     # Citation data, social work journals (Baker 1983); binary,symmetric
#file <- "eies_t1.txt"            # Friendship nomination data, EIES conference; ranked, directional
#file <- "hlebec_notesharing.txt" # Note-borrowing data Slovenian students; valued, directional
#file <- "little_league_TI.txt"   # Baseball friendship nomination; binary, directional
#file <- "zachc.txt"              # Social relations, karate club (valued, symmetric)
file <- "zache.txt"               # Social relations, karate club (binary, symmetric)

# Load text file as matrix
mat <- as.matrix(read.table(file, sep="\t", header=TRUE,row.names = 1))

# Convert to network, check directed parameter
net <- as.network(mat, directed=FALSE)  # Change according to specs

# This object is part of the 'network' library - you can inspect the data by simply:
net


# The network can be plotted directly
#gplot(net)

# For symmetric (non-directional) data, we can remove arrowheads and display the labels to make it easier to interpret
gplot(net,displaylabels = TRUE, arrowhead.cex = FALSE)

# For Zachary's network: zache
# The network has two centers - in fact, these are the cores
# of two separate networks into which the karate club split
# after Zachary's study, when a conflict about how to run the
# club broke out between the club's president, a student
# named John A. (actor 34), and the karate teacher employed, 
# Mr. Hi (actor 1).

# to run a positional analysis, we will use the command
# "blockmodel" in the sna package:
?blockmodel
# as you can see in the help file, this creates a blockmodel using an input graph and a partition

# on the help page, we learn that we need to identify a
# notion of equivalence with the help of this function:
?equiv.clust
# There are distance functions for detecting approximate
# structural and regular equivalence, respectively.
# We treat them one by one.

##################################
# STRUCTURAL EQUIVALENCE ANALYSIS
##################################
?sedist
# When the function "sedist" is evaluated
# on a network object, it indicates structural dissimilarity
# (absence of structural equivalence) of actors.
# In this indrect method, this is measured by comparing the row and column
# vectors of each pair of actors in the sociomatrix
# By default, the "hamming" distance is used:
# sedist(dat,method="hamming")
# but other popular alternatives are "euclidean" and "correlation"

# run a positional analysis based on structural equivalence. Resulting values are Hamming distances,
# i.e. the total number of ties that are different. Also by default: self-ties are ignored (diag=FALSE)
# SE_dist <- 1 - (sedist(net,method = "correlation")/2)
SE_dist <- sedist(net,method = "hamming")

# Check it out
SE_dist

# How to interpret this matrix? In each cell you see the number
# of non-shared contact persons.

# The higher this number, the less
# contact persons they share, and the less structurally
# equivalent the two individuals are.
# In other words: do both know Bono? Great! Do both NOT know Bono? Great!
# All other cases: add 1 to the "penalty" you see here

# For Zache: For instance, lets compare actor 5 and 7 (sedist of 6):
# Both 5 and 7 are connected to 1: great!
# They are also both disconnected to 28 of the alters: great!
# However, 5 is connected to 11; 7 isn't.
# And 7 is connected to 6 and 17, which 5 isn't.
# So a total "penalty" of 3, which is doubled to 6 as this includes both rows and columns

# Okay, back to function "equiv.clust"...
# This function performs a hierarchical clustering analysis:

#net_clustering <- equiv.clust(net,equiv.fun="sedist",cluster.method='complete')

net_clustering <- equiv.clust(net,equiv.fun="sedist",method="hamming",cluster.method='complete')

# You can examine this data object, but you wont see as much detail here:
net_clustering

# But it is far easier to plot the dendrogram of this clustering structure
plot(net_clustering)
# For zache, you note that 33 and 34 are very equivalent, 1 is unique and then there
# are additional sets of students that are very equivalent

# To create a blockmodel, we need a network and a partition
# We use the blockmodel function in the sna package to create blockmodels,
# providing it with the network, the clustering, and number of subsets (aka positions/clusters)
# When working with the zache, I here decided to try two partitions with, respectively
# 3 and 7 positions

se.blockmodel.A <- blockmodel(net, net_clustering, k=2)
se.blockmodel.B <- blockmodel(net, net_clustering, k=5)

# However, whereas this input data was useful for deriving structurally equivalent positions
# ...we want to display the original data in the blockmodels, i.e. what is in 'mat'

# For matching block membership information with the original
# data set, unfortunately some effort needs to be done: 
se.position.A <- 1+se.blockmodel.A$block.membership[match(1:nrow(mat),se.blockmodel.A$order.vector)]
se.position.B <- 1+se.blockmodel.B$block.membership[match(1:nrow(mat),se.blockmodel.B$order.vector)]



# plot the network with the block colouring:
gplot(net, vertex.col=se.position.A, arrowhead.cex = FALSE)
gplot(net, vertex.col=se.position.B, arrowhead.cex = FALSE)
# you can add option "displaylabels=TRUE" if you want to

# For zache: verify that the centres are indeed Mr. Hi (node 1) and John A (node 34)
gplot(net, vertex.col=se.position.A,displaylabels = TRUE, arrowhead.cex = FALSE)
gplot(net, vertex.col=se.position.B,displaylabels = TRUE, arrowhead.cex = FALSE)

# finally, here a plot of the block matrices for the two
# solutions we extracted:
plot(se.blockmodel.A)
# For Zache: the first block consists of John A. and a friend, the second
# block of Mr. Hi and a friend. The third block consists of a
# remaining periphery of club members.
plot(se.blockmodel.B)
# For Zache: the core is now split into single actors - notable here:
#  the periphery is split into two sets of club
# members that - surprising or not - pretty much coincides
# with how the club broke up into two separate clubs later.

###############################
# REGULAR EQUIVALENCE ANALYSIS
###############################

?redist
# The function "redist" indicates approximatively dissimilarity
# in the sense of "absence of regular equivalence" of actors.
# Details are not needed here, can be read in the paper by
# Borgatti & Everett referenced on the R help page.

# We proceed analogous to the structural equivalence analysis above.

# BUT we use different data set - this seems not to work with Zachary's!

# Main reason: the heuristic for determining regular equivalence in R only works for asymmetric data

# With the sna package comes a set of example datasets
# these can be loaded into your workspace with the data command
# So lets load "coleman"

data(coleman)
interact_Coleman <- coleman[2,,]
# Spring 1958 network from James Coleman's "Introduction to
# Mathematical Sociology" boy school cohort.
?coleman
# As you can see, our interest here is on friendship among the 73 boys in the spring of 1958
# This network is directed (reported friendships) and binary

# Plot it (keeping arrowheads, as it is directed)
gplot(interact_Coleman, displaylabels = TRUE)

# calculate approximate distances acc. to regular equivalence:
interact.RE.distance <- redist(interact_Coleman) 
interact.RE.distance
?redist
# Back to function "equiv.clust"...
# This function performs a hierarchical clustering analysis:
interact.clustering <- equiv.clust(interact_Coleman, equiv.fun="redist", method='catrege')
# Let us plot the resulting dendrogram:
plot(interact.clustering)
# Potentially interesting split is into 
# seven groups (at height 0.8)...

# Now get back to function "blockmodel"...
# If we want to extract seven communities, do this:
re.blockmodel.7 <- blockmodel(interact_Coleman, interact.clustering,k=7)

# for matching block membership information with the original
# data set, unfortunately some effort needs to be done: 
re.position.7 <- 1+re.blockmodel.7$block.membership[
	match(1:nrow(coleman[2,,]),re.blockmodel.7$order.vector)]

# plot the network with the block colouring:
gplot(interact_Coleman, vertex.col=re.position.7, displaylabels=T)
# and without arrowheads if that is easier (though it is directed):
gplot(interact_Coleman, vertex.col=re.position.7, displaylabels=T,arrowhead.cex = F)

# >> there is one huge group of red nodes ("average students")
# >> there are four blue complete isolates
# >> there are six green peripherals who try to attach to
#    others but are not accepted / not nominated by others
# >> pink/purple and cyan nodes seem to be boundary spanners
#    in bridge positions, not part of a single cohesive group
# >> yellow nodes seem to be in the center of groups, maybe
#    their leaders?
# >> there is one more grey node who seems popular but does
#    not nominate others

# a plot of the block matrices for the 7 group solutions:
plot(re.blockmodel.7)
# Regular and null blocks! This can be reduced to an image graph

# Would be prudent to also try different number of positions (!=7)