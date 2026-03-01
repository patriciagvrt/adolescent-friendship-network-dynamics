##############################################
# R script illustrating the modelling of
# network-behaviour co-evolution in RSiena.
# version of 4 July 2022
# written by Christian Steglich
##############################################

# The data used in this script were collected
# by Vanina Torlò as part of her PhD study.
# A comprehensive description and analysis can
# be found in: Alessandro Lomi, Tom Snijders, Christian
# Steglich & Vanina Torlò (2011). Why are some more peer
# than others? Evidence from a longitudinal study of
# social networks and individual academic performance.
# Social Science Research 40(6), 1506-1520.

# load data set into the workspace of this R session:
load("MBAData.RData")
ls()

# load SAOM-estimation software into namespace:
library(RSiena)


# Prepare data in "RSiena" format:
# DEPENDENT NETWORK:
advice <- sienaDependent(array(c(advice.net1, advice.net2, advice.net3),
	dim=c(75,75,3)))
# PREDICTOR NETWORK (time-varying):
friendship <- varDyadCovar(array(c(friendship.net1, friendship.net2),
	dim=c(75,75,2)))
# CONSTANT ACTOR ATTRIBUTES:
gender <- coCovar(covariates$sex)
bachelor <- coCovar(covariates$degree)
# DEPENDENT ACTOR ATTRIBUTE:
performance <- sienaDependent(performance.mat, type='behavior')
# JOIN THEM ALL INTO DATA OBJECT FOR SAOM:
coevdata <- sienaDataCreate(advice, friendship, gender, bachelor, performance)


# Specify co-evolution model:
coevmodel <- getEffects(coevdata)

# ADD EFFECTS TO NETWORK OBJECTIVE FUNCTION:
coevmodel <- includeEffects(coevmodel, name='advice',
	transTrip, transRecTrip) # clustering
coevmodel <- includeEffects(coevmodel, name='advice',
	X, interaction1='friendship') # friendship effect
coevmodel <- includeEffects(coevmodel, name='advice',
	egoX, altX, sameX, interaction1='gender') # three gender effects
coevmodel <- includeEffects(coevmodel, name='advice',
	sameX,interaction1='bachelor') # educational background homophily
coevmodel <- includeEffects(coevmodel, name='advice',
	egoX,altX,simX,interaction1='performance') # three performance effects

# ADD EFFECTS TO BEHAVIOUR OBJECTIVE FUNCTION:
coevmodel <- includeEffects(coevmodel, name='performance',
	indeg, outdeg, avSim, interaction1='advice') # three effects conjugate to the above
coevmodel <- includeEffects(coevmodel, name='performance',
	effFrom, interaction1='gender') # effect of gender on performance
# INSPECT MODEL SPECIFICATION:
coevmodel


# Estimate it and inspect results:
controls <- sienaAlgorithmCreate(seed=123456)
(coevresults <- siena07(controls, data=coevdata, effects=coevmodel))

# make HTML table:
siena.table(coevresults, type="html", tstat=TRUE, sig=TRUE, d=2)
browseURL("coevresults.html")



# INTERPRETATION: (possible because of convergence)
# ===============
# - Advice seeking model akin to "network only" results from script 1
# - Performance effects on advice seeking more pronounced here
# - There is peer influence: students' performance moves in the
#   direction of their peers' performance (avSim effect)
# - Students who get asked for advice by more fellow students perform
#   better than students who get asked by few fellow students.
