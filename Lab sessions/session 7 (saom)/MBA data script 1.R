##############################################
# R script illustrating the modelling of
# network dynamics in RSiena.
# version of 8 March 2024
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

# quick overview of actor attribute coding:
head(covariates)
# sex    : 1=female, 2=male
# GPA    : university entry test result
# age    : in years upon study start
# degree : bachelor degree (qualitative variable)
# nat    : nationality (1=Italian, 2=Non-Italian)
head(performance.mat)
# average grade per time point on Italian scale (20-30)


# define function for quantifying change:
Hamming <- function(net1,net2) {
	tbl <- table(c(0,0,1,1,net1),c(0,1,0,1,net2))-1
	return(tbl[1,2]+tbl[2,1])
}
# function expects two binary networks and
# returns their Hamming distance

# How much change is there between observation moments?
Hamming(advice.net1,advice.net2) 
Hamming(advice.net2,advice.net3)
# Quite a bit of change, there should be enough power to detect
# "rules of change" if the model specification is not too big.


# define function for quantifying stability:
Jaccard <- function(net1,net2) {
	tbl <- table(c(0,0,1,1,net1),c(0,1,0,1,net2))-1
	return(tbl[2,2]/(tbl[1,2]+tbl[2,1]+tbl[2,2]))
}
# function expects two binary networks and
# returns their Jaccard stability index

# How much stability is there between observation moments?
Jaccard(advice.net1,advice.net2) 
Jaccard(advice.net2,advice.net3) 
# Both above the manual's rule-of-thumb ">=0.3" criterion;
# enough stability for assessing also the structural network
# effects which condition on the existence of a small subgraph


# load SAOM-estimation software into namespace:
library(RSiena)
packageVersion("RSiena")
# (if not yet installed, please do)

# Prepare data in "RSiena" format:
# DEPENDENT NETWORK:
advice <- sienaDependent(array(c(advice.net1, advice.net2, advice.net3
),
	dim=c(75,75,3)))
# PREDICTOR NETWORK (time-varying):
friendship <- varDyadCovar(array(c(friendship.net1, friendship.net2),
	dim=c(75,75,2)))
# CONSTANT ACTOR ATTRIBUTES:
gender <- coCovar(covariates$sex)
bachelor <- coCovar(covariates$degree)
# TIME-VARYING ACTOR ATTRIBUTE:
performance <- varCovar(performance.mat[,1:2])
# JOIN THEM ALL INTO DATA OBJECT FOR SAOM:
thedata <- sienaDataCreate(advice, friendship, gender, bachelor, performance)


# Specify model for the dynamics of advice seeking:
(themodel <- getEffects(thedata)) # default: rate, outdegree, reciprocity
# ADD EFFECTS TO OBJECTIVE FUNCTION:
themodel <- includeEffects(themodel,
	transTrip, transRecTrip) # clustering
themodel <- includeEffects(themodel,
	X, interaction1='friendship') # friendship main effect
themodel <- includeEffects(themodel,
	egoX,altX, # sender, receiver, and ..
	sameX,interaction1='gender') # .. homophily effects on gender
themodel <- includeEffects(themodel,
	sameX, interaction1='bachelor') # educational background homophily
themodel <- includeEffects(themodel,
	egoX, altX, # sender, receiver, and ..
	simX, interaction1='performance') # .. homophily effects on performance
# INSPECT MODEL SPECIFICATION:
themodel


# Estimate it and inspect results:
controls <- sienaAlgorithmCreate(seed=12345)
# seed=12345 : safeguards we all get same randon numbers & hence results
(theresults <- siena07(controls, data=thedata, effects=themodel))

# make HTML table:
siena.table(theresults, type="html", tstat=TRUE, sig=TRUE, d=2)
browseURL("theresults.html")
# first note that the algorithm has converged according to the two
# criteria that the RSiena manual suggests:
# (1) parameterwise convergence indicator below 0.1 in absolute terms
# (2) max. convergence over all linear combinations of parameters below 0.25

# Interpretation:
# ===============
# - rate of change goes down a bit
# - students seek advice from few other students
# - they exchange advice reciprocally
# - they form transitive (hierarchically ordered) triples
# - reciprocation is weaker when embedded in such a triple 
# - they ask advice from friends
# - males ask fewer others for advice than women
# - there homophily-selection on gender, BA-degree and performance
# - sender and receiver effects of performance in line with hierarchy

	# Exercise: (can be skipped)
	# =========
	# Are hierarchy-related results an artifact of having omitted endogenous
	# effects related to the degree distributions, such as the Matthew effect
	# (Merton 1968)?

	# Model answer:
	# =============
	# Enrich model by endogenous, degree-related terms:
	model2 <- includeEffects(themodel,name='advice',inPop,outAct,inAct) # degree effects
	(results2 <- siena07(controls, data=thedata, effects=model2))

	# again make HTML table:
	siena.table(results2, type="html", tstat=TRUE, sig=TRUE, d=2)
	browseURL("results2.html")

	# Interpretation:
	# ===============
	# - indeed there is evidence for the Matthew effect (indegree popularity)
	# - students also ask fewer others for advice when they get asked by many
	# - interaction reciprocation × transitivity falls from significance
	# - the sender and receiver effects of performance are substantially weaker

# Now let's move towards goodness-of-fit inspection...

# specify the options for the simulation algorithm:
simcontrols <- sienaAlgorithmCreate(n3=250,nsub=0,seed=12345)
# nsub=0 : phase 2 (parameter estimation) will be skipped
# n3=250 : in phase 3, sample 250 independent network evolution processes

# add simulated networks to model results:
thesims <- siena07(simcontrols, data=thedata, effects=themodel,
	returnDeps=TRUE, # simulated networks at ends of observation period made available
	prevAns=theresults) # use estimates from "theresults" for simulations

# goodness of fit for indegree distribution:
(gofIndegrees <- sienaGOF(thesims, varName="advice",
	IndegreeDistribution, # the fit function to be evaluated
	cumulative=FALSE, # (same meaning as in probability distributions)
	levls=0:15)) # evaluate indegrees on these 16 fit dimensions
plot(gofIndegrees) # looks a bit bad at indegrees zero and five
# But p-value is in non-significant region -> fit acceptable

# goodness of fit for outdegree distribution:
(gofOutdegrees <- sienaGOF(thesims, varName="advice",
	OutdegreeDistribution, cumulative=FALSE, levls=0:15)) 
plot(gofOutdegrees) # looks good, and p is really non-significant.

# goodness of fit for triad census:
(gofTriads <- sienaGOF(thesims, varName="advice", TriadCensus, verbose=TRUE))
plot(gofTriads,center=TRUE,scale=TRUE) # looks bad, and p=0
# (exactly zero, this is a Monte-Carlo p-value for 250 simulated data sets)

# Inspect this a bit closer for possible model improvement:
# (1) triad 021U is not simulated often enough,
# (2) triad 111U is simulated too often,
# and two other, less pronounced dimensions: 021C and 111D.

# Conclusion:
# ===========
# Based on the GOF-assessment above, only the triad census gives
# reasons for concern. To improve it, one could include
# (1) effect inPop (indegree popularity = Matthew effect) to improve
#     triad 021U; a positive effect would generate more such triads;
# (2) effect reciAct (reciprocated degree activity) to improve
#     triad 111U; a negative effect would generate less such triads.

	# Exercise: (can be skipped)
	# =========
	# Do the above and check if GOF indeed improves by doing so.

	# Model answer:
	# =============
	# Enrich model as suggested by triad census misfit:
	model3 <- includeEffects(themodel,name='advice',inPop,reciAct)
	(results3 <- siena07(controls, data=thedata, effects=model3))

	# again make HTML table:
	siena.table(results3, type="html", tstat=TRUE, sig=TRUE, d=2)
	browseURL("results3.html")

	# Interpretation of newly added effects:
	# ======================================
	# - estimates of new effects have the expexted signs and are significant
	# - the positive Matthew effect (indegree popularity) indicates that
	#   popular advisors attract additional advice seeking
	# - the negative reciAct effect indicates that students with more
	#   reciprocated advisors tend to have less unreciprocated advisors

	# GOF-assessmemt:
	# ===============
	# Simulate data from this model:
	sims3 <- siena07(simcontrols, data=thedata, effects=model3,
		returnDeps=TRUE, prevAns=results3)

	# goodness of fit for indegree distribution:
	(gofIndegrees3 <- sienaGOF(sims3, varName="advice",
		IndegreeDistribution, cumulative=FALSE, levls=0:15))
	plot(gofIndegrees3) # was okayish fit, is good now

	# goodness of fit for outdegree distribution:
	(gofOutdegrees3 <- sienaGOF(sims3, varName="advice",
		OutdegreeDistribution, cumulative=FALSE, levls=0:15)) 
	plot(gofOutdegrees3) # was good fit, now even better

	# goodness of fit for triad census:
	(gofTriads3 <- sienaGOF(sims3, varName="advice", TriadCensus, verbose=TRUE))
	plot(gofTriads3,center=TRUE,scale=TRUE) 
	# was really bad, now really good 
