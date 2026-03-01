##############################################
# R script for the course "Social Network
# Analysis" illustrating ERGM use,
# version of 21 February 2024
# written by Christian Steglich
##############################################

# The IFM data used in this lab exercise were
# collected by David de Jong at the Institute
# for Finance Management (IFM) in Dar Es Salaam,
# Tanzania, as part of his sociology master
# thesis (2011, University of Groningen).

##############################################

# load data set into the workspace of this R session:
load("IFM.RData")
ls()

# inspect attribute variables:
head(IFM.attributes)

# load required namespace:
library(ergm)

# We formulate a list of ad-hoc hypotheses for
# comprehensive testing on the basis of ERGM
# analysis.

# Hypothesis 1a,b ('tenure and job status receiver')
# ^^^^^^^^^^^^^^^
# IFM employees with longer tenure or higher job
# status receive more advice seeking nominations.

# Hypothesis 2a,b ('tenure and job status sender')
# ^^^^^^^^^^^^^^^
# IFM employees with longer tenure or higher job
# status send fewer advice seeking nominations.

# Hypothesis 3a,b,c ('tenure, job title and gender homophily')
# ^^^^^^^^^^^^
# IFM employees ask each other more for advice when
# they have similar tenure, the same job title or the
# same gender.

# Hypothesis 4 ('preferential attachment')
# ^^^^^^^^^^^^
# IFM employees are asked for advice based on whether
# others also ask them for advice.

# Hypothesis 5 ('transitive closure')
# ^^^^^^^^^^^^
# IFM employees ask the advisors of their advisors
# directly for advice.

# Hypothesis 6 ('friendship closure')
# ^^^^^^^^^^^^
# IFM employees ask the friends of their friends
# for advice.
# >> control for direct friendship, please!

# Prepare data in "network/sna" format:
advice <- network(IFM.advice, directed=TRUE)
advice %v% "sex" <- as.character(IFM.attributes$sex)
advice %v% "tenure" <- IFM.attributes$tenure
advice %v% "jobtitle" <- as.character(IFM.attributes$jobtitle)

# For Hypotheses 1+2 make a numerical status variable
# from the jobtitles:
IFM.attributes$status <- c(2, 3, 4, 1)[as.numeric(as.factor(
	IFM.attributes$jobtitle
))]
# Check how this goes together:
table(IFM.attributes$jobtitle, IFM.attributes$status)
# From low to high status:
# tutor >> assistant >> lecturer >> senior lecturer

advice %v% "status" <- IFM.attributes$status

# For Hypothesis 6 make 'friend of friend' matrix:
friendoffriend <- IFM.friendship %*% IFM.friendship

# Specify model that facilitates all the hypothesis tests:
my.model <- advice ~ edges +		# Intercept
	nodeicov("tenure") +		# Hypothesis 1a +
	nodeicov("status") +		# Hypothesis 1b +
	nodeocov("tenure") +		# Hypothesis 2a -
	nodeocov("status") +		# Hypothesis 2b -
	absdiff("tenure") + 		# Hypothesis 3a -
	nodematch("jobtitle") + 	# Hypothesis 3b +
	nodematch("sex") + 			# Hypothesis 3c +
	mutual +					# Reciprocity effect ?
	edgecov(IFM.friendship) +	# main effect of friendship
	edgecov(friendoffriend) +	# Hypothesis 6 +
	istar(2) +					# Hypothesis 4 
	twopath +					# for stabilisation of model* -
	gwesp(0.5, fixed=TRUE)		# Hypothesis 5 +
# * the twopath effect can also be interpreted as a test
#   for the hypothesis
#   "Are people being asked for advice who themselves
#    ask other people for advice?"

# Describe the (global) advice seeking network in terms of
# how many (local) subgraphs according to the above model
# it contains:
summary(my.model)

# Estimate the model (with all of us having same
# "random" numbers):
my.results <- ergm(my.model, control=control.ergm(888))
# This can take much longer than in our example here!

# FIRST, did the model actually converge?

# MCMC diagnostics show how much the sample of network
# statistics that were simulated deviate from the observed
# network statistics (calculated in code line 95 above):
mcmc.diagnostics(my.results)
# This looks 'reasonably converged' (centered around zero).

# SECOND we can inspect the results:
summary(my.results)

# Hypothesis tests:
# ^^^^^^^^^^^^^^^^^
# H1a (tenure receiver) 	non-significant
# H1b (status receiver) 	supported

# H2a (tenure sender)		non-significant
# H2b (status sender)		non-significant

# H3a (homophily tenure)	supported
# H3b (homophily status)	non-significant
# H3c (homophily gender)	weakly supported

# H4 (hub formation)		non-significant
# H5 (transitive closure)	supported
# H6 (friendship closure)	non-significant

# Other effects:
# ^^^^^^^^^^^^^^
# Reciprocity	happens in advice seeking
# Two-paths		are avoided in advice seeking
# Friendship	facilitates advice seeking


# THIRD let us look at the goodness of fit on dimensions
# other than those included in the model specification.
# Does the model reproduce some global features of the network?

# We check whether the empirical distribution of geodesic distances
# and shared partner counts can be explained by the model:
my.fit <- gof(my.results,
	control = control.gof.ergm(seed = 888, nsim=200))
plot(my.fit)
# This looks like good fit because the empirical data (black
# line) are located well within the cloud of simulated data
# (boxplots, central 95% intervals) on these additional
# dimensions.
# NOTE that model dimensions look not so very good, main reason is
# probably small sample size (nsim=200 is default in the controls).
print(my.fit)
# adds Monte Carlo p-values per box plotted - ideally, they would
# all be non-significant: these values indicate the fraction of
# simulated data sets that are farther away from their own average
# than the observed data are.