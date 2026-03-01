# ===================================================
# This R script illustrates the use of the package
# 'ergm' for the estimation of exponential random
# graph models, simulation of artificial data, and 
# the issue of degenerate model specifications.
# version of 4 March 2023
# written by Christian Steglich 
# ===================================================

# NOTE that the Statnet-team offers a tutorial script
# at https://statnet.org/workshop-ergm
browseURL("https://statnet.org/workshop-ergm")

# Overview: 
# =========
# a) Installation of new packages
# b) ERGM analysis of a directed network (s50)


# a) Installation of new packages
# ===============================

install.packages("latticeExtra") # for next lab only
install.packages("ergm")

# Load the vocabulary of the neccessary packages:
library(network)		# network data storage
library(sna)			# network analysis routines
library(latticeExtra)	# for nicer convergence & gof plots
library(ergm)			# fitting & evaluating ERGMs

# Always check and report software version:
packageVersion("ergm")
# this script was last updated for version 4.1.2


# b) ERGM analysis directed data (s50)
# =========================================
 
# The examples make use of the s50 data set
browseURL("https://www.stats.ox.ac.uk/~snijders/siena/s50_data.htm")
load("s50.RData")

# We will analyse the first of these networks, so let us
# start by preparing it as a network object:

friendnet <- network(s501)
friendnet %v% "alcohol" <- s50a[,1]
friendnet %v% "alcohol.colours" <- 
	c("white","yellow","orange","red","darkred")[s50a[,1]]
friendnet %v% "tobacco" <- s50s[,1]

# Visual inspection always is a good first step:
plot(friendnet,vertex.cex="tobacco",vertex.col="alcohol.colours")
# What do we see?
# visible endogenous dependencies:
# > triangles (clustering tendencies)
# > reciprocation
# > some isolates
# > are there hubs? / does the Matthew effect operate? 
#   (let us check this by testing)
# visible exogenous dependencies:
# > colour (=alcohol) homophily
# > size (=tobacco) homophily
# > maybe larger/darker nodes receive more ties?
#   (let us check this by testing, below)

# Let us inspect the help pages for model estimation...
?ergm
# .. and model specification:
?'ergm-terms'

# Begin with a Bernoulli random graph model / Erdos-Renyi model:
model1 <- friendnet ~ edges
summary(model1)
# The summary method gives the "target statistics", here: the tie count.

# Let's estimate this model:
results1 <- ergm(model1)
# Note that MPLE is reported as final estimate.

# Inspect the results:
summary(results1)
# Not very spectacular yet.

# What kind of networks does this model generate?

# Simulate 100 networks (seed ensures we all get the same):
?simulate.ergm
sims1 <- simulate(results1, nsim = 100, seed = 1234)

# Plot the first one:
plot(sims1[[1]], vertex.cex = "tobacco",
	vertex.col = "alcohol.colours")
# This is a truly random network, all the regularities observed
# for the empirical network above do NOT hold here.

# But the number of edges (and hence density) should be met, on average, right?

# Make histogram of network density of 100 simulated data sets:
hist(gden(sims1))

# Or slightly more beautiful: a density plot
plot(density(gden(sims1)), main = 'network density')

# Add vertical line where empirical data are:
lines(x = rep(gden(friendnet), 2),
	y = c(0, max(density(gden(sims1))$y)),
	col = "red", lwd = 2
)
# looks okay ("converged").

# Make another plot of transitivity index distribution:
plot(density(gtrans(sims1)), xlim=c(0, 0.4), main = 'transitivity index')
lines(x = rep(gtrans(friendnet), 2),
	y = c(0, max(density(gtrans(sims1))$y)),
	col="red", lwd=2
)
# Not well-modelled: simulated networks have too little transitivity!
# In fact, transitivity in simulated networks is in expectation
# the density of the network:
lines(x = rep(gden(friendnet), 2),
	y = c(0, max(density(gtrans(sims1))$y)),
	col = "blue", lwd = 2, lty = "dashed"
)

# Now as we do know that tobacco and alcohol homophily matter
# (from the visualisation), let us add the corresponding effect
# terms to the model. (Note that this still retains the assumption
# of tie independence!)
model2 <- friendnet ~ edges + absdiff("alcohol") + absdiff("tobacco")
summary(model2)

results2 <- ergm(model2)
summary(results2)
# It seems tobacco homophily does not help understanding the data!
# Also the information criterion BIC got worse instead of better, 
# which indicates that model size might be reduced. So, drop tobacco:

model2a <- friendnet ~ edges + absdiff("alcohol")
results2a <- ergm(model2a)
summary(results2a)
# Indeed now AIC and BIC are smaller than for model1.

# Again, let us see if this delivers networks that "look more real":
sims2a <- simulate(results2a, nsim = 100, seed = 1234)
# Plot the first one:
plot(sims2a[[1]], vertex.cex="tobacco", vertex.col="alcohol.colours")
# Well, colours are more separated now - as could be expected.
# Everything else is still very random.


# UNDERSTANDING DEGENERACY: LINEAR TRANSITIVITY SPECIFICATION:
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# 20th century specifications of ERGMs often use "maximum pseudo-
# likelihood estimation", i.e., logistic regression under a
# tie independence assumption. In this tradition, the ERGM was still
# called the p-star-model (p*). But there were inherent problems
# with this approach. It was known already at the beginning that
# the MPL-estimation assumes away network dependence and therefore
# would be prone to give faulty, biased results.

# It turned out to be even more serious: Because logistic
# regression is so reliably producing results, the MPLE glossed over
# serious problems with model specifications everyone thought would
# work well, dating back to Ove Frank & David Strauss' seminal
# paper "Markov Graphs" (1986, Journal of the American Statistical
# Association). This became known as "model degeneracy".

# The models above (1, 2 and 2a) all did not contain terms that
# instantiate tie interdependence. Because of this, the MPL-estimate
# was identical to the proper ML-estimate. But now, we add
# reciprocity and transitivity to the model, and include "twopath"
# as stabilising term for estimating transitivity. Now, MPLE
# and (MC)MLE will differ, and accordingly also only now, proper
# simulation-based inference will be used. 

# FIRST ATTEMPT: linear transitivity specification
model3 <- friendnet ~ edges + absdiff("alcohol") +
	mutual + twopath + ttriple
summary(model3)

results3 <- ergm(model3, control = control.ergm(seed = 1234))
# fails after 3 iterations and signals estimation problems!

# CLOSER INSPECTION: estimate with logistic regression instead
# (this is wrong, as it assumes independence of tie variables -
# but will allow closer inspection):
results3wrong <- ergm(model3, estimate="MPLE")
summary(results3wrong)
# note that AIC and BIC are pseudo-likelihood-based and hence
# faulty, in particular they are incomparable to other models.

# Let us eyeball the model's simulations again:
sims3wrong <- simulate(results3wrong, nsim = 100, seed = 1234)

# Make again a plot of the transitivity index distribution:
plot(density(gtrans(sims3wrong), from = 0, to = 1),
	xlim = c(0, 1), main = 'transitivity index')
lines(x = rep(gtrans(friendnet), 2),
	y = c(0, max(density(gtrans(sims3wrong))$y)),
	col="red", lwd=2
)
# Not well-modelled: simulated networks have way too much transitivity!
# -> not converged by simulation-based inference standards
# Here this is actually caused by a DEGENERATE MODEL SPECIFICATION.

# Plot the first two plus the empirical network:
par(mfrow = c(1, 3)) # makes a screen with three plot areas horizontally ordered
plot(friendnet, vertex.cex="tobacco",
	vertex.col = "alcohol.colours", main = "empirical")
plot(sims3wrong[[1]], vertex.cex = "tobacco",
	vertex.col = "alcohol.colours", main = "sim1")
plot(sims3wrong[[2]], vertex.cex = "tobacco",
	vertex.col = "alcohol.colours", main = "sim2")
# We see what is going wrong here: "Clusters cluster too strongly."
dev.off()


# NON-DEGENERATE TRANSITIVITY SPECIFICATION - HOW TO DO IT SINCE 2006:
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# Now, let us pick a different transitivity effect that
# allows for proper (ML) estimation:
model3good <- friendnet ~ edges + absdiff("alcohol") +
	mutual + twopath + gwesp(0.5,fixed=TRUE)
results3good <- ergm(model3good, control = control.ergm(seed = 1234))
# gives results! Which ones?
summary(results3good)
# Much better  according to information criteria.

# Let us eyeball the model's simulations again:
sims3good <- simulate(results3good, nsim = 100, seed = 1234)

# Make again a plot of the transitivity index distribution:
plot(density(gtrans(sims3good), from = 0, to = 1),
	xlim = c(0, 1), main = 'transitivity index')
lines(x = rep(gtrans(friendnet), 2),
	y = c(0,max(density(gtrans(sims3good))$y)),
	col="red", lwd=2
)
# Finally quite well-modelled! Neat..

# Plot the first two plus the empirical network:
par(mfrow = c(1, 3)) # makes a screen with three plot areas horizontally ordered
plot(friendnet, vertex.cex = "tobacco", 
	vertex.col = "alcohol.colours", main = "empirical")
plot(sims3good[[1]], vertex.cex = "tobacco", 
	vertex.col = "alcohol.colours", main = "sim1")
plot(sims3good[[2]], vertex.cex = "tobacco",
	vertex.col = "alcohol.colours", main = "sim2")
# not entirely bad, this time.
dev.off()


# Final model specification of the day:
# > Add "gwidegree" effect to assess evidence for Matthew effect 
# > Add "nodeicov" effect to assess if drinkers receive more ties

model4 <- friendnet ~ edges + absdiff("alcohol") +
	mutual + twopath + gwesp(0.5,fixed=TRUE) +
	gwidegree(0.5,fixed=TRUE)+nodeicov("alcohol")
summary(model4)
results4 <- ergm(model4, control = control.ergm(seed = 1234))
summary(results4)
# Not better according to the information criteria.
# Results show weak evidence for drinkers being attractive as friends,
# and no evidence for a Matthew effect in friendship.


# NEXT LAB: convergence assessment, goodness of fit, effect sizes.
