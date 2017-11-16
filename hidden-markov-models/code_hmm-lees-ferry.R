# Title: using Hidden Markov models to represent the persistence of flow in the Nile River
# Background: the Nile data used below is a subset of the larger 800+ year record in which hydrologist 
#          H.E. Hurst noticed persistence in streamflow records heavier than average flood years while 
#          below average flood years were typically followed by light flood years.
# Source: http://blog.revolutionanalytics.com/2014/09/intro-to-long-memory-herodotus-hurst-and-h.html

# clear environment
rm(list=ls()) 

# set current working directory
setwd("C:/Users/Billy/Google Drive/Teaching/github/Advanced Data Analysis/hidden-markov-models")

# turn off warnings from packages for printing with knitr
options(warn=-1)

# import packages and user-defined Hidden Markov fucntions
suppressMessages(library(datasets))     # Nile River data
suppressMessages(library(forecast))     # time series analysis
suppressMessages(library(tidyverse))    # set of useful packages: ggplot2, dplyr, etc.
suppressMessages(library(ggfortify))    # visualization for probability distributions
suppressMessages(library(gridExtra))    # arrange multiple ggplot figures on same plot
suppressMessages(library(MASS))         # fitting probability distributions
suppressMessages(library(HiddenMarkov)) # Hidden Markov models
theme_set(theme_gray())   # set theme for ggplot2

# turn warnings back on
options(warn=0)

# import user-defined Hidden Markov functions
Pi_init <- function(n){
  # purpose: set prior for transition matrix
  # input: n = Markov chain order (integer)
  
  matrix(rep(1/n,n^2),n)   
}

delta_init <- function(n){
  # purpose: set prior for marginal distribution of Markov chain
  # input: n = Markov chain order (integer)
  
  d <- rnorm(n)^2
  d/sum(d)
}

# import Lee's Ferry monthly Flow data and transform into annual time series
lees <- read_delim(file="./data/leesferry-mon-data.txt", delim=" ")

lees <- lees %>% gather(month, flow, -year) %>% 
  group_by(year) %>%
  summarize(annual_flow = mean(flow)/10^6)  # annual flow (million acre-feet)

x <- ts(data=lees$annual_flow, start=lees$year[1], end=lees$year[1]+nrow(lees))

# plot annual streamflow at Lee's Ferry
autoplot(x) + 
  ggtitle("Annual flow at Lee's Ferry") +
  ylab("Annual Flow (MAF)") + 
  xlab("Year")

# select number of states
order <- 2 

# set priors for transition matrix (Pi) and marginal distribution of Markov chain (delta)
Pi <- Pi_init(order)  # generate uninformed prior for the transition matrix
delta <- delta_init(order)  # generate prior marginal distribution of Markov chain

# fit lognormal distribution (positive values, fat-tail)
distn.type<- "lognormal"
fit <- fitdistr(x, distn.type)
pm <- list(meanlog = rep(log(mean(x)),order), sdlog = rep(log(sd(x)), order))
distn <- "lnorm"  # string in base R to represent lognormal distribution

# initialize model
init.hmm <- dthmm(x, Pi, delta, distn, pm) 

# fit model using the Baum-Welch algorithm
set.seed(101)  # set pseudorandom seed for reproducibility
fit.hmm <- BaumWelch(init.hmm, bwcontrol(maxiter = 1000,posdiff=TRUE,prt=FALSE))

# ERROR here: could not figure out what was wrong here. 


