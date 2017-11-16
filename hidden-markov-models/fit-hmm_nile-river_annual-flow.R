# Title: using Hidden Markov models to represent the persistence of flow in the Nile River
# Background: the Nile data used below is a subset of the larger 800+ year record in which hydrologist 
#          H.E. Hurst noticed persistence in streamflow records heavier than average flood years while 
#          below average flood years were typically followed by light flood years.
# Source: http://blog.revolutionanalytics.com/2014/09/intro-to-long-memory-herodotus-hurst-and-h.html

# clear environment
rm(list=ls()) 

# set current working directory
setwd("C:/Users/Billy/Google Drive/Teaching/HMM/github")

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

# Nile River flow data
x <- Nile

# plot annual streamflow for the Nile River 
autoplot(x) + 
  ggtitle("Measurements of the annual flow of the Nile River at Aswan") +
  ylab("Annual Flow (10^8 m^3)") + 
  xlab("Year") + 
  geom_vline(xintercept=1898, linetype=4)  # show changepoint

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

# BUG in HiddenMarkov (I have contacted the developer and he will fix it)
fit.hmm$delta <- c(0,1)

# decode the model using the Viterbi algorithm
vit.hmm <- Viterbi(fit.hmm)  

# create a data frame for plotting
x.df <- tibble(year=time(x) %>% as.integer(), flow=x %>% as.double(), state=vit.hmm)

# plot results showing in what state each point belongs (p1)
p1 <- ggplot(aes(x=year, y=flow), data=x.df) + 
  ggtitle("Measurements of the annual flow of the Nile River at Aswan") +
  ylab("Annual Flow (10^8 m^3)") + 
  xlab("Year") + 
  geom_vline(xintercept=1898, linetype=4) + # highlight changepoint
  geom_line() +
  geom_point(aes(color=factor(state))) +  # show state with color
  labs(color="State")  # add title for color legend

# plot fitted distributions (p2)
plot.range <- seq(0, 2000, 1)
n.sim <- 10000

## simulate data from state #1
set.seed(101)  # set pseudorandom seed for reproducibility
sim.state1 <- rlnorm(n.sim, fit.hmm$pm$meanlog[1], fit.hmm$pm$sdlog[1])
df.state1 <- data.frame(flow = sim.state1, state = 1)

## simulate data from state #2
set.seed(101)  # set pseudorandom seed for reproducibility
sim.state2 <- rlnorm(n.sim, fit.hmm$pm$meanlog[2], fit.hmm$pm$sdlog[2])
df.state2 <- data.frame(flow = sim.state2, state = 2)

df.states <- rbind(df.state1, df.state2)

p2 <- ggplot(df.states) +
  geom_density(aes(x = flow, color = factor(state))) +
  theme(legend.position = "none") + # remove legend for this plot
  ggtitle("Fitted lognormal distributions for each state") +
  xlab("Annual Flow (10^8 m^3)") +
  ylab("Density")

# plot time series along with fitted distributions
grid.arrange(p1, p2, ncol=1)

# simulate data
set.seed(201)  # for reproducibility
n.sims <- 3  # number of simulations
sim.length <- 150  # length of simulation
sim.range <- 1:sim.length  # simulation range in years

sim.df <- data.frame(matrix(ncol = 3, nrow = 0))
sim.names <- c("year", "simulation", "flow")
colnames(sim.df) <- sim.names

for (i in 1:n.sims) {
  year <- sim.range
  simulation <- i
  flow <- simulate(fit.hmm, nsim = sim.length)$x 
  
  temp.df <- data.frame(year, simulation, flow)
  colnames(temp.df) <- sim.names
  
  sim.df <- rbind(sim.df, temp.df)
}

# plot simulated data
ggplot(data = sim.df) +
  geom_line(aes(x=year, y=flow, color=factor(simulation))) +
  ylab("Annual Flow (10^8 m^3)") +
  xlab("Simulation Year") +
  scale_color_discrete(name = "Simulation")






