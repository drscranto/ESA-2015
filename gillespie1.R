## Gillespie algorithm with time-independent intensity function
## 
## The implementation of the Gillespie in R here is 
## largely borrowed from Ben Bolker's excellent code, seen here:
## http://kinglab.eeb.lsa.umich.edu/EEID/eeid/gillespie.Rnw
##
## I include all relevant functions in this file to ease
## understanding, but of course you may prefer to have
## individual functions as separate files.
## 
## This code produces a stochastic realization of a simple
## birth-death model of population growth.
##
## `intenfun` contains the intensity functions for birth and death.
## In this case, birth is density-dependent and death is not.
##
## The larger function `gillespie` takes as input the following:
## init     = An array containing all initial conditions
##            In this case, it contains only the starting population size X
## times    = Times over which you wish the function to report population size
## intenfun = Calculates intensities/probabilities of different events
## pproc    = An array containing the state changes caused by the point process
##            In this case, a 2 x 1 array (2 events, 1 compartment/stage)
## param    = An array containing all intrinsic demographic rates required for
##            calculating intensities in `intenfun`
##
## Feel free to send questions to geoffrey.legault@colorado.edu

intenfun <- function(X, param){
    b <- param[1]
    k <- param[2]
    d <- param[3]
    births <- b * X * (1 - X / k)
    deaths <- d * X
    c(births, deaths)
}

gillespie <- function(init, times, intenfun, pproc, param){
    tottime <- times[1]                 # total time of simulation (so far)
    tinc <- length(times)               # for indexing time increments
    N <- init                           # starting state of system
    results <- matrix(nrow = tinc, ncol = length(init))
    for(i in 1:(tinc - 1)){
        results[i, ] <- N               # record state at times[i]
        while(tottime < times[i + 1]){
            inten <- intenfun(N, param) # calculate current intensities
            if(all(inten == 0) ||
               min(inten) < 0) break    # stop if extinction/errors
            allinten <- sum(inten)
            deltat <- rexp(1, allinten) # generate next time step
            which.pproc <- sample(1:nrow(pproc),
                                  size = 1,
                                  prob = inten)
                                        # select point process
            tottime <- tottime + deltat
            N <- N + pproc[which.pproc, ]
        }
    }
    results[tinc, ] <- N                # record last state
    cbind(times, results)               # put things together
}

# Set parameters described above
init <- c(25)
times <- seq(0, 100, by = 2)
pproc <- matrix(c(-1, 1, 1, 0, 0, -1), ncol = 2, nrow = 3)
param <- c(.068, 200, 0.01)

# Run one simulation
gillespie(init, times, intenfun, pproc, param)
