## Gillespie algorithm with time-dependent intensity function
## DRAFT - This is a work in progress.
##
## The implementation of the Gillespie in R here is 
## largely borrowed from Ben Bolker's excellent code, seen here:
## http://kinglab.eeb.lsa.umich.edu/EEID/eeid/gillespie.Rnw
##
## The key difference between this version and gillespie1.R (in this folder)
## is that the intensity function `intenfun2` which outputs the integral of
## the time-dependent intensity functions. In this example, I specify particular
## birth rates on particular days, corresponding to our observed Daphnia data.
## Naturally, this could be replaced with a time-dependent function or other
## array of values.
##
## Note that this may not be the most efficient way of running the algorithm
## in R and I welcome suggestions for speeding it up. Note also that this
## algorithm will run many times faster if implemented in Mathematica or
## a lower-level language such as C

intenfun2 <- function(X, param, t1, t2){
    k <- param[1]
    d <- param[2]                       # note the lack of a general birth term
    births <- function(t){              # enter our observed rates
        if(t < 16) 0
        if(t >= 16 || t < 17) {
            0.093 *  X * (1 - X / k)
        }
        if(t >= 17 || t < 18) {
            0.141 *  X * (1 - X / k)
        }
        if(t >= 18 || t < 19) {
            0.093 *  X * (1 - X / k)
        }
        if(t >= 19 || t < 20) {
            0.125 *  X * (1 - X / k)
        }
        if(t >= 20 || t < 21) {
            0.25 *  X * (1 - X / k)
        }
        if(t >= 21 || t < 22) {
            0.063 *  X * (1 - X / k)
        }
        if(t >= 22 || t < 23) {
            0.063 *  X * (1 - X / k)
        }
        if(t >= 23 || t < 24) {
            0.063 *  X * (1 - X / k)
        }
        if(t >= 24 || t < 25) {
            0.031 *  X * (1 - X / k)
        }
        if(t >= 25 || t < 26) 0
        if(t >= 26 || t < 27) {
            0.016 *  X * (1 - X / k)
        }
        if(t >= 27 || t < 29) 0
        if(t >= 29 || t < 30) {
            0.016 * X * (1 - X / k)
        }
    }
    deaths <- function(t){
        d * X
    }
    c(integrate(Vectorize(births), t1, t2)$value,
      integrate(Vectorize(deaths), t1, t2)$value)
}
   
gillespie2 <- function(init, times, intenfun2, pproc, param){
    tottime <- times[1]                 # total time of simulation (so far)
    tinc <- length(times)               # for indexing time increments
    N <- init                           # starting state of system
    results <- matrix(nrow = tinc, ncol = length(init))
    for(i in 1:(tinc - 1)){
        results[i, ] <- N               # record state at times[i]
        while(tottime < times[i + 1]){
            inten <- intenfun2(N, param, tottime, tottime + 1)
                                        # calculate current intensities
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

# Setting parameters described above
init <- c(25)
times <- seq(0, 30, by = 1)
pproc <- matrix(c(1, -1), ncol = 1, nrow = 2)
param <- c(200, 0.01)

# Run one simulation
gillespie2(init, times, intenfun2, pproc, param)
