# ESA 2015

Included in this repository are two implementations of the Gillespie algorithm, meant to illustrate the differences between time-independent demography (the standard assumption of the Gillespie algorithm) and time-dependent demography.

# Code

## gillespie1.R

This is the standard Gillespie algorithm, which assumes that demographic rates/intensities do not change over time. If this assumption is violated, it does not properly sample the waiting time to the next event.

## gillespie2.R

This is a modified Gillespie algorithm, which allows demographic rates/intensities to change with time and which correctly samples the waiting time to the next event.

# Errors?

This is draft code so I apologize for any errors it may contain. Please let me know about any such cases and I will try to fix them as soon as possible.

Email: geoffrey.legault@colorado.edu
