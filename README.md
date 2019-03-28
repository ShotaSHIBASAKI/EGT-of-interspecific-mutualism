# Shibasaki 2018

The codes here are used in the simulation of the manuscript "The evolutionary game theory 
of interspecific mutualism in the multi-species communities" by S. Shibasaki (2019) in J. Theor. Biol. 
The preprint is available at bioRxiv.

## Description

The codes are for simulating the evolutionary dynamics in the communities where three or 
four species coexist, they mutually interact with one anoter, and their population sizes 
are infinitely large. In addition, I provide a code where populaiton sizes are different 
in each species in the three species model.  The last code can comupte the favorability
of each species in the community where arbitary number of species coexist.
These codes are writen in Python 3.

Note that these programs have large computational costs because they solve the 
ordinary differential equations from various intial conditions. Please check the initial 
conditions and the reputations.
## MC_RQRK.py

This code simulate the evolutionary dynamics of given number of specues (M) with a given value of parameter k
from 1,000 samples of initial conditions with ten replicate.
In other words, this code run the Monte Carlo simulation to estimate the favorability
(probability that each species becomes selfish x=0). Some csv files whose name begins with "MC-" are the results that
are used in the paper.

##  3species_infinite.py

This code is for simulating the evolutionary dynamics in three species community 
where population sizes are infinitely large.  
To analyze the probability that the dynamics from each initial conditions convege to 
stable equilibria, run "3species_infinite_stable_equilibrium_analysis.py". 
In addition, "3species_infinite_analyze.py" enables us to classify the results according 
to the first fixed species.

## 4species_infinite.py

This code is for simulating the evolutionary dynamics in four species community 
where population sizes are infinitely large. 
Note that this code takes long time to finish the simulation.

## 3species_various_population_size.py

This code is for simulating the evolutionary dynamics in four species community 
where population sizes of species are various. By changing the value of input "model", 
we can change the values of populaiton sizes.

## CLT-approximation.py

This code represent whether the central limit theorem offers a good approximation or not when the number of species M is small.

## ideal-condition.py

This code computes the favorability of each species under the ideal conditions 
(you can find the explanaitons for the ideal conditions on the manuscript). 
A function Image can drow  a figure which shows whether a focal species is 
more likely to be generous or selfish.
With a function Analyze (M,k) one can calculate the favorability of each species
given the number of species M and the value of k. 
Analyze2 make the favorability of each species in the three or four species with deifferent values of k. 
This function was used to draw Fig. 4 in the manuscript.

## favorability.py
This code provides the favorabilioty of each species given the value of M.
The favorabilities are calculated by using two values of k (k=0.5 and k=1.5).
With this code, you can find the relationdhip between the evolutionary rates and evolutionary fate with 
an arbitrary value of M, assuming that the difference of evolutionary rates is quite large.
