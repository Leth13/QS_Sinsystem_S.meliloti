# Quorum-sensing simulation of the Sin system of Shinorizobium meliloti

This project contains three algorythms to simulate on Matlab the quorum sensing Sin system of Sinorhizobium meliloti .

## Model

blah blah

## Reaction Rates

The reaction rates are in the parameters.txt file. They were extracted from the data of [Bettenworth et al.](https://doi.org/10.1038/s41467-022-30307-6) or, where information was lacking, assumed from the simulation of the LuxR/LuxI quorum sensing system of Vibrio fisheri in [Buceta et al.](https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-7-6). If new literature comes out about the actual rates of the S. meliloti species, it's advised to update them.

## Deterministic simulations

The deterministic simulation is based on a set of ODES and the [ODEset solver of Matlab](https://it.mathworks.com/help/matlab/ref/odeset.html).

myode_full is blah blah

myode_expr has the internal variable 'DNAexpr = 0' to simulate a mutant lacking the expr gene.

  example code

## Stochastic simulation

The SSA_QS.m contains the gillespie style simulation (insert link) for all tot reaction of the Sin QS system. 
It needs 
  blah blah 
and spouts 
  blah blah

## Confronting simulations

See examples

## How to count pulses

We used the [peakfinder](https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder-x0-sel-thresh-extrema-includeendpoints-interpolate) function to calculate the frequency of pulses of expression of SinI.

example code.

As seen in example X, we used the [peakfinder.m]() algorythm to count the maximums.

## credit


