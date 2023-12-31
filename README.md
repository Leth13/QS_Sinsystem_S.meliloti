# Quorum-sensing simulation of the Sin system of _Sinorhizobium meliloti_

This project contains three algorythms to simulate on Matlab the quorum sensing Sin system in _Sinorhizobium meliloti_.

## Model

There are a total of 23 reactions characterizing the model.

The transcription and translation reactions were simplified in one reaction with one rate, beta.

## Reaction Rates

The reaction rates are in the [parameters.txt](parameters.txt) file. They were extracted from the data of [Bettenworth et al.](https://doi.org/10.1038/s41467-022-30307-6) or, when information was lacking, assumed identical to the simulation of the LuxR/LuxI quorum sensing system of _Vibrio fisheri_ in [Buceta et al.](https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-7-6). If new literature comes out about the actual rates of the _S. meliloti_ species, it's advised to update them.

## Deterministic simulations

The deterministic simulation is based on a set of ODES resolved by the [ODEset solver of Matlab](https://it.mathworks.com/help/matlab/ref/odeset.html).

It needs all 15 reactions rates and the starting number of all metabolite/genes, except ` DNAexpr `, which is an internal variable.

[myode_full.m](myode_full.m) is the wild-type model of the Sin QS system.

[myode_expr.m](myode_expr.m) simulates a mutant lacking the _expr_ gene by setting the internal variable ` DNAexpr = 0 `.

## Stochastic simulation

[SSA_QS.m](SSA_QS.m) contains a single-cell [gillespie style simulation](https://en.wikipedia.org/wiki/Gillespie_algorithm) of the Sin QS system.

It needs all 15 reactions rates and the starting number of all metabolite/genes. The internal variable V, extracted from [Xiongfeng et al.](https://journals.asm.org/doi/abs/10.1128/msphere.00567-18), accounts for the statistical effect of the volume of the cell on reactions of one order or more.

`
  function [SINIplot, SINRplot, EXPRstarplot, EXPRplot, timeplot, meantime] = SSA_QS(y, parameters, numberofcells, numberofreactions)
`

Mutants can be simulated by changing the values of the reactions rates or the genes.

## Confronting simulations

See the file [example_confronting_simulations.m](example_confronting_simulations.m).

## Frequency of SinI expression

We used the [peakfinder](https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder-x0-sel-thresh-extrema-includeendpoints-interpolate) function to calculate the frequency of pulses of expression of SinI.

See [example_stat_analysis_pulses.m](example_stat_analysis_pulses.m).

## Credit


