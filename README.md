# Quorum-sensing simulation of the Sin system of Shinorizobium meliloti

This project contains three algorythms to simulate on Matlab the quorum sensing Sin system of Sinorhizobium meliloti .

## Model

blah blah

## Reaction Rates

The reaction rates are in the parameters.txt file. They come from the simulation of X species [Buceta et al.](http://example.com/ "Title") or from the data of [Bettenworth et al.](). If new literature comes out about the actual rates of the S. meliloti species, it's advised to update them.

## Deterministic simulations

The deterministic simulation is based on a set of ODES and the [ODEset solver of Matlab](insert link to that documentation).

myode_full is blah blah

myode_expr has the internal variable DNAexpr set to 0 to simulate a mutant.

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

As seen in example X, we used the [peakfinder.m]() algorythm to count the maximums.

## credit


