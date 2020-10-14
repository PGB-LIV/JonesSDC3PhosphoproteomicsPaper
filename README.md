# JonesSDC3PhosphoproteomicsPaper

## Requirements
 - [Julia language v0.5.2](https://julialang.org/downloads/oldreleases/#v052_may_6_2017)
 - [CmdStan](https://mc-stan.org/users/interfaces/cmdstan)

## Julia Package Requirements
 - FileIO
 - JLD
 - DataFrames
 - Stan
 - Mamba
 - Distributions
 
## Usage
Running:
```
> julia importing.jl
```
will split the Progenesis.csv data file into subsets for each of the four initial comparisons

Then:
```
> julia progenesis-processing.jl
```
will perform linear regression for missing value imputation and differential expression for each of the four initial comparisons.

Finally:
```
> julia fcdifference.jl
```
uses the results form the previous step to calculate the results for comparison #5 (difference in fold-changes).


 
  