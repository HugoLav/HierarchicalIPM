# HierarchicalIPM
Repository for the article "Hierarchical Integral Probability Metrics: A distance on random probability measures with low sample complexity", accepeted for presentation in ICML. 
 
The article can be found on [arxiv](https://arxiv.org/abs/2402.00423).  

Implementation in Julia of the Wasserstein over Wasserstein distance and of the Lipschitz HIPM dlip. Details of the algorithms can be found in the article. 

The code for the Wasserstein over Wasserstein distance can be found in  `distance_Wasserstein.jl` while the code for the distance dlip can be found in `new_distance.jl`. The data for the figures of the article can be regenerated by running the files `Figurek_*****.jl` for k=1,2,3.

The code requires the following Julia packages: DelimitedFiles, Distributions, ExactOptimalTransport, LinearAlgebra, Statistics, Tulip.