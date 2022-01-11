# On the efficiency of blind and non-blind estimation for coupled LL1 tensor models using the randomly-constrained Cramér-Rao bound

Copyright (c) 2022 Clémence Prévost, Eric Chaumette, Konstantin Usevich, Pierre Comon, David Brie <br>
Contact: ```clemence.prevost@univ-lille.fr```

This MATLAB software reproduces the results from the following article:
```
@unpublished{prevost:hal-03504402,
  TITLE = {{On the efficiency of blind and non-blind estimation for coupled LL1 tensor models using the randomly-constrained Cramér-Rao bound}},
  AUTHOR = {Pr{\'e}vost, Cl{\'e}mence and Usevich, Konstantin and Chaumette, Eric and Brie, David and Comon, Pierre},
  URL = {https://hal.archives-ouvertes.fr/hal-03504402},
  NOTE = {working paper or preprint},
  YEAR = {2021},
  MONTH = Dec,
  KEYWORDS = {Cram{\'e}r-Rao bounds ; random equality constraints ; tensor models ; low-rank approximations},
  PDF = {https://hal.archives-ouvertes.fr/hal-03504402/file/rccrb_tsp_v2.pdf},
  HAL_ID = {hal-03504402},
  HAL_VERSION = {v1},
}
```
<br><br>
[Link to the project](https://github.com/cprevost4/RCCRB_Software)

## Content

 - demo.m: demo file with minimal requirements 
 - /demos : contains demo files that produce tables and figures (including ```main.m```
 - /data : contains data for synthetic examples (Section VI.D)
 - /figures : where the figures are saved
 - /src : contains helpful files to run the demos

## Minimal requirements

 In order to run the demo file ```demo.m```, you will need to:
 - Download and install Tensorlab 3.0: https://www.tensorlab.net

 ## How it works
 
 ### Generate coupled tensor model
 
Every code starts by generating a coupled tensor model admitting LL1 decomposition. The entries of the LL1 factors are i.i.d. Gaussian variables with zero mean and unit variance.
The degradation in the first and second modes are blurring and downsampling operators. The degradation in the third mode is computed from the ```SRF_S2``` file, that corresponds to the spectral response function of the Sentinel-2A imaging sensor. 
 
 ### Compute the bounds
 
 The next step is to compute the bounds. In this software, we compute the standard Constrained Cramér-Rao bound (CCRB) and random CCRB (RCCRB) accounting for random equality constraints.
 
 ### Run estimators
 
 We evaluate the performance of the estimators using the above bounds. The MSE on the reconstructed tensor is computed by averaging the squared errors between the reference and low-rank estimate over a given number of trials.


## Reproduce figures and tables from the paper

To do so, you need to run the ```demos/main.m``` file.
Here, a menu is available and allows you to choose which figure or table you want to generate.
Each number in the table below corresponds to a set of figures, and/or tables.

| Number | Content                                        |
|--------|------------------------------------------------|
| 1      | produces Figure 1                              |
| 2      | produces Figure 2 and 3                        |
| 3      | produces Figure 4                              |
| 4      | produces Figure 5                              |
