%           On the efficiency of blind and non-blind estimation
%                for coupled LL1 tensor models using the
%                 randomly-constrained Cramér-Rao bound                   %
%-------------------------------------------------------------------------%

% Copyright (c) 2022 Clemence Prevost, Konstantin Usevich, Eric Chaumette,
% Pierre Comon, David Brie
% https://github.com/cprevost4/RCCRB_Software
% Contact: clemence.prevost@univ-lille.fr

% This software reproduces the results from the preprint called:
% (1) "On the efficiency of blind and non-blind estimation
% for coupled LL1 tensor models using the randomly-constrained
% Cramér-Rao bound- C.Prévost, K.Usevich, E. Chaumette, 
% D.Brie, P.Comon.

% In order to run the demo, you will need to add to your MATLAB path:
% - Tensorlab 3.0: https://www.tensorlab.net

%-------------------------------------------------------------------------%
%                              CONTENT
% - /data : contains data for degradation matrices and random vectors
% - /demos : contains demo files that produce tables and figures
% - /figures : where the tables and figures are saved
% - /src : contains helpful files to run the demos
%-------------------------------------------------------------------------%
%                                MENU
% You can launch a specified demo by typing its number. The resulting tables
% and figures produced will be stored in the figures folder.
%
% 1:  produces Figure 1 
% 2:  produces Fig. 2 and 3
% 3:  produces Fig. 4
% 4:  produces Fig. 5

%-------------------------------------------------------------------------%

list_demos = ["deterministic" "conditional_rccrb" "mismatch" "recap"];

prompt = "Which file do you want to run ?";
num = input(prompt);
eval(list_demos(num));