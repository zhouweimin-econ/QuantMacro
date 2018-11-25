function [c_guess, inc]=policyguess(meshes,WW,par)
%policyguess returns autarky policy guesses (in the first period only):
% Consumption is compositite leisure and physical consumption (x_it) in the
% paper, therefore labor income is reduced by the fraction of leisure
% consumed.

inc.labor   = WW.*meshes.h;
inc.money   = par.RB.*meshes.m;
inc.profits = 0;%lump sum profits

%% Initial policy guesses:  Autarky policies as guess

% Consumption guess
c_guess = inc.labor + max(inc.money,0) + inc.profits;


end
