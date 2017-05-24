function r_mean = MD_Dir_mean(alpha)
% _
% Expected Frequencies for Dirichlet-distributed random variables
% FORMAT r_mean = MD_Dir_mean(alpha)
% 
%     alpha  - a 1 x K vector with parameters of the distribution
% 
%     r_mean - a 1 x K vector with mean or expected frequencies
% 
% FORMAT r_mean = MD_Dir_mean(alpha) returns the means, i.e. expected
% frequencies <r>, given that the frequencies r follow a Dirichlet
% distribution, i.e. r ~ Dir(alpha) with sum(r) = 1. Formulas are given
% in standard textbooks [1].
% 
% References:
% [1] Gelman A, Carlin JB, Stern HS, Dunson DB, Vehtari A, Rubin DB (2013):
%     "Bayesian Data Analysis". Chapman & Hall, 3rd edition, pp. 578-579. 
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 21/10/2014, 18:10 (V0.2/V6)
%  Last edit: 15/09/2016, 09:15 (V0.9a/V13a)


% Calculate expected frequencies
%-------------------------------------------------------------------------%
r_mean = alpha./sum(alpha);