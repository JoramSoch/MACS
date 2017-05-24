function r_mode = MD_Dir_mode(alpha)
% _
% Likeliest Frequencies for Dirichlet-distributed random variables
% FORMAT r_mode = MD_Dir_mode(alpha)
% 
%     alpha  - a 1 x K vector with parameters of the distribution
% 
%     r_mean - a 1 x K vector with mode or likeliest frequencies
% 
% FORMAT r_mode = MD_Dir_mode(alpha) returns the modes, i.e. likeliest
% frequencies [r], given that the frequencies r follow a Dirichlet
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
% First edit: 04/12/2014, 14:30 (V0.2/V8)
%  Last edit: 15/09/2016, 09:15 (V0.9a/V13a)


% Calculate likeliest frequencies
%-------------------------------------------------------------------------%
r_mode = (alpha-1)./(sum(alpha)-numel(alpha));