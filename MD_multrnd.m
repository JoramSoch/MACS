function r = MD_multrnd(p,N)
% _
% Multinomially distributed random numbers
% FORMAT r = MD_multrnd(p,N)
% 
%     p - a  1 x K vector with multinomial probabilities
%     N - an integer specifying the number of samples
% 
%     r - an N x 1 vector with multinomial random numbers
% 
% FORMAT r = MD_multrnd(p,N) returns N samples from the multinomial
% distribution with probablities p. Each sample is a number between and
% including 1 and K. Random numbers are generated using a cumulative
% probability approach based on the standard uniform distribution [1].
% 
% References:
% [1] Gelman A, Carlin JB, Stern HS, Dunson DB, Vehtari A, Rubin DB (2013):
%     "Bayesian Data Analysis". Chapman & Hall, 3rd edition, pp. 583-584.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 20/11/2014, 17:40 (V0.2/V8)
%  Last edit: 09/03/2018, 11:25 (V1.2/V18)


% Get dimensionality
%-------------------------------------------------------------------------%
K = numel(p);

% Calculate cumulative probabilities
%-------------------------------------------------------------------------%
c = zeros(1,K+1);
for j = 1:K
    c(j+1) = sum(p(1:j));
end;

% Generate multinomial random numbers
%-------------------------------------------------------------------------%
r = zeros(N,1);
for i = 1:N
    r(i) = sum(rand > c);
end;