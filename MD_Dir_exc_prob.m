function exc_p = MD_Dir_exc_prob(alpha)
% _
% Exceedance Probabilities for Dirichlet-distributed random variables
% FORMAT exc_p = MD_Dir_exc_prob(alpha)
% 
%     alpha - a 1 x K vector with parameters of the distribution
% 
%     exc_p - a 1 x K vector with Dirichlet exceedance probabilities
% 
% FORMAT exc_p = MD_Dir_exc_prob(alpha) returns exceedance probabilities,
% given that that the frequencies r follow a Dirichlet distribution,
% i.e. r ~ Dir(alpha) with sum(r) = 1. The exceedance probability for X,
% given that (X,Y) follows a multivariate distribution, is defined as the
% probability that X is larger than any element of Y [1].
% 
% Exceedance probabilities occur in Bayesian model selection [1,2] where
% one arrives at a posterior distribution over model frequencies that is
% Dirichlet. In this situation, it is highly useful to know the probability
% that a particular model X is more frequent than any of the other models Y.
% 
% If K = 2, the multivariate Dirichlet distribution becomes a univariate
% Beta distribution. In this case, we can use the Beta CDF to calculate
% exceedance probabilities.
% 
% If K > 2, no closed-form expression for exceedance probabilities exists.
% One option would be to use a sampling approach. Here, we use a numerical
% integral over several Gamma CDFs instead [3].
% 
% References:
% [1] Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ (2009):
%     "Bayesian model selection for group studies".
%     NeuroImage, vol. 46, pp. 1004-1017.
% [2] Penny WD, Stephan KE, Daunizeau J, Rosa MJ, Friston KJ, Schofield TM,
%     Leff AP (2010): "Comparing Families of Dynamic Causal Models".
%     PLoS ONE, vol. 6, iss. 3, e1000709.
% [3] Soch J & Allefeld C (2015): "Non-Critical Comments on
%     Bayesian Model Selection". Internal Report, June 2015.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 21/10/2014, 18:20 (V0.2/V6)
%  Last edit: 15/09/2016, 09:20 (V0.9a/V13a)


% Get dimensionality
%-------------------------------------------------------------------------%
K = numel(alpha);

% Analytical computation, if bivariate Dirichlet
%-------------------------------------------------------------------------%
if K == 2
    % using the Beta CDF
    exc_p(1) = 1 - betainc(1/2,alpha(1),alpha(2));
    exc_p(2) = 1 - exc_p(1);
end;

% Numerical integration, if multivariate Dirichlet
%-------------------------------------------------------------------------%
if K > 2
    % using Gamma CDFs
    exc_p = zeros(1,K);
    for j = 1:K
        f = @(x) integrand(x,alpha(j),alpha([1:K]~=j));
        % exc_p(j) = integral(f,0,Inf);
        % exc_p(j) = integral(f,0,Inf,'WayPoints',[alpha(j)]);
        exc_p(j) = integral(f,0,alpha(j)) + integral(f,alpha(j),Inf);
    end;
end;

% Integrand function for numerical integration
%-------------------------------------------------------------------------%
function p = integrand(x,aj,ak)

% product of Gamma CDFs
p = ones(size(x));
for k = 1:numel(ak)
    p = p .* gammainc(x,ak(k));
end;
% times a Gamma PDF
p = p .* exp((aj-1).*log(x) - x - gammaln(aj));
% p = p .* x.^(aj-1) .* exp(-x) ./ gamma(aj);