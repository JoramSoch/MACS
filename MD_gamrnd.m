function R = MD_gamrnd(a,b)
% _
% Gamma-distributed random numbers
% FORMAT R = MD_gamrnd(a,b)
% 
%     a - a matrix of Gamma shape parameters ("alpha")
%     b - a matrix of Gamma rate  parameters ("beta")
% 
%     R - a matrix with Gamma-distributed random numbers
% 
% FORMAT R = MD_gamrnd(a,b) returns a matrix of random numbers for the
% Gamma distributions specified by a and b. The matrices a and b need to
% have the same size, unless one of them is scalar. The matrix R will have
% the same size as a or b. Random numbers are generated using a rejection
% method [1].
% 
% References:
% [1] Marsaglia G, Tsang WW (2000): "A Simple Method for Generating Gamma
%     Variables". ACM Transactions Mathematical Software, vol. 26, no. 3,
%     pp. 363-372. URL: http://portal.acm.org/citation.cfm?id=358414.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 27/11/2014, 18:20 (V0.2/V8)
%  Last edit: 04/12/2014, 17:20 (V0.2/V8)


% Return NaN for invalid parameter values
%-------------------------------------------------------------------------%
a(a < 0) = NaN;
b(b < 0) = NaN;

% Return random numbers for the gamma distribution Gamma(a,b)
%-------------------------------------------------------------------------%
R = sgamrnd(a)./b;

% Get random numbers for standard gamma distribution Gamma(a,1)
%-------------------------------------------------------------------------%
function r = sgamrnd(a)

r = zeros(size(a));
for i = 1:numel(a)
    % (1) Setup
    d = a(i)-1.0/3;
    c = 1.0/sqrt(9*d);
    while true
        % (2) Generate normal
        v = 0;
        while v <= 0
            x = randn;
            v = 1+c*x;
        end;
        v = v^3;
        x = x^2;
        % (3) Generate uniform
        u = rand;
        % (4)/(5) Random number
        if u < (1-0.0331*x^2) || log(u) < (0.5*x + d*(1-v+log(v)))
            break
        end;
    end;
    r(i) = d*v;
end;