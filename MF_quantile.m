function q = MF_quantile(y,p)
% _
% Quantiles for univariate data
% FORMAT q = MF_quantile(y,p)
% 
%     y - a vector of sample data from a random variable
%     p - a k x 1 vector of percentages between 0 and 1
% 
%     q - a k x 1 vector of quantiles for these percentages
% 
% FORMAT q = MF_quantile(y,p) returns quantiles for data y at cut-offs p.
% For example, the q for p = 0.5 would be the median, a value chosen such
% that 50% of the data are smaller than this threshold.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 31/07/2014, 11:50 (V0.1/V4)
%  Last edit: 13/11/2014, 12:00 (V0.2/V7)


% Sort data
%-------------------------------------------------------------------------%
ys = sort(y);

% Get quantiles
%-------------------------------------------------------------------------%
q = NaN(size(p));
for j = 1:length(p)
    i = 1;
    n = numel(ys);
    while (i/n) < p(j)
        i = i + 1;
    end;
    q(j) = ys(i);
end;