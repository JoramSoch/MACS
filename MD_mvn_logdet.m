function LD = MD_mvn_logdet(V,CF)
% _
% Log-Determinant for Multivariate Normal Covariance Matrix
% FORMAT LD = MD_mvn_logdet(V,CF)
% 
%     V  - an n x n multivariate normal covariance matrix
%     CF - a logical indicating use of Cholesky factorization
% 
%     LD - the log-determinant of that covariance matrix
% 
% FORMAT LD = MD_mvn_logdet(V) computes the logarithm of the determinant
% of the matrix V where V is an n x n square matrix. If V is singular,
% such that det(V) = 0, then it returns -Inf.
% 
% FORMAT LD = MD_mvn_logdet(V,true) computes the log-determinant using
% Cholesky factorization. This is applicable if V is positive definite,
% e.g. a covariance matrix.
% 
% Calling MD_mvn_logdet(V) is equivalent to calling log(det(V)). However,
% this function avoids the overflow or underflow problems that are likely
% to occur when applying det to large matrices.
% 
% The key idea of this algorithm is that the determinant of a triangular
% matrix equals the product of its diagonal elements. Therefore, the log-
% determinant of such a matrix is equal to the sum of their logarithms. By
% keeping all computations in log-space, overflow and underflow problems
% caused by products of many numbers can be effectively circumvented.
% 
% This implementation is based on LU factorization. When the matrix of
% interest is positive definite, such as a covariance matrix, Cholesky
% factorization can be selected which is typically more efficient.
%
% Log-Determinants of matrices occur in the context of multivariate
% probability distributions. The log-pdf, entropy and KL-divergence of the
% multivariate normal distribution include terms that have the form of log-
% determinants. Especially in high-dimensional space, computing the
% log-determinant directly can be very useful.
% 
% This function is based on code by Dahua Lin <dhlin@mit.edu> (2008) [1].
% 
% References:
% [1] http://www.mathworks.com/matlabcentral/fileexchange/22026-safe-
%     computation-of-logarithm-determinat-of-large-matrix/content/logdet.m
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 12/02/2015, 01:30 (V0.3/V10)
%  Last edit: 17/08/2015, 21:55 (V0.3/V11)


% Set Cholesky factorization
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(CF), CF = false; end;

% Compute log-determinant
%-------------------------------------------------------------------------%
if CF
    v = 2 * sum(log(diag(chol(V))));
else
    [L, U, P] = lu(V);
    d = diag(U);
    c = det(P) * prod(sign(d));
    v = log(c) + sum(log(abs(d)));
end;
LD = v;