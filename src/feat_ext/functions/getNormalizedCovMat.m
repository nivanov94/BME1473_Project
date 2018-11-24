function C = getNormalizedCovMat(E)

% Function returns the normalized covariance matrix
% for a matrix E where each row represents a different
% measurement and each column is a sample in time

C = E*E'/trace(E*E');

