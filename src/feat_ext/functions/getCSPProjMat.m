function W = getCSPProjMat(Cc,C1)

% decompose the composite covariance matrix to get the whitening matrix
[V, eigen] = eig(Cc);
P = sqrt(eigen)\V';

% calculate the rotation matrix that maximizes the differentiation between
% the two whitened classes
C1h = P*C1*P';
[U, ~] = eig(C1h);

W = U'*P; 
