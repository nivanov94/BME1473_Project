function W = getCSPProjMat(E,labels,normalized_passband)

% Function returns a matrix W that maximizes the 
% differentiation between two sets of whitened samples
% E - 3D struct with each matrix containing EEG recording of distinct sample
% labels - Vector of labels representing the trial captured by the EEG data

% begin by filter the each epoch
[b,a] = butter(10,normalized_passband);
Efilt = filter(b,a,E,[],2);


[C1, C2] = getAverageClassCovMats(Efilt,labels);

% sum the average covariance matrices to get composite cov mat
Cc = C1 + C2;

% decompose the composite covariance matrix to get the whitening matrix
[V, eigen] = eig(Cc);
P = sqrt(eigen)\V';

% calculate the rotation matrix that maximizes the differentiation between
% the two whitened classes
C1h = P*C1*P';
[U, ~] = eig(C1h);

W = U'*P; 
