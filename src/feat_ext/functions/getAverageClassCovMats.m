function [C1,C2] = getAverageClassCovMats(E,labels)

% Function returns the average normalized covariance
% matrices for EEG recordings of two different classes.
% E - 3D structure with each matrix representing EEG
%     recordings from a different trial
% labels - Vector containing the class labels for each
%          trial (labels should be 0 or 1)


C1_samples = 0;
C2_samples = 0;

C1 = zeros(size(E{1},1));
C2 = zeros(size(E{1},1));

% iterate through all samples
for i = 1:length(E)
  if labels(i) == 0
    C1 = C1 + getNormalizedCovMat(E{i});
    C1_samples = C1_samples + 1;
  else
    C2 = C2 + getNormalizedCovMat(E{i});
    C2_samples = C2_samples + 1;
  end
end

C1 = C1/C1_samples;
C2 = C2/C2_samples;
