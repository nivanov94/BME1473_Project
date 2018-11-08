function feat_vector = getMomentFeatures(x,features)
% Function calculates p-order moments for a given signal x.
%    x        - The vector of samples for the signal of interest
%    features - Vector containing the moments to calculate
%

% remove any repeated feature requests
features = unique(features,'stable');

% calcualate the length of the feature vector
feat_vector = zeros(1,length(features));

for i = 1:length(features)
  feat_vector(i) = sum(x.^features(i));
end

feat_vector = feat_vector/length(x);
