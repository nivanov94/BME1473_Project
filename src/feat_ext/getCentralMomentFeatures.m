function feat_vector = getCentralMomentFeatures(x,features)
% Function calculates p-order central moments for a given signal x.
%    x        - The vector of samples for the signal of interest
%    features - Vector containing the central moments to calculate
%

% remove any repeated feature requests
features = unique(features,'stable');

% calcualate the length of the feature vector
feat_vector = zeros(1,length(features));

for i = 1:length(features)
  feat_vector(i) = moment(x,features(i));
end
