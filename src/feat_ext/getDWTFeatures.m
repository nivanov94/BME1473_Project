function [feat_vector,l] = getDWTFeatures(x, wname, levels, features, Fs)

% Function performs DWT transformation decomposition on a vector of samples
% x and returns the decomposition coefficients
%   x         - The vector of samples
%   wname     - The type of wavelet to use
%   levels    - The number of decompositon levels to use
%   features  - Vector of features to return. Encoded as:
%                   1: Raw DWT coefficients
%                   2: RMS of DWT coefficients at each decomposition level
%                   3: SD of DWT coefficients at each decomposition level
%                   4: Energy of the DWT coefficients at each decompositon level
%                   5: Maximum value of DWT coefficients at each decomposition level
%                   6: Minimum value of DWT coefficients at each decomposition level
%  Fs         - The sampling frequency used to obtain the signal (only needed if you
%               want the energy of the coefficients.

if ((sum(features==4) > 0) && nargin < 5) || (nargin < 4)
  help getDWTFeatures
else

  [all_coeffs,l] = wavedec(x,levels,wname);
  
  % remove any repeated features requests
  features = unique(features,'stable');
  
  % separate the decomposition levels into different cell arrays
  dlevel_coeffs = cell(1,levels+1);
  dlevel_coeffs{1} = appcoef(all_coeffs,l,wname);
  for i = 2:(levels+1)
    dlevel_coeffs{i} = detcoef(all_coeffs,l,i-1);
  end
  
  % compute the length of the feature vector and the starting index of each feature set
  l = zeros(1,length(features));
  l(1) = 1;
  feat_vec_len = 0;
  for i = 1:length(features)
    if features(i) == 1
      feat_len = length(all_coeffs);
    else
      feat_len = levels + 1;
    end
  
    if i < length(features)
      l(i+1) = l(i) + feat_len;
    else
      feat_vec_len = l(i) + feat_len - 1;
    end
  end
  
  feat_vector = zeros(1,feat_vec_len);
  for i = 1:length(features)
    if (i==length(features)) end_index = feat_vec_len; else; end_index = l(i+1)-1; end;
    switch features(i)
      case 1
        % add the raw DWT coefficients to the feature vector
        feat_vector(l(i):end_index) = all_coeffs';
      case 2
        % add the RMS of each DWT decomposition level to the feature vector
        cellrms = @(x) rms(cell2mat(x));
        feat_vector(l(i):end_index) = arrayfun(cellrms,dlevel_coeffs);
      case 3
        % add the SD of each DWT decomposition to the feature vector
        cellstd = @(x) std(cell2mat(x));
        feat_vector(l(i):end_index) = arrayfun(cellstd,dlevel_coeffs);
      case 4
        % add the total energy of the DWT coeffeicients to the feature vector
        dt = 1/Fs;
        cellenergy = @(x) trapz((1:length(cell2mat(x)))*dt,abs(cell2mat(x)).^2);
        feat_vector(l(i):end_index) = arrayfun(cellenergy,dlevel_coeffs);
      case 5
        % add the minimum value of each DWT coefficient to the feature vector
        cellmin = @(x) min(cell2mat(x));
        feat_vector(l(i):end_index) = arrayfun(cellmin,dlevel_coeffs);
      case 6
        % add the maximum value of each DWT coefficient to the feature vector
        cellmax = @(x) max(cell2mat(x));
        feat_vector(l(i):end_index) = arrayfun(cellmax,dlevel_coeffs);
    end
  end
end
