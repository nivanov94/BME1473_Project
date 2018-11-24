function [feat_vector,l] = getARFeatures(x,order,features,bands,Fs)

% Function computes auto-regressive coefficients for the signal x using the
% Yule-Walker equation.
%   x        - The vector containing the signal of interest
%   order    - The order of the AR model to compute
%   features - Vector of features to return. Encoded as:
%                 1: Raw AR coefficients
%                 2: Energy of AR coefficients in different frequency bands
%   bands    - Matrix defining the frequency bands of interest. Column 1
%              is starting frequency (Hz), column 2 is ending frequency (Hz).
%   Fs       - The sampling frequency used to obtain x

if ((sum(features==2) > 0) && nargin < 5) || (nargin < 3)
  help getARFeatures
else

  coeffs = ar(x,order,'yw');
  coeffs = coeffs.A;

  % remove any repeated feature requests
  features = unique(features,'stable');

  % compute the length of the feature vector
  l = ones(1,length(features));
  feat_vec_len = any(features==1)*(order+1) + any(features==2)*size(bands,1);

  for i = 1:(length(features)-1)
    if features(i) == 1
      feat_len = order+1;
    else
      feat_len = size(bands,1);
    end

    if i < length(features)
      l(i+1) = l(i) + feat_len;
    end
  end

  % generate the feature vector
  feat_vector = zeros(1,feat_vec_len);
  for i = 1:length(features)
    if (i==length(features)); end_index = feat_vec_len; else; end_index = l(i+1)-1; end;
    switch features(i)
      case 1
	% add the AR coefficients to the feature vector
	feat_vector(l(i):end_index) = coeffs;
      case 2
        % compute the energy of the AR coefficients within different frequency bands
        df = 0.01; % not sure what to use for this value yet
        band_energies = zeros(1,size(bands,1));
        for j = 1:size(bands,1)
          f = bands(j,1):df:bands(j,2);
          h = freqz(1,coeffs,f,Fs);
          h = abs(h).^2;
          band_energies(j) = trapz(f,h);
        end
        feat_vector(l(i):end_index) = band_energies;
    end
  end
end

