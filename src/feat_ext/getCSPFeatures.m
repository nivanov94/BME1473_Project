function features = getCSPFeatures(E,labels,bands,m,Fs)

% Function returns matrix of common spatial pattern scores
% with each row containing the scores for a single trial.
% E - multidimensional array containing matrices of different EEG samples
% labels - The numerical class label for each sample
% bands - List of frequency band tuples containing frequency range for bandpass filters
% m - Number of scores to extract
% Fs - The sampling frequency

% output 2m features for each class in each sample in each frequency band
unique_labels = unique(labels);

features = zeros(size(E,1),size(bands,1),length(unique_labels),2*m);

% iterate over all frequency bands
for band = 1:size(bands,1)
  [b,a] = butter(10,bands(band,:)/(Fs/2));
  E_filt = filter(b,a,E,[],2); % filter the readings from each electrode from each sample
  % iterate over all classes
  for class = 1:length(unique_labels)
    % set class labels equal to the current class to 1, others to 0
    class_labels = labels==unique_labels(class);

    % get the CSP projection matrices for these classes
    Wcsp = getCSPProjMat(E_filt,class_labels);

    % reduce the projection matrix to the requrested number of rows
    Wcsp2m = zeros(2*m,size(Wcsp,2));
    Wcsp2m(1:m,:) = Wcsp(1:m,:);
    Wcsp2m(m+1:end,:) = Wcsp(end-m+1:end,:);

    % iterate over all samples
    for s = 1:length(labels)
      Z = Wcsp2m*squeeze(E_filt(s,:,:));

      % now calculate the variance of each row and generate the 2m scores
      v = var(Z');
      features(s,band,class,:) = log10(v/sum(v));
    end
  end
end

