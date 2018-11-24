%
% [spectral_features] = compute_spectral_features(eeg_data)
%
% input:
%   eeg_data            - A 1xN cell array containing the extracted EEG segments
%                       corresponding to imagined speech.
% output:
%   spectral_features   - A 1xN cell array containing the FFT features computed
%                       for each segment of EEG data.
%
function [spectral_features] = compute_spectral_features(eeg_data)
    
    fs = 1000; % sampling rate
    
    spectral_features = cell(size(eeg_data));
    for i=1:length(eeg_data)
        mat = eeg_data{i};
        spec_mat = zeros(size(mat, 1), 49);
        
        N = size(mat, 2); % number of samples.
        f = 0:fs/N:fs/2; % frequencies vector.
        rel_f = f >= 1 & f <= 50; % The band-passed filtered frequencies
        
        for j=1:size(mat, 1)
            x = mat(j, :);
            xdft = rfft(x);
            xdft = xdft(rel_f);
            bins = round(linspace(1, length(xdft), 50));
            for k=2:length(bins)
                inds = bins(k-1):bins(k);
                spec_mat(j, k-1) = sum(abs(xdft(inds)));
            end
        end
        spectral_features{i} = spec_mat;
    end
end