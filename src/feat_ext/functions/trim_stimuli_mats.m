%
% [stimuli_mats] = trim_stimuli_mats(speaking_mats)
%
% input:
%   speaking_mats   - A 1xN cell array, where N is number of EEG samples where
%                   the participant was presented the stimulus.
function [stimuli_mats] = trim_stimuli_mats(speaking_mats)
    
    Fs = 1000; % The sampling rate for the EEG data, which is 1 kHz
    % Sometimes, the number of samples is less than 2000, which is weird.
    s = 2; % The number of non-stimuli seconds.
    stimuli_mats = {};
    
    for i=1:size(speaking_mats, 2)
        mat = speaking_mats{i};
        num_samples = size(mat, 2) - s*Fs;
        mat = mat(:, 1:num_samples);
        stimuli_mats = [stimuli_mats mat];
    end
end