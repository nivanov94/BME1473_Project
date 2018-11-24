%
% [spk_mats] = trim_speaking_mats(speaking_mats, wav_folder)
%
% input:
%   speaking_mats   - A 1xN cell array, where N is the number of matrices
%   wav_folder      - The folder containg the wav files, and there should be a
%                   total of N wav files.
% output:
%   spk_mats        - The same as speaking_mats, except with the matrices
%                   trimmed to contain only the "speaking" segments of the EEG
%                   data.
function [spk_mats] = trim_speaking_mats(speaking_mats, wav_folder)
    
    wav_fns = dir([wav_folder '/*.wav']);
    eeg_fs = 1000; % The sampling rate for EEG data, which is 1 kHz
    assert(size(wav_fns, 1) == size(speaking_mats, 2));
    spk_mats = {};
    num_files = size(wav_fns, 1);
    j = 0;
    
    while j < num_files
        mat = speaking_mats{j+1};
        [w, Fs] = audioread([wav_folder '/' num2str(j) '.wav']);
        s = length(w) / Fs;
        num_samples = floor(s * eeg_fs);
        e_ind = min(num_samples, size(mat, 2));
        mat = mat(:, 1:e_ind);
        spk_mats = [spk_mats mat];
        
        j = j + 1;
    end 
end