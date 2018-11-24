% [wav_features, feature_labels] = get_acoustic_features(dirname)
%
% input:
%   dirname         - A string containing the name of the directory containing 
%                   the wav files.
%   simple          - A boolean value that is true if we are just using simple
%                   descriptive feature and is false otherwise.
% output:
%   wav_features    - A 1xN cell array, where N is the number of wav files and
%                   element i contains the computed features for wav file i.wav
%   feature_labels  - A 1xD cell array, where D is the number of features and
%                   element i is the label for feature i.
function [wav_features, feature_labels] = get_acoustic_features(dirname, simple)
    
    % Getting the parameters.
    params = load_params('audio');
    
    % Selecting the function used to extract the features.
    if simple
        feat_func = @extract_simple_features;
    else
        feat_func = @extract_features;
    end

    % Getting the files.
    wav_fns = dir([dirname '/*.wav']);
    num_files = size(wav_fns, 1);
    i = 0;
    wav_features = cell(1, num_files);

    while i < num_files
        w = audioread([dirname '/' num2str(i) '.wav']);
        % Extracting the features.
        features = feat_func(w', params);
        wav_features{i+1} = features;
            
        i = i + 1;
    end

    if simple
        feature_labels = extract_face_labels();
    else
        feature_labels = extract_labels(params, 19);
    end
end