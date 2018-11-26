% [eeg_features, feature_labels] = get_eeg_features(epoch_data, W, art_inds)
%
% input:
%   epoch_data      - A struct containing the epoched data for the 4 different
%                   states: thinking, clearing, speaking, and stimuli. The data
%                   may contain the raw EEG signals or the ICs.
%   W               - The 64x64 mixing matrix. W can be the identity matrix if
%                   we only want to use the original EEG signals.
%   num_prompts     - A numeric class label
% output:
%   eeg_features    - A struct with a similar format to epoch_data containing
%                   the computed features for the data in epoch_data.
%   feature_labels  - A 1xD cell array, where D is the number of features and
%                   element i is the label for feature i.
function [eeg_features, feature_labels] = get_eeg_features(epoch_data, W)
    
    % loading the parameters.
    params = load_params('eeg');
    
    % Selecting the function used to extract the features.
    feat_func = @extract_features;
    
    % Sanity check.
    num_epochs = length(epoch_data.thinking_mats);
    assert(num_epochs == length(epoch_data.thinking_mats));
    assert(num_epochs == length(epoch_data.clearing_mats));
    
    thinking_feats = cell(1, num_epochs);
    clearing_feats = cell(1, num_epochs);
    
    for i=1:num_epochs
        disp(['Epoch ' num2str(i) ' of ' num2str(num_epochs)]);
        
        thinking_mat = W * epoch_data.thinking_mats{i};
        clearing_mat = W * epoch_data.clearing_mats{i};
        
        num_chans = size(thinking_mat, 1);
        % Initiliazing the feature matrices.
        empty = true;
        thinking_feat = [];
        clearing_feat = [];
        
        for j=1:num_chans
            % A small speedup
            if empty           
                t_f = feat_func(thinking_mat(j, :), params);
                num_feats = length(t_f);
                thinking_feat = zeros(num_chans, num_feats);
                clearing_feat = zeros(num_chans, num_feats);
                
                thinking_feat(j, :) = t_f;
                clearing_feat(j, :) = feat_func(clearing_mat(j, :), params);

                empty = false;
            else
                thinking_feat(j, :) = feat_func(thinking_mat(j, :), params);
                clearing_feat(j, :) = feat_func(clearing_mat(j, :), params);
            end
        end
        
        thinking_feats{i} = thinking_feat;
        clearing_feats{i} = clearing_feat;
    end
    
    eeg_features = struct();
    eeg_features.thinking_feats = thinking_feats;
    eeg_features.clearing_feats = clearing_feats;
    
    feature_labels = extract_labels(params, 19);
end
