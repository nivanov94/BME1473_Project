
% Add the feature extraction functions to the path 
addpath(genpath('D:\School\GRAD\BME1473\Project\code\src\feat_ext\functions'));
% Add the EEGLAB functions to the path
addpath(genpath('C:\Users\Nick\Documents\MATLAB\eeglab\eeglab14_1_2b\functions'));

% The folders for which we want to compute the features.
% folders = {'P02' 'MM05' 'MM08' 'MM09' 'MM10' 'MM11' 'MM12' 'MM14' 'MM15' 'MM16' 'MM18' 'MM19' 'MM20' 'MM21'};
folders = {'P02'};
data_path = 'D:\School\GRAD\BME1473\Project\data';
misc_folder = [data_path '\misc'];
ICA = true;

for j=1:length(folders)
    disp(['Computing features for folder ' folders{j}]);
    folder = [data_path '\' folders{j}];
    if ~ICA
        D = dir([folder '\' folders{j} '.set']);
    else
        D = dir([folder '\' folders{j} '-ICA.set']);
    end
    % there should only be 1 set file.
    set_path = [folder '\'];
    %set_fn = [folder '\' D(1).name];
    
    kinect_folder = [folder '\kinect_data'];
    labels_fn = [kinect_folder '\labels.txt'];

    % split the file, create the epochs, and save.
    EEG = pop_loadset(D(1).name,set_path);
    data = EEG.data;
    
    % Applying the ICA decomposition.
    if ICA
        W = EEG.icaweights * EEG.icasphere;
    else
        W = eye(62, 62);
    end
    
    disp('Splitting the data');
    load([folder '\epoch_inds.mat']);
    % splitting the data itself.
    clearing_mats = split_data(clearing_inds, data);
    thinking_mats = split_data(thinking_inds, data);
    num_e = length(clearing_mats);

    epoch_data.thinking_mats = thinking_mats;
    epoch_data.clearing_mats = clearing_mats;

    % computing the EEG features.
    disp('Computing the various features');
    [eeg_features, feature_labels] = get_eeg_features(epoch_data, W);

    % Getting the prompts.
    prompts = get_prompts(labels_fn);

    % Creating and saving the struct.
    disp('saving the struct');
    all_features = struct();
    all_features.eeg_features = eeg_features;
    all_features.feature_labels = feature_labels;
    all_features.prompts = prompts;
    if ICA
        save([folder '\all_features_ICA.mat'], 'all_features');
    else
        save([folder '\all_features_noICA.mat'], 'all_features','-v7.3');
    end
end
