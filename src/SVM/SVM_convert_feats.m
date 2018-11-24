clear all
clc

% The folders containing feature files.
folders = {'P02' 'MM05' 'MM08' 'MM09' 'MM10' 'MM11' 'MM12' 'MM14' 'MM15' 'MM16' 'MM18' 'MM19' 'MM20' 'MM21'};
data_path = 'F:\BME1473\Project\downsampled_data';

for j = 1:length(folders)
    clear class_num
    clear feature
    
    disp(['Convert extrated features for folder ' folders{j}]);
    folder = [data_path '\' folders{j}];
%     D = dir([folder '\*.mat']);
    set_fn = [folder '\' 'all_features_noICA-downsampled-baisc.mat'];

% Read the extracted features
    load(set_fn);
    think_features = all_features.eeg_features.thinking_feats;
    feature_labels = all_features.feature_labels;
    if size(feature_labels,2) ~= size(think_features{1,1},2)
        disp('Features do not match with labels');
        disp('Please check');
        exit
    end
    %Convert words to numbers
    class_word = all_features.prompts;
    for i = 1:length(class_word)
        switch class_word{1,i}
            case '/iy/'
                class_num(1,i) = 1;
            case '/uw/'
                class_num(1,i) = 2;
            case '/piy/'
                class_num(1,i) = 3;
            case '/tiy/'
                class_num(1,i) = 4;
            case '/diy/'
                class_num(1,i) = 5;
            case '/m/'
                class_num(1,i) = 6;
            case '/n/'
                class_num(1,i) = 7;
            case 'pat'
                class_num(1,i) = 8;
            case 'pot'
                class_num(1,i) = 9;
            case 'knew'
                class_num(1,i) = 10;
            case 'gnaw'
                class_num(1,i) = 11;
            otherwise
                disp('unknown class exists');
                class_num(1,i) = 0;
        end
    end
    %% Reshape the extracted features
    % Rows of TRAINING correspond to observations; columns correspond to features.
    num_trials = size(think_features,2);
    num_channels = size(think_features{1,num_trials},1);
    num_features_per_channels = size(think_features{1,num_trials},2);
    for r = 1:num_trials % r is the number of trials
        current_array = think_features{1,r};
        for c = 1:num_channels % c is the number of channels, should be 62 here
            for f = 1:num_features_per_channels % f is the number of features in each channels
                feature(f+c*num_channels,r) = current_array(c,f);
            end
        end
    end
    %% Save as mat
    disp('saving the converted features and labels');
    save(['Features P' num2str(j), '.mat'], 'feature');
    save(['Class P' num2str(j) '.mat'], 'class_num');
end
    