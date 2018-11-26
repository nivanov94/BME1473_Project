clear all; clc;

% Add the feature extraction functions to the path 
addpath(genpath('D:\School\GRAD\BME1473\Project\code\src\feat_ext\functions'));
addpath(genpath('D:\School\GRAD\BME1473\Project\code\src'));
% Add the EEGLAB functions to the path
addpath(genpath('C:\Users\Nick\Documents\MATLAB\eeglab\eeglab14_1_2b'));
addpath(genpath('C:\Users\Nick\Documents\MATLAB\eeglab\eeglab14_1_2b\functions\popfunc'));
rmpath C:\Users\Nick\Documents\MATLAB\eeglab\eeglab14_1_2b\functions\octavefunc
rmpath C:\Users\Nick\Documents\MATLAB\eeglab\eeglab14_1_2b\functions\octavefunc\optim
rmpath C:\Users\Nick\Documents\MATLAB\eeglab\eeglab14_1_2b\functions\octavefunc\signal

% The folders for which we want to compute the features.
data_path = 'D:\School\GRAD\BME1473\Project\data';
misc_folder = [data_path '\misc'];

data_set = {'P02' 'MM05' 'MM08' 'MM09' 'MM10' 'MM11' 'MM12' 'MM14' 'MM15' 'MM16' 'MM18' 'MM19' 'MM20' 'MM21'};
subjects = size(data_set,2);

folds = 10;

train_set_sz = 8;
val_set_sz = 3;
test_set_sz = 3;

training_set = cell(1,folds);
test_set = cell(1,folds);
validation_set = cell(1,folds);

for k = 1:folds

  tr_set = cell(1,train_set_sz);
  te_set = cell(1,test_set_sz);
  va_set = cell(1,val_set_sz);

  train_start_index = floor((k-1)*subjects/folds);

  for i = 0:(subjects-1)
    if i < train_set_sz
      tr_set{i+1} = data_set{mod((train_start_index + i),subjects)+1};
    elseif (i - train_set_sz) < val_set_sz
      va_set{i+1 - train_set_sz} = data_set{mod((train_start_index + i),subjects)+1};
    else
      te_set{i+1 - train_set_sz - val_set_sz} = data_set{mod((train_start_index + i),subjects)+1};
    end
  end

  training_set{k}   = tr_set;
  validation_set{k} = va_set;
  test_set{k}       = te_set;
end

% clean up variables
clear tr_set; clear te_set; clear va_set; clear k; clear i; clear train_start_index;


ICA = false;
CSP = false;

params = load_params('eeg');

if CSP
    % learn the common spatial pattern projection matrices here
    % these need to be calculated separately because they are
    % dependent on which values are in the training set
    
    empty = true;

    for j=1:length(data_set)
        
        disp(['Calcuating the class covariance matrices for subject ' data_set{j}]);
        folder = [data_path '\' data_set{j}];
        if ~ICA
            D = dir([folder '\' data_set{j} '.set']);
        else
            D = dir([folder '\' data_set{j} '-ICA.set']);
        end

        kinect_folder = [folder '\kinect_data'];
        labels_fn = [kinect_folder '\labels.txt'];

        % split the file, create the epochs, and save.
        EEG = pop_loadset(D(1).name,folder);
        data = EEG.data;
        
        % Applying the ICA decomposition.
        if ICA
            ICA_Xform = EEG.icaweights * EEG.icasphere;
        else
            ICA_Xform = eye(62, 62);
        end

        disp('Splitting the data');
        load([folder '\epoch_inds.mat']);
        % splitting the data itself.
        clearing_mats = split_data(clearing_inds, data);
        thinking_mats = split_data(thinking_inds, data);
        num_e = length(clearing_mats);

        % Getting the prompts.
        prompts = get_prompts(labels_fn);
        num_prompts = [getNumLabel(prompts) zeros(1,num_e)]; % add zeros for the clearing samples

        unique_classes = unique(num_prompts);
        filters = params.feature.csp.bands;

        % computing the EEG features.
        disp('Computing the class covariance matrices');

        % combine the thinking and clearing data
        all_trial_eegs = [thinking_mats clearing_mats];
        % garbage collection
        clear clearing_mats; clear thinking_mats; clear EEG; clear data;
 
        if empty
            disp('Allocating composite covariance matrices');
            tot_cc = zeros(folds,length(unique_classes),size(filters,1),params.num_electrodes,params.num_electrodes);
            tot_c1 = zeros(folds,length(unique_classes),size(filters,1),params.num_electrodes,params.num_electrodes);
            empty = false;
        end

        C1 = zeros(length(unique_classes),size(filters,1),params.num_electrodes,params.num_electrodes);
        C2 = zeros(length(unique_classes),size(filters,1),params.num_electrodes,params.num_electrodes);
        for class = 1:length(unique_classes)
            disp(['  class: ',num2str(class)]);
            for band = 1:size(filters,1)
                disp(['    filter: ',num2str(band)]);
                [C1(class,band,:,:), C2(class,band,:,:)] = getCSPCompCovMat(all_trial_eegs,...
                                                                 num_prompts==unique_classes(class),...
                                                                 filters(band,:)/(params.fs/2),...
                                                                 ICA_Xform);
            end
        end
        % garbage collection
        clear all_trial_eegs;

        % apply the matrices to the folds where this subject is in the training set
        for k = 1:folds
            set_subjects = training_set{k};
            for i = 1:length(set_subjects)
                if strcmp(data_set{j},set_subjects{i})
                    tot_c1(k,:,:,:,:) = squeeze(tot_c1(k,:,:,:,:)) + C1;
                    tot_cc(k,:,:,:,:) = squeeze(tot_cc(k,:,:,:,:)) + C1 + C2;
                    break;
                end
            end
        end
    end
    % garbage collection
    clear C1; clear C2; clear set_subjects;

    % Calculate all the CSP projection matrices for each fold
    Wcsp = cell(1,folds);
    for k = 1:folds
        fold_Wcsp = cell(length(unique_classes),size(filters,1));
        for class = 1:length(unique_classes)
            for band = 1:size(filters,1)
                % calculate Wcsp matrix
                Wcsp_mat = getCSPProjMat(squeeze(tot_cc(k,class,band,:,:)),squeeze(tot_c1(k,class,band,:,:)));
                % save to appropriate cell
                fold_Wcsp{class,band} = Wcsp_mat;
             end
         end
         % save these matrices for this fold
         Wcsp{k} = fold_Wcsp;
    end
    % clean up variables
    clear tot_cc; clear tot_c1; clear Wcsp_mat; clear fold_Wcsp;


    % calculate the features for each trial using the projection matrices for each fold
    for j=1:length(data_set)
        folder = [data_path '\' data_set{j}];
        if ~ICA
            D = dir([folder '\' data_set{j} '.set']);
        else
            D = dir([folder '\' data_set{j} '-ICA.set']);
        end

        % split the file, create the epochs, and save.
        EEG = pop_loadset(D(1).name,folder);
        data = EEG.data;
        
        % Applying the ICA decomposition.
        if ICA
            ICA_Xform = EEG.icaweights * EEG.icasphere;
        else
            ICA_Xform = eye(62, 62);
        end

        disp('Splitting the data');
        load([folder '\epoch_inds.mat']);
        % splitting the data itself.
        clearing_mats = split_data(clearing_inds, data);
        thinking_mats = split_data(thinking_inds, data);
        num_e = length(clearing_mats);

        unique_classes = unique(num_prompts);
        filters = params.feature.csp.bands;

        % computing the EEG features.
        % combine the thinking and clearing data
        all_trial_eegs = [thinking_mats clearing_mats];

        CSP_feats = cell(1,folds);
        for k = 1:folds
            fold_Wcsp = Wcsp{k};
            trial_csp_feats = cell(1,length(all_trial_eegs));
            disp(['fold: ' num2str(k)]);
            for t = 1:length(all_trial_eegs)
                disp(['  trial ' num2str(t)]);
                trial_eeg = ICA_Xform*double(cell2mat(all_trial_eegs(t))); % note application of ICA here
                class_band_feats = cell(length(unique_classes),size(filters,1));
                for class = 1:length(unique_classes)
                    for band = 1:size(filters,1)
                        class_band_feats{class,band} = getCSPFeatures(trial_eeg,...
                                                                      fold_Wcsp{class,band},...
                                                                      filters(band,:),...
                                                                      params.fs);
                    end
                end
                trial_csp_feats{t} = class_band_feats;
            end
            CSP_feats{k} = trial_csp_feats;
        end

        % save the data
        if ICA
             save([folder '\CSP_features_ICA.mat'],'CSP_feats');
        else
             save([folder '\CSP_features_noICA.mat'],'CSP_feats');
        end
    end
end


% calcualte the features that are independent of the training set
for j=1:length(data_set)
    disp(['Computing features for folder ' data_set{j}]);
    folder = [data_path '\' data_set{j}];
    if ~ICA
        D = dir([folder '\' data_set{j} '.set']);
    else
        D = dir([folder '\' data_set{j} '-ICA.set']);
    end
    % there should only be 1 set file.
    %set_path = [folder '\'];
    %set_fn = [folder '\' D(1).name];
    
    kinect_folder = [folder '\kinect_data'];
    labels_fn = [kinect_folder '\labels.txt'];

    % split the file, create the epochs, and save.
    EEG = pop_loadset(D(1).name,folder);
    data = EEG.data;
    
    % Applying the ICA decomposition.
    if ICA
        ICA_Xform = EEG.icaweights * EEG.icasphere;
    else
        ICA_Xform = eye(62, 62);
    end
    
    disp('Splitting the data');
    load([folder '\epoch_inds.mat']);
    % splitting the data itself.
    clearing_mats = split_data(clearing_inds, data);
    thinking_mats = split_data(thinking_inds, data);
    num_e = length(clearing_mats);

    epoch_data.thinking_mats = thinking_mats;
    epoch_data.clearing_mats = clearing_mats;

    % Getting the prompts.
    prompts = get_prompts(labels_fn);
    num_prompts = getNumLabel(prompts);

    % computing the EEG features.
    disp('Computing the various features');
    [eeg_features, feature_labels] = get_eeg_features(epoch_data, ICA_Xform);

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
