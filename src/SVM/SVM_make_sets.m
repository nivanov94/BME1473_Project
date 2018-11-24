clear all
clc

all_feature = [];
all_class = [];

for j = 1:14
    load(['Features P' num2str(j), '.mat']);
    feature_new = feature;
    all_feature = [all_feature feature_new];
    load(['Class P' num2str(j), '.mat']);
    class_new = class_num;
    all_class = [all_class class_new];
end

%% Here you can define the criteria to select the traning sets

% Version 1, pick 0 - 75% as traning, 75% - 100% as test
num_trial = size(all_class,2);

training_feature = all_feature(:,1:round(0.75*num_trial))';
training_class = all_class(:,1:round(0.75*num_trial))';

test_feature = all_feature(:,round(0.75*num_trial)+1:num_trial)';
test_class = all_class(:,round(0.75*num_trial)+1:num_trial)';

%% a quick SVM for benchmarking 
Model = svm.train(training_feature,training_class);
Predict = svm.predict(Model,test_feature);
Accuracy = mean(test_class==Predict)*100;
fprintf('\nAccuracy =%d\n',Accuracy)

