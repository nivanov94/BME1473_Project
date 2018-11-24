clear all; clc; close all;

%% Split the data into training, test, and validation sets for 10-fold CV

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



%% Iterate over all folds
for k = 1:folds


end
