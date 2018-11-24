%
% [r_mat] = perform_regression(r_mat, train_data, test_data, classifier)
%
% inputs:
%   r_mat           - Matrix to store the results.
%   train_data      - Training data for the model.
%   train_labels    - Training labels.
%   test_data       - Test data for the model.
%   test_labels     - test_labels
%   i               - The subject number.
%   N               - Feature size number.
%
% outputs:
%   r_mat           - The same matrix, except populated with results.
%
function [r_mat] = perform_regression(r_mat, train_data, train_labels, test_data, test_labels, i, n, classifier)

    test_mels = test_labels{1};
    test_mfccs = test_labels{2};
    train_mels = train_labels{1};
    train_mfccs = train_labels{2};
    results_mels = zeros(size(test_mels));
    results_mfccs = zeros(size(test_mfccs));
        
    if strcmp(classifier, 'lr') || strcmp(classifier, 'lr_poly')
        % Using a simple linear regression model or one with a
        % second-degree polynomial basis function.

        for j=1:size(test_mels, 2)
            if strcmp(classifier, 'lr')
                mels_model = polyfitn(train_data, train_mels(:, j), 1);
                mfccs_model = polyfitn(train_data, train_mfccs(:, j), 1);
            else
                modelterms = eye(size(train_data, 2)) * 2;
                modelterms = [zeros(1, size(train_data, 2)); modelterms];
                mels_model = polyfitn(train_data, train_mels(:, j), modelterms);
                mfccs_model = polyfitn(train_data, train_mfccs(:, j), modelterms);
            end

            results_mels(:, j) = polyvaln(mels_model, test_data);
            results_mfccs(:, j) = polyvaln(mfccs_model, test_data);
        end
    elseif strcmp(classifier, 'nn')
        % Using a neural network.
        
        net = fitnet([round(n * 0.8) round(n * 0.8)]);
        [net_mels, ~] = train(net, train_data', train_mels');
        [net_mfccs, ~] = train(net, train_data', train_mfccs');
        results_mels = net_mels(test_data')';
        results_mfccs = net_mfccs(test_data')';
    end
    
    % Having computed the predictions, we now evaluate the models.
    RMSE_mels = sqrt(mean((results_mels - test_mels) .^ 2, 1));
    RMSE_mfccs = sqrt(mean((results_mfccs - test_mfccs) .^ 2, 1));

    mel_corr = zeros(1, size(test_mels, 2));
    mfcc_corr = zeros(1, size(test_mels, 2));
    for j=1:size(test_mels, 2)
        mel_corr(j) = corr(results_mels(:, j), test_mels(:, j));
        mfcc_corr(j) = corr(results_mfccs(:, j), test_mfccs(:, j));
    end
    
    % If we are not doing user-dependent stuff. In this case, we are just doing
    % cross validation and so we'll need to accumulate values.
    if i == -1
        r_mat(n, 3, 1) = r_mat(n, 3, 1) + mean(mel_corr);
        r_mat(n, 3, 2) = r_mat(n, 3, 2) + mean(mfcc_corr);
        r_mat(n, 1, 1) = r_mat(n, 1, 1) + mean(RMSE_mels);
        r_mat(n, 1, 2) = r_mat(n, 1, 2) + mean(RMSE_mfccs);
        r_mat(n, 2, 1) = r_mat(n, 2, 1) + mean(RMSE_mels ./ sqrt(mean(test_mels .^ 2, 1)) * 100);
        r_mat(n, 2, 2) = r_mat(n, 2, 2) + mean(RMSE_mfccs ./ sqrt(mean(test_mfccs .^ 2, 1)) * 100);
    else
        r_mat(i, n, 3, 1) = mean(mel_corr);
        r_mat(i, n, 3, 2) = mean(mfcc_corr);
        r_mat(i, n, 1, 1) = mean(RMSE_mels);
        r_mat(i, n, 1, 2) = mean(RMSE_mfccs);
        r_mat(i, n, 2, 1) = mean(RMSE_mels ./ sqrt(mean(test_mels .^ 2, 1)) * 100);
        r_mat(i, n, 2, 2) = mean(RMSE_mfccs ./ sqrt(mean(test_mfccs .^ 2, 1)) * 100);
    end
end