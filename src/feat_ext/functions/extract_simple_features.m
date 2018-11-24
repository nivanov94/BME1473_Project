% Similar to extract_features.m, except that these features are much more simple
% and I don't window the data. We also take as input the struct params, which
% we do not use.

function [features] = extract_simple_features(data, params)
    
    delta = get_delta(data, false);
    ddelta = get_delta(data, true);
    features = [get_features(data) get_features(delta) get_features(ddelta)];
end

% Given the facial information, compute a few simple features. Variable dat is
% just a simple vector.
function [features] = get_features(dat)
    
    features = zeros(1, 14);
    % Mean and absolute mean
    features(1:2) = [mean(dat) mean(abs(dat))];
    % min and max features
    features(3:6) = [max(dat) max(abs(dat)) min(dat) min(abs(dat))];
    features(7:8) = [max(dat) + min(dat), max(dat) - min(dat)];
    % sixth power
    features(9) = [mean(dat.^6)];
    % standard deviation and variance.
    features(10:11) = [std(dat) var(dat)];
    % kurtosis, skewness, and median
    features(12:15) = [kurtosis(dat) skewness(dat) sum(dat) median(dat)];
    % RQA features.
    % y = crqa(dat, 1, 1, 0.1, 2000, 1000, 'euc', 'nogui');
    % features = [features mean(y)];
end

% Given the vector of facial information for one of the dimensions, compute the
% delta.
function [delta] = get_delta(dat, double_delta)
    
    dat_copy = dat(3:end);
    delta = dat_copy - dat(1:end-2);
end