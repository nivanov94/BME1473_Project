% [face_features] = get_face_features(kinect_folder)
%
% input:
%   kinect_folder       - The filename of the directory containing the kinect
%                       data.
% output:
%   face_features       - 1xN cell array, where element i contains a Dx6 
%                       matrix, where D is the number of data points.

function [face_features] = get_face_features(kinect_folder)
    
    % We look at the wav files, since they are complete, but for some folders,
    % there are missing 'face' files.
    D = dir([kinect_folder '/*.wav']);
    num_files = size(D, 1);
    face_features = cell(1, num_files);
    
    for i=1:num_files
        face_fn = [kinect_folder '/' num2str(i-1) '.AnimU'];
        if exist(face_fn) == 2
            data = dlmread(face_fn);
            features = [];
            % computing the features, if there is enough data.
            if size(data, 1) >= 5
                for j=1:size(data, 2)
                    dat = data(:, j)';
                    delta = get_delta(dat, false);
                    ddelta = get_delta(delta, true);
                    features = [features; [get_features(dat) get_features(delta) get_features(ddelta)]];
                end
            end
            face_features{i} = features;
        end
    end
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