%
% [params] = load_params(data_type)
%
% input:
%   data_type   - A string denoting the type of data: EEG or acoustic
% output:
%   params      - The parameters used to train the features for both the EEG
%               data and speech data.
function [params] = load_params(data_type)

    % The parameters sections.
    params = struct();

    % additional feature selection options
    params.ica_used = 0;
    params.windowed = 1;
    params.num_electrodes = 62;
 
    params.simplefeature.type = 'fft'; % 'fft' / 'logfft' / 'mfcc'
    
    if strcmp(data_type, 'eeg')
        params.fs = 200; % The sampling rate of the EEG data.
        params.data_type = 'eeg';
    else
        disp('we have a problem');
        return;
    end

    % type of window for windowing : rect/hamming/hann/triang
    params.feature.wintype      = 'hamming';

    params.feature.selection.dimension = 100; % Top n-dimensions chosen

    % FFT based features
    params.feature.fft.enable    = 0;
    params.feature.logfft.enable = 0;

    % Simple measurement features
    params.feature.mean.enable        = 1;
    params.feature.abs_mean.enable    = 1;
    params.feature.max.enable         = 1;
    params.feature.min.enable         = 1;
    params.feature.abs_max.enable     = 1;
    params.feature.abs_min.enable     = 1;
    params.feature.max_min_sum.enable = 1;
    params.feature.max_min_dif.enable = 1;
    params.feature.std.enable         = 1;
    params.feature.var.enable         = 1;
    params.feature.skew.enable        = 1;
    params.feature.kurt.enable        = 1;
    params.feature.sum.enable         = 1;
    params.feature.median.enable      = 1;

    % AR features
    params.feature.ar.enable = 0;
    params.feature.ar.order = 3;
    params.feature.ar.option = 1;
    params.feature.ar.band = [[1:49]', [2:50]'];

    % Moment features
    params.feature.moment.enable = 0;
    params.faeture.moment.values = [1:10];

    % Central moment features
    params.feature.central_moment.enable = 1;
    params.feature.central_moment.values = [1:10];

    % DWT features
    params.feature.dwt.enable = 0;
    params.feature.dwt.option = [2:6];
    params.feature.dwt.wavelet = 'db4';
    params.feature.dwt.level = 5;

    % SPF features
    params.feature.spf.enable = 1;
    params.feature.spf.values = [1:50];

    % CSP features
    params.feature.csp.enable = 1;
    params.feature.csp.bands  = [[0.1, 4:4:32]' [4:4:36]'];
    params.feature.csp.m      = 1;
end
