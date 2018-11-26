%
% [features] = extract_features(data, params)
%
% input:
%    data   - 1xN matrix containing the audio or EEG samples for one channel.
%    params - EEG general parameters
% output:
%    features       - (1xF) matrix of features, where F is the number of 
%                     features per epoch
%    feature_labels - {F} strings describing each feature dimension
%
function [features] = extract_features(data, params)
    global feature_fft
    global feature_wavelet
    global feature_logfft
    global feature_spf
    global feature_ar
    global feature_dwt
    global feature_moment
    global feature_central_moment
    
    features = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters used
    datatype       = params.data_type;
    wintype        = params.feature.wintype;
    samplefreq     = params.fs;

    fft_enable    = params.feature.fft.enable;
    logfft_enable = params.feature.logfft.enable;

    mean_enable        = params.feature.mean.enable;
    abs_mean_enable    = params.feature.abs_mean.enable;
    max_enable         = params.feature.max.enable;
    abs_max_enable     = params.feature.abs_max.enable;
    min_enable         = params.feature.min.enable;
    abs_min_enable     = params.feature.abs_min.enable;
    max_min_sum_enable = params.feature.max_min_sum.enable;
    max_min_dif_enable = params.feature.max_min_dif.enable;

    std_enable    = params.feature.std.enable;
    var_enable    = params.feature.var.enable;
    skew_enable   = params.feature.skew.enable;
    kurt_enable   = params.feature.kurt.enable;
    sum_enable    = params.feature.sum.enable;
    median_enable = params.feature.median.enable;
    
    %% Nick's parameters
    % AR
    ar_enable = params.feature.ar.enable;
    ar_order = params.feature.ar.order;
    ar_option = params.feature.ar.option;
    ar_band = params.feature.ar.band;
    % DWT
    dwt_enable = params.feature.dwt.enable;
    dwt_option = params.feature.dwt.option;
    dwt_wavelet = params.feature.dwt.wavelet;
    dwt_level = params.feature.dwt.level;
    % central moment
    central_moment_enable = params.feature.central_moment.enable;
    central_moment_values = params.feature.central_moment.values;
    % Moment
    moment_enable = params.feature.moment.enable;
    moment_values = params.faeture.moment.values;
    % SPF
    spf_enable = params.feature.spf.enable;
    spf_values = params.feature.spf.values;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract features
            
    current_data = data;
    epoch_features = [];
    
    % Computing the window size so that it corresponds to roughly 10 percent
    % of the data.
    window_size = ceil(size(current_data, 2) / 10);
    
    if mod(window_size, 2) == 1
        window_size = window_size + 1;
    end
    overlap_size = window_size / 2;

    % Create the overlapping frames (dimensions: numframes x framesize)
    x_framed = buffer(current_data, window_size, overlap_size, 'nodelay').';

    numframes = size(x_framed, 1); % Number of window frames
    framesize = size(x_framed, 2); % Number of samples per window frame
    
    % Each frame corresponds to ~10% of the data, and they overlap by a half,
    % so there should be 19 frames.
    assert(numframes == 19);
    
    % Choose windowing function
    switch wintype
        case 'hamming'
        w = hamming( framesize );
        case 'hann'
        w = hann( framesize );
        case 'triang'
        w = triang( framesize );
        case 'rect'
        w = rectwin( framesize );
        otherwise
        w = rectwin( framesize );
    end

    % Window every frame
    for f=1:numframes
        x_framed(f,:) = x_framed(f,:) .* w';
    end

    if fft_enable
        % Calculate FFT from the frame
        NFFT        = 2^nextpow2(framesize);
        fftfeat     = abs(fft(x_framed, NFFT, 2));
        feature_fft = fftfeat(:, 1:(NFFT/2+1));
    end

    if logfft_enable
        % Calculate log FFT from the frame
        NFFT           = 2^nextpow2(framesize);
        logfftfeat     = abs(fft(x_framed, NFFT, 2));
        feature_logfft = log(logfftfeat(:, 1:(NFFT/2+1)));
    end
    
    % For every frame
    feature_ar = [];
    feature_dwt = [];
    feature_spf = [];
    feature_moment = [];
    feature_central_moment = [];
    for f= 1:numframes

        %% Nick's functions
        if ar_enable
           feature_ar(f,:) = getARFeatures(x_framed(f,:),ar_order,ar_option,ar_band,samplefreq)';
        end
        
        if moment_enable
           feature_moment(f,:) = getMomentFeatures(x_framed(f,:),moment_values)';
        end
        
        if central_moment_enable
           feature_central_moment(f,:) = getCentralMomentFeatures(x_framed(f,:),central_moment_values)';
        end
        
        if dwt_enable
           feature_dwt(f,:) = getDWTFeatures(x_framed(f,:),dwt_wavelet,dwt_level,dwt_option,samplefreq)';
        end
        
        if spf_enable
           feature_spf(f,:) = getSpectralPowerFeatures(x_framed(f,:),spf_values,samplefreq)';
        end
        
    end

    %% Assemble data
    if fft_enable
        epoch_features = [epoch_features reshape(feature_fft, 1, size(feature_fft,1)*size(feature_fft,2))];
    end

    if logfft_enable
        epoch_features = [epoch_features reshape(feature_logfft, 1, size(feature_logfft,1)*size(feature_logfft,2))];
    end

    if mean_enable
        feats = mean(x_framed', 1);
        epoch_features = [epoch_features feats];
    end

    if abs_mean_enable
        feats = mean(abs(x_framed'), 1);
        epoch_features = [epoch_features feats];
    end

    if max_enable
        feats = max(x_framed');
        epoch_features = [epoch_features feats];
    end

    if abs_max_enable
        feats = max(abs(x_framed'));
        epoch_features = [epoch_features feats];
    end

    if min_enable
        feats = min(x_framed');
        epoch_features = [epoch_features feats];
    end

    if abs_min_enable
        feats = min(abs(x_framed'));
        epoch_features = [epoch_features feats];
    end

    if max_min_sum_enable
        feats = (max(x_framed') + min(x_framed'));
        epoch_features = [epoch_features feats];
    end

    if max_min_dif_enable
        feats = (max(x_framed') - min(x_framed'));
        epoch_features = [epoch_features feats];
    end

    if std_enable
      feats = std(x_framed, 0, 2)';
      epoch_features = [epoch_features feats];
    end

    if var_enable
      feats = var(x_framed, 0, 2)';
      epoch_features = [epoch_features feats];
    end

    if skew_enable
      feats = skewness(x_framed, 0, 2)';
      epoch_features = [epoch_features feats];
    end

    if kurt_enable
      feats = kurtosis(x_framed, 0, 2)';
      epoch_features = [epoch_features feats];
    end

    if sum_enable 
      feats = sum(x_framed, 2)';
      epoch_features = [epoch_features feats];
    end

    if median_enable
      feats = median(x_framed, 2)';
      epoch_features = [epoch_features feats];
    end
%% Nick's function
% AR
    if ar_enable
        feats = reshape(feature_ar, 1, size(feature_ar,1)*size(feature_ar,2));
        epoch_features = [epoch_features feats];
    end
% MOM
    if moment_enable
        feats = reshape(feature_moment, 1, size(feature_moment,1)*size(feature_moment,2));
        epoch_features = [epoch_features feats];
    end
% CENTRAL MOM
    if central_moment_enable
        feats = reshape(feature_central_moment, 1, size(feature_central_moment,1)*size(feature_central_moment,2));
        epoch_features = [epoch_features feats];
    end
% DWT
    if dwt_enable
        feats = reshape(feature_dwt, 1, size(feature_dwt,1)*size(feature_dwt,2));
        epoch_features = [epoch_features feats];
    end
% SPF
    if spf_enable
        feats = reshape(feature_spf, 1, size(feature_spf,1)*size(feature_spf,2));
        epoch_features = [epoch_features feats];
    end
features = epoch_features;
end
