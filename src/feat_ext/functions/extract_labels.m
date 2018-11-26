%
% [feature_labels] = extract_labels(params, numframes)
%
% input:
%    params - EEG general parameters
% output:
%    feature_labels - {F} strings describing each feature dimension
%

function [feature_labels] = extract_labels(params, numframes)
    global feature_fft
    global feature_logfft
    global feature_spf
    global feature_ar
    global feature_dwt
    global feature_moment
    global feature_central_moment

    
    feature_labels = {};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters used
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
    % Create labels 
    
    % Add all FFT labels
    if fft_enable
        for d=1:size(feature_fft,2)
            for w=1:numframes
                feature_labels{end+1} = ['FFT:W', int2str(w), ',D', int2str(d)];
            end
        end
    end

    % Add all log FFT labels
    if logfft_enable
        for d=1:size(feature_logfft,2)
          for w=1:numframes
            feature_labels{end+1} = ['logFFT:W', int2str(w), ',D', int2str(d)];
          end
        end
    end

    % Add mean window value label
    if mean_enable
        for w=1:numframes
          feature_labels{end+1} = ['Mean:W', int2str(w)];
        end
    end

    % Add absolute mean window value label
    if abs_mean_enable
        for w=1:numframes
          feature_labels{end+1} = ['Absmean:W', int2str(w)];
        end
    end

    % Add maximum window value label
    if max_enable
        for w=1:numframes
          feature_labels{end+1} = ['Max:W', int2str(w)];
        end
    end

    % Add absolute maximum window value label
    if abs_max_enable
        for w=1:numframes
          feature_labels{end+1} = ['Absmax:W', int2str(w)];
        end
    end

    % Add minimum window value label
    if min_enable
        for w=1:numframes
          feature_labels{end+1} = ['Min:W', int2str(w)];
        end
    end

    if abs_min_enable
        for w=1:numframes
          feature_labels{end+1} = ['Absmin:W', int2str(w)];
        end
    end

    % Add minimum maximum sum window value label
    if max_min_sum_enable
        for w=1:numframes
          feature_labels{end+1} = ['Min+Max:W', int2str(w)];
        end
    end

    % Add maximum-minimum difference window value label
    if max_min_dif_enable
        for w=1:numframes
          feature_labels{end+1} = ['Max-Min:W', int2str(w)];
        end
    end

    % Add standard deviation window value label
    if std_enable
        for w=1:numframes
          feature_labels{end+1} = ['STD:W', int2str(w)];
        end
    end

    % Add variance window value label
    if var_enable
        for w=1:numframes
          feature_labels{end+1} = ['VAR:W', int2str(w)];
        end
    end

    % Add skewness window value label
    if skew_enable
        for w=1:numframes
          feature_labels{end+1} = ['SKEW:W', int2str(w)];
        end
    end

    % Add kurtosis window value label
    if kurt_enable
        for w=1:numframes
          feature_labels{end+1} = ['KURT:W', int2str(w)];
        end
    end

    % Add sum window value label
    if sum_enable
        for w=1:numframes
          feature_labels{end+1} = ['Sum:W', int2str(w)];
        end
    end

    % Add median window value label
    if median_enable
        for w=1:numframes
          feature_labels{end+1} = ['MEDIAN:W', int2str(w)];
        end
    end
    
    %% Nick's functions
    if ar_enable
        for d=1:size(feature_ar,2)
          for w=1:numframes
            feature_labels{end+1} = ['AR:W', int2str(w), ',D', int2str(d)];
          end
        end
    end
    
    if moment_enable
        for d=1:size(feature_moment,2)
          for w=1:numframes
            feature_labels{end+1} = ['MOMENT:W', int2str(w), ',D', int2str(d)];
          end
        end
    end
    
    if central_moment_enable
        for d=1:size(feature_central_moment,2)
          for w=1:numframes
            feature_labels{end+1} = ['CENTRALMOMENT:W', int2str(w), ',D', int2str(d)];
          end
        end
    end
    
    if dwt_enable
        for d=1:size(feature_dwt,2)
          for w=1:numframes
            feature_labels{end+1} = ['DWT:W', int2str(w), ',D', int2str(d)];
          end
        end
    end
    
    if spf_enable
        for d=1:size(feature_spf,2)
          for w=1:numframes
            feature_labels{end+1} = ['SPF:W', int2str(w), ',D', int2str(d)];
          end
        end
    end
end