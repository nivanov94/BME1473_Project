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
    global feature_mom
    
    features = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters used
    datatype       = params.data_type;
    wintype        = params.feature.wintype;
    samplefreq     = params.fs;

    fft_enable    = params.feature.fft.enable;
    logfft_enable = params.feature.logfft.enable;

    rpde_enable  = params.feature.rpde.enable;
    rpde_m       = params.feature.rpde.m;
    rpde_tau     = params.feature.rpde.tau;
    rpde_epsilon = params.feature.rpde.epsilon;

    mfcc_enable    = params.feature.mfcc.enable;
    mfcc_dimension = params.feature.mfcc.dimension;
    mfcc_pivot     = params.feature.mfcc.pivot;

    lpc_enable    = params.feature.lpc.enable;
    lpc_dimension = params.feature.lpc.dimension;
    lpc_order     = params.feature.lpc.order;

    wavelet_enable = params.feature.wavelet.enable;
    wavelet_levels = params.feature.wavelet.levels;
    wavelet_window = params.feature.wavelet.window;

    mean_enable        = params.feature.mean.enable;
    abs_mean_enable    = params.feature.abs_mean.enable;
    max_enable         = params.feature.max.enable;
    abs_max_enable     = params.feature.abs_max.enable;
    min_enable         = params.feature.min.enable;
    abs_min_enable     = params.feature.abs_min.enable;
    max_min_sum_enable = params.feature.max_min_sum.enable;
    max_min_dif_enable = params.feature.max_min_dif.enable;

    n400_energy_windows = params.feature.n400_energy_windows;

    ehf_enable                  = params.feature.ehf.enable;
    curve_length_enable         = params.feature.curve_length.enable;
    energy_enable               = params.feature.energy.enable;
    nonlinear_energy_enable     = params.feature.nonlinear_energy.enable;
    avg_nonlinear_energy_enable = params.feature.avg_nonlinear_energy.enable;
    spectral_entropy_enable     = params.feature.spectral_entropy.enable;
    sixth_power_enable          = params.feature.sixth_power.enable;
    integral_enable             = params.feature.integral.enable;

    std_enable    = params.feature.std.enable;
    var_enable    = params.feature.var.enable;
    skew_enable   = params.feature.skew.enable;
    kurt_enable   = params.feature.kurt.enable;
    sum_enable    = params.feature.sum.enable;
    median_enable = params.feature.median.enable;

    vel_enable   = params.feature.vel.enable;
    accel_enable = params.feature.accel.enable;
    
    %% Nick's parameters
    % AR
    nick_ar_enable = params.feature.nick.ar.enable;
    nick_ar_order = params.feature.nick.ar.order;
    nick_ar_option = params.feature.nick.ar.option;
    nick_ar_band = params.feature.nick.ar.band;
    % Moment
    nick_mom_enable = params.feature.nick.mom.enable;
    nick_mom_vector = params.faeture.nick.mom.vector;
    % DWT
    nick_dwt_enable = params.feature.nick.dwt.enable;
    nick_dwt_option = params.feature.nick.dwt.option;
    nick_dwt_wavelet = params.feature.nick.dwt.wavelet;
    nick_dwt_level = params.feature.nick.dwt.level;
    % SPF
    nick_spf_enable = params.feature.nick.spf.enable;
    
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

    if mfcc_enable
        feature_mfcc = melcepst(current_data, samplefreq, 'e0', mfcc_dimension, floor(3*log(samplefreq)), window_size, window_size - overlap_size );
    end

    if lpc_enable
        [lpc_ar, lpc_e, lpc_k] = lpcauto(current_data, lpc_order, [window_size, window_size, 0]);
        feature_lpc = [lpc_ar, lpc_e, lpc_k];
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
    
    %% Nick's function (non-frame-based)
    

    % For every frame
    feature_rpde    = [];
    feature_wavelet = [];
    feature_ar = [];
    feature_mom = [];
    feature_dwt = [];
    feature_spf = [];
    for f= 1:numframes

        if rpde_enable
            % Calculate RPDE
            [H, rpd] = rpde(double(x_framed(f,:))', rpde_m, rpde_tau, rpde_epsilon);
            H(isnan(H)) = 0;
            rpd(isnan(rpd)) = 0;
            feature_rpde(f,:) = [H, rpd'];
        end

        if wavelet_enable
            % Calculate wavelet
            [feature_wavelet(f,:), ~] = wavedec(double(x_framed(f,:)), wavelet_levels, wavelet_window );
        end
        
        %% Nick's functions
        if nick_ar_enable
           [feature_ar(f,:),~] = getARFeatures(x_framed(f,:),nick_ar_order,nick_ar_option,nick_ar_band,samplefreq)';
        end
        
        if nick_mom_enable
           feature_mom(f,:) = getMomentFeatures(x_framed(f,:),nick_mom_vector)';
        end
        
        if nick_dwt_enable
           feature_dwt(f,:) = getDWTFeatures(x_framed(f,:),nick_dwt_wavelet,nick_dwt_level,nick_dwt_option,samplefreq)';
        end
        
        if nick_spf_enable
            % Calculate spf
            feature_spf(f,:) = (10*log10(pwelch(x_framed(f,:),[],[],2,samplefreq)))';
        end
        
    end

    %% Assemble data
    if fft_enable
        epoch_features = [epoch_features reshape(feature_fft, 1, size(feature_fft,1)*size(feature_fft,2))];
        if vel_enable || accel_enable
            vels = compute_delta(feature_fft);
            if vel_enable
                epoch_features = [epoch_features reshape(vels, 1, size(vels,1)*size(vels,2))];
            end
            if accel_enable
                accels = compute_delta(vels);
                epoch_features = [epoch_features reshape(accels, 1, size(accels,1)*size(accels,2))];
            end
        end
    end

    if logfft_enable
        epoch_features = [epoch_features reshape(feature_logfft, 1, size(feature_logfft,1)*size(feature_logfft,2))];

        if vel_enable || accel_enable
            vels = compute_delta(feature_logfft);

            if vel_enable
                epoch_features = [epoch_features reshape(vels, 1, size(vels,1)*size(vels,2))];
            end
            if accel_enable
                accels = compute_delta(vels);
                epoch_features = [epoch_features reshape(accels, 1, size(accels,1)*size(accels,2))];
            end
        end
    end

    if rpde_enable
        epoch_features = [epoch_features reshape(feature_rpde, 1, size(feature_rpde,1)*size(feature_rpde,2))];

        if vel_enable || accel_enable
            vels = compute_delta(feature_rpde);

            if vel_enable
                epoch_features = [epoch_features reshape(vels, 1, size(vels,1)*size(vels,2))];
            end
            if accel_enable
                accels = compute_delta(vels);
                epoch_features = [epoch_features reshape(accels, 1, size(accels,1)*size(accels,2))];
            end
        end
    end

    if mfcc_enable
        epoch_features = [epoch_features reshape(feature_mfcc, 1, size(feature_mfcc,1)*size(feature_mfcc,2))];

        if vel_enable || accel_enable
            vels = compute_delta(feature_mfcc);

            if vel_enable
                epoch_features = [epoch_features reshape(vels, 1, size(vels,1)*size(vels,2))];
            end
            if accel_enable
                accels = compute_delta(vels);
                epoch_features = [epoch_features reshape(accels, 1, size(accels,1)*size(accels,2))];
            end
        end
    end

    if lpc_enable
        epoch_features = [epoch_features reshape(feature_lpc, 1, size(feature_lpc,1)*size(feature_lpc,2))];

        if vel_enable || accel_enable
            vels = compute_delta(feature_lpc);

            if vel_enable
                epoch_features = [epoch_features reshape(vels, 1, size(vels,1)*size(vels,2))];
            end
            if accel_enable
                accels = compute_delta(vels);
                epoch_features = [epoch_features reshape(accels, 1, size(accels,1)*size(accels,2))];
            end
        end
    end

    if wavelet_enable
        epoch_features = [epoch_features reshape(feature_wavelet, 1, size(feature_wavelet,1)*size(feature_wavelet,2))];

        if vel_enable || accel_enable
            vels = compute_delta(feature_wavelet);

            if vel_enable
                epoch_features = [epoch_features reshape(vels, 1, size(vels,1)*size(vels,2))];
            end
            if accel_enable
                accels = compute_delta(vels);
                epoch_features = [epoch_features reshape(accels, 1, size(accels,1)*size(accels,2))];
            end
        end
    end

    if mean_enable
        feats = mean(x_framed', 1);
        epoch_features = [epoch_features feats];

        if vel_enable || accel_enable
            vels = compute_delta(feats');

            if vel_enable
                epoch_features = [epoch_features vels'];
            end
            if accel_enable
                accels = compute_delta(vels);
                epoch_features = [epoch_features accels'];
            end
        end
    end

    if abs_mean_enable
        feats = mean(abs(x_framed'), 1);
        epoch_features = [epoch_features feats];

        if vel_enable || accel_enable
            vels = compute_delta(feats');

            if vel_enable
                epoch_features = [epoch_features vels'];
            end
            if accel_enable
                accels = compute_delta(vels);
                epoch_features = [epoch_features accels'];
            end
        end
    end

    if max_enable
        feats = max(x_framed');
        epoch_features = [epoch_features feats];

        if vel_enable || accel_enable
            vels = compute_delta(feats');

            if vel_enable
                epoch_features = [epoch_features vels'];
            end
            if accel_enable
                accels = compute_delta(vels);
                epoch_features = [epoch_features accels'];
            end
        end
    end

    if abs_max_enable
        feats = max(abs(x_framed'));
        epoch_features = [epoch_features feats];

        if vel_enable || accel_enable
            vels = compute_delta(feats');

            if vel_enable
                epoch_features = [epoch_features vels'];
            end
            if accel_enable
                accels = compute_delta(vels);
                epoch_features = [epoch_features accels'];
            end
        end
    end

    if min_enable
        feats = min(x_framed');
        epoch_features = [epoch_features feats];

        if vel_enable || accel_enable
            vels = compute_delta(feats');

            if vel_enable
                epoch_features = [epoch_features vels'];
            end
            if accel_enable
                accels = compute_delta(vels);
                epoch_features = [epoch_features accels'];
            end
        end
    end

    if abs_min_enable
        feats = min(abs(x_framed'));
        epoch_features = [epoch_features feats];

        if vel_enable || accel_enable
            vels = compute_delta(feats');

            if vel_enable
                epoch_features = [epoch_features vels'];
            end
            if accel_enable
                accels = compute_delta(vels);
                epoch_features = [epoch_features accels'];
            end
        end
    end

    if max_min_sum_enable
        feats = (max(x_framed') + min(x_framed'));
        epoch_features = [epoch_features feats];

        if vel_enable || accel_enable
            vels = compute_delta(feats');

            if vel_enable
                epoch_features = [epoch_features vels'];
            end
            if accel_enable
                accels = compute_delta(vels);
                epoch_features = [epoch_features accels'];
            end
        end
    end

    if max_min_dif_enable
        feats = (max(x_framed') - min(x_framed'));
        epoch_features = [epoch_features feats];

        if vel_enable || accel_enable
            vels = compute_delta(feats');

            if vel_enable
                epoch_features = [epoch_features vels'];
            end
            if accel_enable
                accels = compute_delta(vels);
                epoch_features = [epoch_features accels'];
            end
        end
    end

    if ehf_enable
        tmp = x_framed';
        feats = max(tmp) ./ ( 2 * sqrt(2 * sqrt(sum(tmp.^2) / size(tmp,1))) );
        epoch_features = [epoch_features feats];

        if vel_enable || accel_enable
            vels = compute_delta(feats');

            if vel_enable
                epoch_features = [epoch_features vels'];
            end
            if accel_enable
                accels = compute_delta(vels);
                epoch_features = [epoch_features accels'];
            end
        end
    end

    if curve_length_enable
      feats = sum(abs(x_framed - [zeros(numframes, 1) x_framed(:, 2:end)]), 2)';
      epoch_features = [epoch_features feats];

      if vel_enable || accel_enable
        vels = compute_delta(feats');

        if vel_enable
          epoch_features = [epoch_features vels'];
        end
        if accel_enable
          accels = compute_delta(vels);
          epoch_features = [epoch_features accels'];
        end
      end
    end

    if energy_enable
      feats = sum(x_framed.^2, 2)' / numframes;
      epoch_features = [epoch_features feats];

      if vel_enable || accel_enable
        vels = compute_delta(feats');

        if vel_enable
          epoch_features = [epoch_features vels'];
        end
        if accel_enable
          accels = compute_delta(vels);
          epoch_features = [epoch_features accels'];
        end
      end
    end

    if nonlinear_energy_enable || avg_nonlinear_energy_enable
      tmp = x_framed.^2 - [ones(numframes,1) x_framed(:,2:end)] .* [x_framed(:,1:(end-1)) ones(numframes,1)];
      if nonlinear_energy_enable
        feats = reshape(tmp, 1, numframes*framesize);
        epoch_features = [epoch_features feats];

        if vel_enable || accel_enable
          vels = compute_delta(feats');

          if vel_enable
            epoch_features = [epoch_features vels'];
          end
          if accel_enable
            accels = compute_delta(vels);
            epoch_features = [epoch_features accels'];
          end
        end
      end

      if avg_nonlinear_energy_enable
        feats = mean(tmp, 2)';
        epoch_features = [epoch_features feats];

        if vel_enable || accel_enable
          vels = compute_delta(feats');

          if vel_enable
            epoch_features = [epoch_features vels'];
          end
          if accel_enable
            accels = compute_delta(vels);
            epoch_features = [epoch_features accels'];
          end
        end
      end
    end

    if spectral_entropy_enable
      tmp = abs(fft(x_framed')').^2 * samplefreq / numframes;
      feats = -sum(tmp .* log2(tmp), 2)';
      epoch_features = [epoch_features feats];

      if vel_enable || accel_enable
        vels = compute_delta(feats');

        if vel_enable
          epoch_features = [epoch_features vels'];
        end
        if accel_enable
          accels = compute_delta(vels);
          epoch_features = [epoch_features accels'];
        end
      end
    end

    if sixth_power_enable
      feats = mean(x_framed.^6,2)';
      epoch_features = [epoch_features feats];

      if vel_enable || accel_enable
        vels = compute_delta(feats');

        if vel_enable
          epoch_features = [epoch_features vels'];
        end
        if accel_enable
          accels = compute_delta(vels);
          epoch_features = [epoch_features accels'];
        end
      end
    end

    if integral_enable
      feats = trapz(0:(1/samplefreq):((framesize-1)/samplefreq), x_framed');
      epoch_features = [epoch_features feats];

      if vel_enable || accel_enable
        vels = compute_delta(feats');

        if vel_enable
          epoch_features = [epoch_features vels'];
        end
        if accel_enable
          accels = compute_delta(vels);
          epoch_features = [epoch_features accels'];
        end
      end
    end

    if std_enable
      feats = std(x_framed, 0, 2)';
      epoch_features = [epoch_features feats];

      if vel_enable || accel_enable
        vels = compute_delta(feats');

        if vel_enable
          epoch_features = [epoch_features vels'];
        end
        if accel_enable
          accels = compute_delta(vels);
          epoch_features = [epoch_features accels'];
        end
      end
    end

    if var_enable
      feats = var(x_framed, 0, 2)';
      epoch_features = [epoch_features feats];

      if vel_enable || accel_enable
        vels = compute_delta(feats');

        if vel_enable
          epoch_features = [epoch_features vels'];
        end
        if accel_enable
          accels = compute_delta(vels);
          epoch_features = [epoch_features accels'];
        end
      end
    end

    if skew_enable
      feats = skewness(x_framed, 0, 2)';
      epoch_features = [epoch_features feats];

      if vel_enable || accel_enable
        vels = compute_delta(feats');

        if vel_enable
          epoch_features = [epoch_features vels'];
        end
        if accel_enable
          accels = compute_delta(vels);
          epoch_features = [epoch_features accels'];
        end
      end
    end

    if kurt_enable
      feats = kurtosis(x_framed, 0, 2)';
      epoch_features = [epoch_features feats];

      if vel_enable || accel_enable
        vels = compute_delta(feats');

        if vel_enable
          epoch_features = [epoch_features vels'];
        end
        if accel_enable
          accels = compute_delta(vels);
          epoch_features = [epoch_features accels'];
        end
      end
    end

    if sum_enable 
      feats = sum(x_framed, 2)';
      epoch_features = [epoch_features feats];

      if vel_enable || accel_enable
        vels = compute_delta(feats');

        if vel_enable
          epoch_features = [epoch_features vels'];
        end
        if accel_enable
          accels = compute_delta(vels);
          epoch_features = [epoch_features accels'];
        end
      end
    end

    if median_enable
      feats = median(x_framed, 2)';
      epoch_features = [epoch_features feats];

      if vel_enable || accel_enable
        vels = compute_delta(feats');

        if vel_enable
          epoch_features = [epoch_features vels'];
        end
        if accel_enable
          accels = compute_delta(vels);
          epoch_features = [epoch_features accels'];
        end
      end
    end
%% Nick's function
% AR
    if nick_ar_enable
        feats = reshape(feature_ar, 1, size(feature_ar,1)*size(feature_ar,2));
        epoch_features = [epoch_features feats];
    end
% MOM
    if nick_mom_enable
        feats = reshape(feature_mom, 1, size(feature_mom,1)*size(feature_mom,2));
        epoch_features = [epoch_features feats];
    end
% DWT
    if nick_dwt_enable
        feats = reshape(feature_dwt, 1, size(feature_dwt,1)*size(feature_dwt,2));
        epoch_features = [epoch_features feats];
    end
% SPF
    if nick_spf_enable
        feats = reshape(feature_spf, 1, size(feature_spf,1)*size(feature_spf,2));
        epoch_features = [epoch_features feats];
    end
features = epoch_features;
end
  
function [delta_feat] = compute_delta(windowed_feat)
    
    W = size(windowed_feat,1);
    D = size(windowed_feat,2);

    delta_feat = zeros(W, D);

    % Compute deltas
    delta_feat(1,:) = windowed_feat(2,:) - windowed_feat(1,:);
    for w=2:(W-1)
        delta_feat(w,:) = windowed_feat(w+1,:) - windowed_feat(w-1,:);
    end
    if W > 1
        delta_feat(end,:) = windowed_feat(end,:) - windowed_feat(end-1,:);
    end
end