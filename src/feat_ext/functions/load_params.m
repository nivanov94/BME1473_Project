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
    params.simplefeature.type = 'fft'; % 'fft' / 'logfft' / 'mfcc'
    
    if strcmp(data_type, 'eeg')
        params.fs = 200; % The sampling rate of the EEG data.
        params.data_type = 'eeg';
    elseif strcmp(data_type, 'audio')
        params.fs = 16000; % The sampling rate of the audio data.
        params.data_type = 'audio';
    else
        disp('Shit, we have a problem');
        return;
    end

    % type of window for windowing : rect/hamming/hann/triang
    params.feature.wintype      = 'rect';

    params.feature.selection.test      = 'ttest'; % ttest / kl / cohen / all
    params.feature.selection.dimension = 100; % Top n-dimensions chosen

    % Number of histogram bins for KL divergence
    params.feature.selection.kl.numbins = 50;
    % 1 for using shrinkage covariance estimates, else 0
    params.feature.selection.cohen.shrinkage = 0;
    params.feature.selection.cohen.nboot     = 0; % Number of bootstrap resamples

    params.feature.fft.enable    = 0;
    params.feature.logfft.enable = 0;

    params.feature.rpde.enable  = 0;
    params.feature.rpde.m       = 10;   % embedding dimension
    params.feature.rpde.tau     = 0;    % embedding time delay
    params.feature.rpde.epsilon = 0.12; % recurrence neighborhood radius

    params.feature.mfcc.enable    = 0;
    params.feature.mfcc.dimension = 15; % dimension of mfcc features per channel
    params.feature.mfcc.pivot     = 100; % pivot frequency for mel-scale

    params.feature.lpc.enable    = 0;
    params.feature.lpc.order     = 12;
    params.feature.lpc.dimension = 10; % dimension of LPC features per channel

    params.feature.wavelet.enable = 0;
    params.feature.wavelet.levels = 4;      % # of levels of wavelet decomposition
    params.feature.wavelet.window = 'sym4'; % mother wavelet

    params.feature.mean.enable        = 1;
    params.feature.abs_mean.enable    = 1;
    params.feature.max.enable         = 1;
    params.feature.min.enable         = 1;
    params.feature.abs_max.enable     = 1;
    params.feature.abs_min.enable     = 1;
    params.feature.max_min_sum.enable = 1;
    params.feature.max_min_dif.enable = 1;

    params.feature.n400_energy_windows = [4 8 16 32 64 128];

    % See "Epileptic Seizure Prediction Using Hybrid Feature Selection Over Multiple Intracranial EEG Electrode Contacts: A Report of Four Patients" 
    % by D'Alessandro et al. for a description of some of these features 
    params.feature.ehf.enable                  = 0;
    params.feature.curve_length.enable         = 0;
    params.feature.energy.enable               = 0;
    params.feature.nonlinear_energy.enable     = 0;
    params.feature.avg_nonlinear_energy.enable = 0;
    params.feature.sixth_power.enable          = 0;
    
    if strcmp(data_type, 'face')
        % Because we don't really know the sampling rate, so computing these
        % features feels way too hacky.
        params.feature.integral.enable = 0;
        params.feature.spectral_entropy.enable = 0;
    else
        params.feature.spectral_entropy.enable     = 0;
        params.feature.integral.enable             = 0;
    end

    params.feature.std.enable    = 1;
    params.feature.var.enable    = 1;
    params.feature.skew.enable   = 1;
    params.feature.kurt.enable   = 1;
    params.feature.sum.enable    = 1;
    params.feature.median.enable = 1;

    params.feature.vel.enable   = 0;
    params.feature.accel.enable = 0;
    
    %% Nick's functions
    % AR
    params.feature.nick.ar.enable = 0;
    params.feature.nick.ar.order = 1;
    params.feature.nick.ar.option = 1;
    params.feature.nick.ar.band = [2; 30];
    % Moment
    params.feature.nick.mom.enable = 0;
    params.faeture.nick.mom.vector = [];
    % DWT
    params.feature.nick.dwt.enable = 0;
    params.feature.nick.dwt.option = 1;
    params.feature.nick.dwt.wavelet = 'sym4';
    params.feature.nick.dwt.level = 4;
    % SPF
    params.feature.nick.spf.enable = 0;
end