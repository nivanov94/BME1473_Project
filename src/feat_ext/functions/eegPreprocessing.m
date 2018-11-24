% Function to perform preprocessing.
function [EEG] = eegPreprocessing(cnt_fn, params)
    
    % location of the channel locations file.
    chan_loc_fn = params.chan_loc_fn;

    EEG = pop_loadcnt(cnt_fn);
    % First removing the non-data channels: M1, M2, EMG, EKG, and Trigger.
    EEG.data = EEG.data(1:66, :);
    EEG.data(33, :) = [];
    EEG.data(43, :) = [];
    EEG.nbchan = 64;
    
    % Performing eye-movement correction.
    if params.occular
        EEG = pop_lms_regression(EEG, params.occular_chans, 3, 1e-10, []);
        
        % Removing the occular channels and adding the channel location file.
        EEG.data(params.occular_chans, :) = [];
        EEG.nbchan = 62;
        eloc = readlocs(chan_loc_fn);
        EEG.chanlocs = eloc;
    end
    
    % Performing band-pass filtering.
    EEG = pop_eegfilt(EEG, params.bandpass.low, 0);
    EEG = pop_eegfilt(EEG, 0, params.bandpass.high);

    % Removing the mean values from the channels.
    data = EEG.data';
    data = detrend(data, 'constant');
    EEG.data = data';

    % Performing spacial filtering.
    if params.laplacian
        data = laplacian_filter(EEG.data, EEG.chanlocs, params.nbrs);
        EEG.data = data;
    end

    % Perform ICA
    if params.ica
        EEG = pop_runica(EEG, 'chanind', 1:size(EEG.data, 1));
    end
end
