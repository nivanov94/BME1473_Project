%
% writeFeatures(feat_params)
%
% input:
%   feat_params     - A struct containing the information needed to write the
%                   features.
% output:
%   None. Just write the features to some file.
function writeFeatures(feat_params)
    
    % The feature indices.
    EEG_feats = feat_params.EEG_feats;
    WAV_feats = feat_params.WAV_feats;
    FAC_feats = feat_params.FAC_feats;
    EW_feats = feat_params.EW_feats;
    EF_feats = feat_params.EF_feats;
    ALL_feats = feat_params.ALL_feats;
    % Experiment information.
    t = feat_params.t;
    class_i = feat_params.class_i;
    logPath = feat_params.logPath;
    
    if class_i ~= 0
        fid = fopen([logPath 't' num2str(t) '_class' num2str(class_i) '.log'], 'w');
    else
        fid = fopen([logPath 't' num2str(t) '.log'], 'w');
    end
    
    fprintf(fid, 'Feature information\n');
    
    if ~isempty(EEG_feats)
        ep = feat_params.ep;
        EEG_labels = feat_params.EEG_labels;
        EEG_lab_chans = feat_params.EEG_lab_chans;
        fprintf(fid, '\nWriting out the EEG features\n');
        
        for i=1:length(EEG_feats) 
            fprintf(fid, ['Channel ' num2str(EEG_lab_chans(EEG_feats(i))) ' ' ...
            EEG_labels{EEG_feats(i)} 9 num2str(ep(i)) '\n']);
        end
    end
    
    if ~isempty(WAV_feats)
        wp = feat_params.wp;
        WAV_labels = feat_params.WAV_labels;
        fprintf(fid, '\nWriting out the WAV features\n');
        
        for i=1:length(WAV_feats) 
            fprintf(fid, [WAV_labels{WAV_feats(i)} 9 num2str(wp(i)) '\n']);
        end
    end
    
    if ~isempty(FAC_feats)
        fp = feat_params.fp;
        FAC_labels = feat_params.FAC_labels;
        FAC_lab_AUs = feat_params.FAC_lab_AUs;
        fprintf(fid, '\nWriting out the FAC features\n');
        
        for i=1:length(FAC_feats) 
            fprintf(fid, ['AU ' num2str(FAC_lab_AUs(FAC_feats(i))) ' ' ...
            FAC_labels{FAC_feats(i)} 9 num2str(fp(i)) '\n']);
        end
    end
    
    if ~isempty(EW_feats)
        EEG_labels = feat_params.EEG_labels;
        EEG_lab_chans = feat_params.EEG_lab_chans;
        WAV_labels = feat_params.WAV_labels;
        max_eeg_feat = length(EEG_labels);
        fprintf(fid, '\nWriting out the EEG and WAV features\n');
        
        for i=1:length(EW_feats)
            f = EW_feats(i);
            if f <= max_eeg_feat
                fprintf(fid, ['Channel ' num2str(EEG_lab_chans(f)) ' ' ...
                EEG_labels{f} '\n']);
            else
                f = f - max_eeg_feat;
                fprintf(fid, [WAV_labels{f} '\n']);
            end
        end
    end
    
    if ~isempty(EF_feats)
        EEG_labels = feat_params.EEG_labels;
        EEG_lab_chans = feat_params.EEG_lab_chans;
        FAC_labels = feat_params.FAC_labels;
        FAC_lab_AUs = feat_params.FAC_lab_AUs;
        max_eeg_feat = length(EEG_labels);
        fprintf(fid, '\nWriting out the EEG and FAC features\n');
        
        for i=1:length(EF_feats)
            f = EF_feats(i);
            if f <= max_eeg_feat
                fprintf(fid, ['Channel ' num2str(EEG_lab_chans(f)) ' ' ...
                EEG_labels{f} '\n']);
            else
                f = f - max_eeg_feat;
                fprintf(fid, ['AU ' num2str(FAC_lab_AUs(f)) ' ' ...
                FAC_labels{f} '\n']);
            end
        end
    end
    
    if ~isempty(ALL_feats)
        EEG_labels = feat_params.EEG_labels;
        EEG_lab_chans = feat_params.EEG_lab_chans;
        FAC_labels = feat_params.FAC_labels;
        FAC_lab_AUs = feat_params.FAC_lab_AUs;
        WAV_labels = feat_params.WAV_labels;
        
        max_eeg_feat = length(EEG_labels);
        max_wav_feat = length(EEG_labels) + length(WAV_labels);
        fprintf(fid, '\nWriting out the EEG, WAV, and FAC features\n');
        
        for i=1:length(ALL_feats)
            f = ALL_feats(i);
            if f <= max_eeg_feat
                fprintf(fid, ['Channel ' num2str(EEG_lab_chans(f)) ' ' ...
                EEG_labels{f} '\n']);
            elseif f <= max_wav_feat
                f = f - max_eeg_feat;
                fprintf(fid, [WAV_labels{f} '\n']);
            else
                f = f - max_wav_feat;
                fprintf(fid, ['AU ' num2str(FAC_lab_AUs(f)) ' ' ...
                FAC_labels{f} '\n']);
            end
        end
    end
    
    fclose(fid);
end