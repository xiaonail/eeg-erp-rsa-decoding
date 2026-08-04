%% step2_preprocess.m - Clean and epoch raw EEG files (one per subject).
%
% Per subject: add channel locations, drop non-EEG channels, band-pass 0.1-40 Hz,
% clean_rawdata bad-channel removal, average reference, ICA (fit on a 1 Hz copy)
% with ICLabel eye/muscle rejection at p>=0.80, interpolate, epoch, baseline,
% and reject epochs exceeding +/-150 uV.


project_root = fileparts(fileparts(mfilename('fullpath')));
load(fullfile(project_root, 'pipeline_setup.mat'));
addpath(functions_dir);

eeglab;
close;

if isempty(which('iclabel'))
    error('ICLabel is required for this preprocessing script. Install the EEGLAB ICLabel plugin.');
end
if isempty(which('pop_clean_rawdata')) && isempty(which('clean_rawdata'))
    error('clean_rawdata is required for this preprocessing script. Install the EEGLAB clean_rawdata plugin.');
end

%% Event-code mapping
FACE_CODES     = 1:40;
CAR_CODES      = 41:80;
SCRAMBLED_FACE = 101:140;
SCRAMBLED_CAR  = 141:180;
ALL_STIM_CODES = [FACE_CODES, CAR_CODES, SCRAMBLED_FACE, SCRAMBLED_CAR];
ALL_STIM_CODES_STR = arrayfun(@num2str, ALL_STIM_CODES, 'UniformOutput', false);

%% Preprocessing parameters
hp_freq = 0.1;              % high-pass cutoff for ERP-preserving data
lp_freq = 40;               % low-pass cutoff
ica_hp_freq = 1.0;          % high-pass cutoff for ICA training copy
epoch_win = [-0.2, 0.8];    % seconds around stimulus onset
baseline_ms = [-200, 0];    % baseline interval in milliseconds
reject_threshold_uV = 150;  % simple absolute-amplitude epoch rejection
iclabel_threshold = 0.80;   % reject eye/muscle ICs at probability >= threshold

non_eeg_tokens = {'EOG', 'HEOG', 'VEOG', 'EXG', 'STATUS', 'STI', 'TRIGGER'};

%% Open preprocessing log
log_file = fullfile(results_dir, 'preprocessing_log.txt');
fid = fopen(log_file, 'w');
if fid == -1
    error('Cannot open preprocessing log for writing: %s', log_file);
end
cleanup_log = onCleanup(@() fclose(fid));

fprintf(fid, 'ERP CORE N170 preprocessing log\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));
fprintf(fid, 'hp=%.2f Hz, lp=%.2f Hz, ICA hp=%.2f Hz, epoch=[%.1f %.1f] s, baseline=[%d %d] ms\n\n', ...
    hp_freq, lp_freq, ica_hp_freq, epoch_win(1), epoch_win(2), baseline_ms(1), baseline_ms(2));

%% Subject loop
for s = 1:n_subjects
    subject_id = subjects{s};
    filename = sprintf(filename_pattern, subject_id);
    fprintf('\n=== Subject %s (%d/%d) ===\n', subject_id, s, n_subjects);
    fprintf(fid, 'Subject %s\n', subject_id);

    EEG = pop_loadset('filename', filename, 'filepath', raw_dir);
    EEG = eeg_checkset(EEG);
    EEG.setname = subject_id;

    % Standard channel positions.
    lookup_file = fullfile(fileparts(which('eeglab')), ...
        'plugins', 'dipfit', 'standard_BEM', 'elec', 'standard_1005.elc');
    if exist(lookup_file, 'file') == 2
        EEG = pop_chanedit(EEG, 'lookup', lookup_file);
    else
        warning('Standard channel-location file not found; using existing locations.');
    end

    % Remove non-EEG/status channels.
    remove_channels = {};
    for ch = 1:EEG.nbchan
        label_upper = upper(EEG.chanlocs(ch).labels);
        if any(contains(label_upper, non_eeg_tokens))
            remove_channels{end+1} = EEG.chanlocs(ch).labels; %#ok<SAGROW>
        end
    end
    if ~isempty(remove_channels)
        EEG = pop_select(EEG, 'nochannel', remove_channels);
    end

    original_chanlocs = EEG.chanlocs;

    % Band-pass filter.
    EEG = pop_eegfiltnew(EEG, 'locutoff', hp_freq, 'hicutoff', lp_freq);

    % Bad-channel detection; burst/window cleaning disabled.
    EEG = pop_clean_rawdata(EEG, ...
        'FlatlineCriterion', 5, ...
        'ChannelCriterion', 0.8, ...
        'LineNoiseCriterion', 4, ...
        'Highpass', 'off', ...
        'BurstCriterion', 'off', ...
        'WindowCriterion', 'off');

    removed_bad_channels = setdiff({original_chanlocs.labels}, {EEG.chanlocs.labels});
    if isempty(removed_bad_channels)
        removed_bad_channels = {'none'};
    end
    fprintf('  Removed channels: %s\n', strjoin(removed_bad_channels, ', '));
    fprintf(fid, '  Removed channels: %s\n', strjoin(removed_bad_channels, ', '));

    EEG = pop_reref(EEG, []);

    % ICA on a 1 Hz copy; transfer weights back to the 0.1 Hz data.
    EEG_ica = pop_eegfiltnew(EEG, 'locutoff', ica_hp_freq);
    EEG_ica = pop_runica(EEG_ica, 'icatype', 'runica', 'extended', 1);
    EEG.icaweights = EEG_ica.icaweights;
    EEG.icasphere  = EEG_ica.icasphere;
    EEG = eeg_checkset(EEG, 'ica');
    clear EEG_ica;

    EEG = iclabel(EEG);
    % ICLabel class order: Brain, Muscle, Eye, Heart, Line Noise, Channel Noise, Other.
    EEG = pop_icflag(EEG, [NaN NaN; ...
                           iclabel_threshold 1; ...
                           iclabel_threshold 1; ...
                           NaN NaN; NaN NaN; NaN NaN; NaN NaN]);
    rejected_ics = find(EEG.reject.gcompreject);
    fprintf('  Rejected ICs: %d / %d\n', numel(rejected_ics), size(EEG.icaweights, 1));
    fprintf(fid, '  Rejected ICs: %d / %d\n', numel(rejected_ics), size(EEG.icaweights, 1));
    EEG = pop_subcomp(EEG, rejected_ics, 0);

    EEG = pop_interp(EEG, original_chanlocs, 'spherical');
    EEG = pop_reref(EEG, []);

    % Match event-type format: char or numeric, whichever this file uses.
    if isnumeric(EEG.event(1).type)
        epoch_selector = num2cell(ALL_STIM_CODES);
        n_matching = sum(ismember([EEG.event.type], ALL_STIM_CODES));
    else
        epoch_selector = ALL_STIM_CODES_STR;
        n_matching = sum(ismember({EEG.event.type}, ALL_STIM_CODES_STR));
    end

    if n_matching == 0
        error(['No events matched ALL_STIM_CODES for subject %s. Run ' ...
               'step1_discover.m and correct the code mapping at the top ' ...
               'of this script.'], subject_id);
    end

    EEG = pop_epoch(EEG, epoch_selector, epoch_win);
    n_epochs_before_rejection = EEG.trials;
    EEG = pop_rmbase(EEG, baseline_ms);
    EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan, ...
        -reject_threshold_uV, reject_threshold_uV, EEG.xmin, EEG.xmax, 0, 1);

    retention_percent = EEG.trials / max(n_epochs_before_rejection, 1) * 100;
    fprintf('  Epochs retained: %d / %d (%.1f%%)\n', EEG.trials, n_epochs_before_rejection, retention_percent);
    fprintf(fid, '  Epochs retained: %d / %d (%.1f%%)\n\n', EEG.trials, n_epochs_before_rejection, retention_percent);

    EEG.setname = sprintf('%s_clean', subject_id);
    clean_name = sprintf('%s_clean.set', subject_id);
    pop_saveset(EEG, 'filename', clean_name, 'filepath', proc_dir);

    % Rewrite as v7 so scipy/MNE can read the .set header.
    clean_path = fullfile(proc_dir, clean_name);
    S = load('-mat', clean_path);
    save(clean_path, '-struct', 'S', '-v7');
end

save(fullfile(proc_dir, 'event_code_mapping.mat'), ...
    'FACE_CODES', 'CAR_CODES', 'SCRAMBLED_FACE', 'SCRAMBLED_CAR', 'ALL_STIM_CODES');

fprintf('\nPreprocessing complete. Next: code/step3_erp.m\n');
