%% step4_timefreq.m - Morlet time-frequency analysis (4-40 Hz), face vs car.
%
% Requires FieldTrip. Per subject: wavelet power (5 cycles), dB baseline
% correction (-200..0 ms), extracted at PO8 (fallback P8).

project_root = fileparts(fileparts(mfilename('fullpath')));
load(fullfile(project_root, 'pipeline_setup.mat'));
addpath(functions_dir);

load(fullfile(proc_dir, 'event_code_mapping.mat'));

eeglab;
close;
if isempty(which('ft_defaults'))
    error('FieldTrip is not on the MATLAB path. Run step0_setup.m first.');
end
ft_defaults;

%% TFR parameters
freqs = 4:1:40;
n_cycles = 5;
tfr_toi = -0.2:0.01:0.8;
baseline_window_sec = [-0.2 0];
preferred_channel = 'PO8';
fallback_channels = {'P8', 'PO7', 'P7'};
plot_chan = '';

n_freqs = numel(freqs);
n_times = numel(tfr_toi);
tfr_face_all = NaN(n_subjects, n_freqs, n_times);
tfr_car_all = NaN(n_subjects, n_freqs, n_times);
channel_used = cell(n_subjects, 1);

%% Subject-level FieldTrip TFR
for s = 1:n_subjects
    subject_id = subjects{s};
    fprintf('TFR: subject %s (%d/%d)\n', subject_id, s, n_subjects);

    clean_filename = sprintf('%s_clean.set', subject_id);
    EEG = pop_loadset('filename', clean_filename, 'filepath', proc_dir);
    epoch_codes = get_epoch_codes(EEG);

    face_trials = find(ismember(epoch_codes, FACE_CODES));
    car_trials = find(ismember(epoch_codes, CAR_CODES));
    if isempty(face_trials) || isempty(car_trials)
        warning('Subject %s has missing face or car trials. Skipping TFR for this subject.', subject_id);
        continue;
    end

    data_ft = eeglab2fieldtrip(EEG, 'preprocessing');

    cfg = [];
    cfg.trials = face_trials;
    data_face = ft_selectdata(cfg, data_ft);

    cfg = [];
    cfg.trials = car_trials;
    data_car = ft_selectdata(cfg, data_ft);

    cfg = [];
    cfg.method = 'wavelet';
    cfg.output = 'pow';
    cfg.foi = freqs;
    cfg.width = n_cycles;
    cfg.toi = tfr_toi;
    cfg.pad = 'nextpow2';
    tfr_face = ft_freqanalysis(cfg, data_face);
    tfr_car = ft_freqanalysis(cfg, data_car);

    cfg = [];
    cfg.baseline = baseline_window_sec;
    cfg.baselinetype = 'db';
    tfr_face = ft_freqbaseline(cfg, tfr_face);
    tfr_car = ft_freqbaseline(cfg, tfr_car);

    [channel_idx, this_channel] = get_preferred_channel(tfr_face.label, preferred_channel, fallback_channels);
    channel_used{s} = this_channel;
    if isempty(plot_chan)
        plot_chan = this_channel;
    elseif ~strcmp(plot_chan, this_channel)
        warning('Different TFR channels were used across subjects. Check channel_used in tfr_results.mat.');
    end

    tfr_face_all(s,:,:) = squeeze(tfr_face.powspctrm(channel_idx,:,:));
    tfr_car_all(s,:,:) = squeeze(tfr_car.powspctrm(channel_idx,:,:));
end

%% Group averages
grand_tfr_face = squeeze(mean(tfr_face_all, 1));
grand_tfr_car = squeeze(mean(tfr_car_all, 1));
grand_tfr_diff = grand_tfr_face - grand_tfr_car;

%% Figure 4: descriptive TFR panels
figure('Position', [50 100 1400 400]);
subplot(1,3,1);
imagesc(tfr_toi*1000, freqs, grand_tfr_face);
axis xy;
colorbar;
caxis([-3 3]);
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title(sprintf('Faces (dB) - %s', plot_chan));

subplot(1,3,2);
imagesc(tfr_toi*1000, freqs, grand_tfr_car);
axis xy;
colorbar;
caxis([-3 3]);
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title(sprintf('Cars (dB) - %s', plot_chan));

subplot(1,3,3);
imagesc(tfr_toi*1000, freqs, grand_tfr_diff);
axis xy;
colorbar;
caxis([-2 2]);
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title(sprintf('Face minus car - %s', plot_chan));

sgtitle(sprintf('Group time-frequency power (N=%d)', n_subjects));
saveas(gcf, fullfile(results_dir, 'fig4_group_TFR.png'));

save(fullfile(results_dir, 'tfr_results.mat'), ...
    'tfr_face_all', 'tfr_car_all', 'grand_tfr_face', 'grand_tfr_car', ...
    'grand_tfr_diff', 'freqs', 'tfr_toi', 'baseline_window_sec', ...
    'preferred_channel', 'fallback_channels', 'plot_chan', 'channel_used');

fprintf('\nTFR analysis complete. Next: code/step5_rsa.m\n');
