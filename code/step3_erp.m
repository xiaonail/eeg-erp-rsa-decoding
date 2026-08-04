%% step3_erp.m - ERP grand averages and N170 face-vs-car statistics.
%
% Subject and grand-average ERPs for faces and cars; paired t-test on mean N170
% amplitude (140-200 ms, P7/P8/PO7/PO8). 

project_root = fileparts(fileparts(mfilename('fullpath')));
load(fullfile(project_root, 'pipeline_setup.mat'));
addpath(functions_dir);

load(fullfile(proc_dir, 'event_code_mapping.mat'));

eeglab;
close;

%% Load first cleaned subject to define common dimensions
EEG_template = pop_loadset('filename', sprintf('%s_clean.set', subjects{1}), 'filepath', proc_dir);

n_channels = EEG_template.nbchan;
n_times = EEG_template.pnts;
times = EEG_template.times;
chanlocs = EEG_template.chanlocs;
clear EEG_template;

%% N170 region/time-window definition
n170_channel_names = {'P7', 'P8', 'PO7', 'PO8'};
all_channel_names = {chanlocs.labels};
n170_channel_idx = find(ismember(all_channel_names, n170_channel_names));
if isempty(n170_channel_idx)
    error('No predefined N170 channels found. Expected one or more of: %s.', strjoin(n170_channel_names, ', '));
end

n170_window_ms = [140, 200];
n170_time_idx = times >= n170_window_ms(1) & times <= n170_window_ms(2);

%% Preallocate outputs
erps_face = NaN(n_subjects, n_channels, n_times);
erps_car = NaN(n_subjects, n_channels, n_times);
n170_amp = NaN(n_subjects, 2);       % columns: face, car
trial_counts = zeros(n_subjects, 2); % columns: face, car

%% Subject-level ERP extraction
for s = 1:n_subjects
    subject_id = subjects{s};
    fprintf('ERP: subject %s (%d/%d)\n', subject_id, s, n_subjects);

    clean_filename = sprintf('%s_clean.set', subject_id);
    EEG = pop_loadset('filename', clean_filename, 'filepath', proc_dir);

    % Dimension check against the template.
    if EEG.nbchan ~= n_channels || EEG.pnts ~= n_times
        error(['Subject %s has %d channels x %d samples, expected %d x %d. ' ...
               'Check the interpolation step in step2_preprocess.m.'], ...
               subject_id, EEG.nbchan, EEG.pnts, n_channels, n_times);
    end

    epoch_codes = get_epoch_codes(EEG);
    face_trials = ismember(epoch_codes, FACE_CODES);
    car_trials = ismember(epoch_codes, CAR_CODES);
    trial_counts(s,:) = [sum(face_trials), sum(car_trials)];
    if trial_counts(s,1) == 0 || trial_counts(s,2) == 0
        error(['Subject %s has %d face and %d car trials after rejection. ' ...
               'A condition with no trials produces NaN averages.'], ...
               subject_id, trial_counts(s,1), trial_counts(s,2));
    end

    erps_face(s,:,:) = mean(EEG.data(:,:,face_trials), 3);
    erps_car(s,:,:) = mean(EEG.data(:,:,car_trials), 3);

    n170_amp(s,1) = mean(EEG.data(n170_channel_idx, n170_time_idx, face_trials), 'all');
    n170_amp(s,2) = mean(EEG.data(n170_channel_idx, n170_time_idx, car_trials), 'all');
end

%% Grand averages
grand_face = squeeze(mean(erps_face, 1));
grand_car = squeeze(mean(erps_car, 1));
grand_diff = grand_face - grand_car;

% Auto-detect volts vs microvolts for display.
if max(abs(grand_face(:))) < 1e-3
    erp_scale = 1e6;
else
    erp_scale = 1;
end
erp_units = 'uV';

%% Figure 1: grand-average waveform at N170 electrodes
face_wave = mean(grand_face(n170_channel_idx,:), 1) * erp_scale;
car_wave = mean(grand_car(n170_channel_idx,:), 1) * erp_scale;
y_min = min([face_wave, car_wave]);
y_max = max([face_wave, car_wave]);
y_pad = max(1, 0.10 * (y_max - y_min));

figure('Position', [100 100 800 500]);
h_patch = patch([n170_window_ms(1) n170_window_ms(2) n170_window_ms(2) n170_window_ms(1)], ...
    [y_min-y_pad y_min-y_pad y_max+y_pad y_max+y_pad], ...
    [0.5 0.5 0.5], 'FaceAlpha', 0.10, 'EdgeColor', 'none');
hold on;
plot(times, face_wave, '-', 'LineWidth', 2);
plot(times, car_wave, '--', 'LineWidth', 2);
xline(0, 'k--');
yline(0, 'k--');
set(gca, 'YDir', 'reverse');
xlabel('Time (ms)');
ylabel(sprintf('Amplitude (%s)', erp_units));
legend('N170 window', 'Faces', 'Cars', 'Location', 'northwest');
title(sprintf('Grand-average ERP at posterior N170 electrodes (N=%d)', n_subjects));
saveas(gcf, fullfile(results_dir, 'fig1_grand_average_N170.png'));

%% Figure 2: face-minus-car difference topographies
figure('Position', [100 100 1000 300]);
topography_times_ms = [100, 170, 250, 400];
for i = 1:numel(topography_times_ms)
    [~, time_idx] = min(abs(times - topography_times_ms(i)));
    subplot(1, numel(topography_times_ms), i);
    topoplot(grand_diff(:, time_idx) * erp_scale, chanlocs, ...
        'maplimits', 'maxmin', 'electrodes', 'on');
    title(sprintf('%d ms', topography_times_ms(i)));
    colorbar;
end
sgtitle(sprintf('Face minus car topography (N=%d)', n_subjects));
saveas(gcf, fullfile(results_dir, 'fig2_topographies.png'));

%% Figure 3: subject-level N170 amplitudes
figure('Position', [100 100 600 400]);
bar(n170_amp * erp_scale);
xlabel('Subject');
ylabel(sprintf('Mean amplitude (%s)', erp_units));
legend('Face', 'Car', 'Location', 'best');
title('N170 amplitude by subject (140-200 ms)');
saveas(gcf, fullfile(results_dir, 'fig3_individual_N170.png'));

%% Group-level paired test
diff_face_minus_car = n170_amp(:,1) - n170_amp(:,2);
[~, p_n170, ci_n170, stats_n170] = ttest(n170_amp(:,1), n170_amp(:,2));
cohens_dz = mean(diff_face_minus_car) / std(diff_face_minus_car, 0);

fprintf('\nN170 face-car contrast (N=%d):\n', n_subjects);
fprintf('  Face mean: %.2f %s\n', mean(n170_amp(:,1)) * erp_scale, erp_units);
fprintf('  Car mean:  %.2f %s\n', mean(n170_amp(:,2)) * erp_scale, erp_units);
fprintf('  t(%d)=%.2f, p=%.4f, dz=%.2f\n', stats_n170.df, stats_n170.tstat, p_n170, cohens_dz);

save(fullfile(results_dir, 'erp_results.mat'), ...
    'erps_face', 'erps_car', 'grand_face', 'grand_car', 'grand_diff', ...
    'n170_amp', 'trial_counts', 'times', 'chanlocs', 'n170_channel_idx', ...
    'n170_channel_names', 'n170_window_ms', 'p_n170', 'ci_n170', ...
    'stats_n170', 'cohens_dz', 'erp_scale', 'erp_units');

fprintf('\nERP analysis complete. Next: code/step4_timefreq.m or code/step5_rsa.m\n');
