%% step5_rsa.m - Representational similarity analysis of face/car EEG patterns.
%
% Tests whether condition-level EEG channel patterns are organised by a
% face-vs-car category model. Neural RDM = 1 - Pearson r between the 80
% condition patterns; model RDM = 0 within category, 1 between. Inference:
% one-sample t-test of subject rho vs 0 with Benjamini-Hochberg FDR.
%
% Caution: shows representational alignment, not causal category coding, and
% does not separate category from low-level visual differences.

project_root = fileparts(fileparts(mfilename('fullpath')));
load(fullfile(project_root, 'pipeline_setup.mat'));
addpath(functions_dir);
load(fullfile(proc_dir, 'event_code_mapping.mat'));
eeglab; close;

load(fullfile(results_dir, 'erp_results.mat'), 'times');

%% Part A: windowed RSA
rsa_windows = {[80 130], [140 200], [200 300], [300 500]};
rsa_labels  = {'P1 (80-130)', 'N170 (140-200)', 'N250 (200-300)', 'Late (300-500)'};
n_windows   = numel(rsa_windows);
rho_windowed = zeros(n_subjects, n_windows);

for s = 1:n_subjects
    fprintf('RSA windowed: subject %s (%d/%d)\n', subjects{s}, s, n_subjects);
    EEG = pop_loadset('filename', sprintf('%s_clean.set', subjects{s}), 'filepath', proc_dir);
    epoch_codes = get_epoch_codes(EEG);

    stim_codes = unique(epoch_codes);
    stim_codes = stim_codes(ismember(stim_codes, [FACE_CODES, CAR_CODES]));
    categories = 1 + ~ismember(stim_codes, FACE_CODES);  % 1 = face, 2 = car
    model_rdm = make_categorical_rdm(categories);

    cond_erps = zeros(numel(stim_codes), EEG.nbchan, EEG.pnts);
    for c = 1:numel(stim_codes)
        cond_erps(c,:,:) = mean(EEG.data(:,:,epoch_codes == stim_codes(c)), 3);
    end

    for w = 1:n_windows
        t_idx = EEG.times >= rsa_windows{w}(1) & EEG.times <= rsa_windows{w}(2);
        patterns = squeeze(mean(cond_erps(:,:,t_idx), 3));
        rho_windowed(s,w) = compare_rdms(compute_neural_rdm(patterns), model_rdm);
    end
end

p_windowed = zeros(1, n_windows);
for w = 1:n_windows
    [~, p_windowed(w), ~, st] = ttest(rho_windowed(:,w), 0);
    fprintf('  %-16s rho=%+.3f  t(%d)=%+.2f  p=%.4f\n', ...
        rsa_labels{w}, mean(rho_windowed(:,w)), st.df, st.tstat, p_windowed(w));
end
[sig_windowed_fdr, pcrit_windowed] = fdr_bh(p_windowed, 0.05);

%% Figure 5: windowed RSA
figure('Position', [100 100 700 450]);
bar_means = mean(rho_windowed, 1);
bar_sems  = std(rho_windowed, 0, 1) / sqrt(n_subjects);
bar(bar_means); hold on;
errorbar(1:n_windows, bar_means, bar_sems, 'k.', 'LineWidth', 1.5);
offsets = linspace(-0.06, 0.06, n_subjects)';  % deterministic, reproducible
for w = 1:n_windows
    scatter(w + offsets, rho_windowed(:,w), 40, 'filled', 'MarkerFaceAlpha', 0.5);
end
yline(0, 'k--');
set(gca, 'XTick', 1:n_windows, 'XTickLabel', rsa_labels);
ylabel('Spearman rho');
title(sprintf('RSA: face/car category model (N=%d)', n_subjects));
for w = 1:n_windows
    if sig_windowed_fdr(w)
        text(w, bar_means(w) + bar_sems(w) + 0.03, '*', 'FontSize', 16, 'HorizontalAlignment', 'center');
    end
end
saveas(gcf, fullfile(results_dir, 'fig5_windowed_RSA.png'));

%% Part B: time-resolved RSA (20 ms window, 10 ms step)
win_ms = 20; step_ms = 10;
t_centres = (times(1) + win_ms/2) : step_ms : (times(end) - win_ms/2);
n_t = numel(t_centres);
rho_ts_all = zeros(n_subjects, n_t);

for s = 1:n_subjects
    fprintf('RSA time-resolved: subject %s (%d/%d)\n', subjects{s}, s, n_subjects);
    EEG = pop_loadset('filename', sprintf('%s_clean.set', subjects{s}), 'filepath', proc_dir);
    epoch_codes = get_epoch_codes(EEG);
    stim_codes = unique(epoch_codes);
    stim_codes = stim_codes(ismember(stim_codes, [FACE_CODES, CAR_CODES]));
    categories = 1 + ~ismember(stim_codes, FACE_CODES);
    model_rdm = make_categorical_rdm(categories);

    cond_erps = zeros(numel(stim_codes), EEG.nbchan, EEG.pnts);
    for c = 1:numel(stim_codes)
        cond_erps(c,:,:) = mean(EEG.data(:,:,epoch_codes == stim_codes(c)), 3);
    end

    for t = 1:n_t
        t_idx = EEG.times >= (t_centres(t) - win_ms/2) & EEG.times <= (t_centres(t) + win_ms/2);
        patterns = squeeze(mean(cond_erps(:,:,t_idx), 3));
        rho_ts_all(s,t) = compare_rdms(compute_neural_rdm(patterns), model_rdm);
    end
end

grand_rho = mean(rho_ts_all, 1);
sem_rho   = std(rho_ts_all, 0, 1) / sqrt(n_subjects);
[~, p_ts] = ttest(rho_ts_all, 0);
[sig_fdr, pcrit_ts] = fdr_bh(p_ts, 0.05);
fprintf('Time-resolved RSA: %d/%d windows significant after FDR q<.05\n', sum(sig_fdr), n_t);

%% Figure 6: time-resolved RSA
figure('Position', [100 100 900 450]);
y_min = min(grand_rho - sem_rho); y_max = max(grand_rho + sem_rho);
y_pad = 0.10 * max(eps, y_max - y_min);
y_limits = [y_min - y_pad, y_max + y_pad];
h_n170 = patch([140 200 200 140], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
    [0.5 0.5 0.5], 'FaceAlpha', 0.08, 'EdgeColor', 'none'); hold on;
h_sem = fill([t_centres, fliplr(t_centres)], [grand_rho + sem_rho, fliplr(grand_rho - sem_rho)], ...
    [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
h_mean = plot(t_centres, grand_rho, 'k-', 'LineWidth', 2);
xline(0, 'k--'); yline(0, 'k--');
handles = [h_sem, h_mean, h_n170]; labels = {'SEM', 'Mean rho', 'N170 window'};
if any(sig_fdr)
    sig_y = y_limits(2) - 0.02 * (y_limits(2) - y_limits(1));
    h_sig = plot(t_centres(sig_fdr), repmat(sig_y, 1, sum(sig_fdr)), 'k.', 'MarkerSize', 10);
    handles(end+1) = h_sig; labels{end+1} = 'FDR q<.05';
end
ylim(y_limits);
xlabel('Time (ms)'); ylabel('Spearman rho (+/- SEM)');
title(sprintf('Time-resolved RSA: face/car category model (N=%d)', n_subjects));
legend(handles, labels, 'Location', 'northwest');
saveas(gcf, fullfile(results_dir, 'fig6_time_resolved_RSA.png'));

%% Part C: group-average RDM at the N170 window
EEG = pop_loadset('filename', sprintf('%s_clean.set', subjects{1}), 'filepath', proc_dir);
epoch_codes = get_epoch_codes(EEG);
stim_codes_all = unique(epoch_codes);
stim_codes_all = stim_codes_all(ismember(stim_codes_all, [FACE_CODES, CAR_CODES]));
n_cond = numel(stim_codes_all);
all_rdms = zeros(n_cond, n_cond, n_subjects);

for s = 1:n_subjects
    EEG = pop_loadset('filename', sprintf('%s_clean.set', subjects{s}), 'filepath', proc_dir);
    epoch_codes = get_epoch_codes(EEG);
    stim_codes = unique(epoch_codes);
    stim_codes = stim_codes(ismember(stim_codes, [FACE_CODES, CAR_CODES]));
    cond_erps = zeros(numel(stim_codes), EEG.nbchan, EEG.pnts);
    for c = 1:numel(stim_codes)
        cond_erps(c,:,:) = mean(EEG.data(:,:,epoch_codes == stim_codes(c)), 3);
    end
    t_idx = EEG.times >= 140 & EEG.times <= 200;
    all_rdms(:,:,s) = compute_neural_rdm(squeeze(mean(cond_erps(:,:,t_idx), 3)));
end

group_rdm = mean(all_rdms, 3);
categories = 1 + ~ismember(stim_codes_all, FACE_CODES);
model_rdm_fig = make_categorical_rdm(categories);
[rho_group, p_group] = compare_rdms(group_rdm, model_rdm_fig);

figure('Position', [100 100 900 400]);
subplot(1,2,1); imagesc(group_rdm); axis square; colorbar;
title(sprintf('Group neural RDM (N=%d, 140-200 ms)', n_subjects));
xlabel('Stimulus condition'); ylabel('Stimulus condition');
subplot(1,2,2); imagesc(model_rdm_fig); axis square; colorbar;
title('Model RDM: face vs car');
xlabel('Stimulus condition'); ylabel('Stimulus condition');
sgtitle(sprintf('Group RDM alignment: Spearman rho = %.3f, p = %.4f', rho_group, p_group));
saveas(gcf, fullfile(results_dir, 'fig7_group_RDMs.png'));

save(fullfile(results_dir, 'rsa_results.mat'), ...
    'rho_windowed', 'rsa_windows', 'rsa_labels', 'p_windowed', ...
    'sig_windowed_fdr', 'pcrit_windowed', ...
    'rho_ts_all', 't_centres', 'grand_rho', 'sem_rho', 'win_ms', 'step_ms', ...
    'sig_fdr', 'p_ts', 'pcrit_ts', 'group_rdm', 'rho_group', 'p_group');

fprintf('\nRSA complete. Next: code/step6_summary.m\n');
