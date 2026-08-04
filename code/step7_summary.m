%% step7_summary.m - Summary of ERP, RSA, TFR, and decoding results.
%
% Reads the .mat files from steps 3-6 and prints one report. 

project_root = fileparts(fileparts(mfilename('fullpath')));
load(fullfile(project_root, 'pipeline_setup.mat'));
addpath(functions_dir);

erp = load(fullfile(results_dir, 'erp_results.mat'));
rsa = load(fullfile(results_dir, 'rsa_results.mat'));
tfr = load(fullfile(results_dir, 'tfr_results.mat'), ...
    'freqs', 'tfr_toi', 'plot_chan', 'channel_used');

fprintf('\n===== ERP CORE N170 pipeline summary (N=%d) =====\n', n_subjects);

%% Trial counts
fprintf('\nAccepted trial counts after preprocessing:\n');
fprintf('  Subject | Face | Car\n');
fprintf('  --------|------|----\n');
for s = 1:n_subjects
    fprintf('  %-7s | %4d | %3d\n', subjects{s}, erp.trial_counts(s,1), erp.trial_counts(s,2));
end

%% ERP summary
fprintf('\nERP N170 contrast, %d-%d ms, posterior N170 electrodes:\n', ...
    erp.n170_window_ms(1), erp.n170_window_ms(2));
fprintf('  Face mean: %.2f %s\n', mean(erp.n170_amp(:,1)) * erp.erp_scale, erp.erp_units);
fprintf('  Car mean:  %.2f %s\n', mean(erp.n170_amp(:,2)) * erp.erp_scale, erp.erp_units);
fprintf('  t(%d)=%.2f, p=%.4f, dz=%.2f\n', ...
    erp.stats_n170.df, erp.stats_n170.tstat, erp.p_n170, erp.cohens_dz);

%% RSA summary
fprintf('\nWindowed RSA against face/car category model:\n');
for w = 1:numel(rsa.rsa_windows)
    if rsa.sig_windowed_fdr(w)
        sig_label = 'FDR significant';
    else
        sig_label = '';
    end
    fprintf('  %-16s mean rho=%+.3f, p=%.4f  %s\n', ...
        rsa.rsa_labels{w}, mean(rsa.rho_windowed(:,w)), rsa.p_windowed(w), sig_label);
end

fprintf('\nTime-resolved RSA:\n');
fprintf('  FDR-significant windows: %d / %d\n', sum(rsa.sig_fdr), numel(rsa.sig_fdr));
if any(rsa.sig_fdr)
    fprintf('  Significant time range: %.0f to %.0f ms\n', ...
        min(rsa.t_centres(rsa.sig_fdr)), max(rsa.t_centres(rsa.sig_fdr)));
end
[peak_rho, peak_idx] = max(rsa.grand_rho);
fprintf('  Peak mean rho: %.3f at %.0f ms\n', peak_rho, rsa.t_centres(peak_idx));
fprintf('  Group N170-window RDM alignment: rho=%.3f, p=%.4f\n', rsa.rho_group, rsa.p_group);

%% TFR summary
fprintf('\nTime-frequency output:\n');
fprintf('  Figure available: fig4_group_TFR.png\n');
fprintf('  Plot channel: %s\n', tfr.plot_chan);
fprintf('  Frequency range: %.0f-%.0f Hz\n', min(tfr.freqs), max(tfr.freqs));

%% Decoding summary
dec_file = fullfile(results_dir, 'decoding_results.mat');
fprintf('\nMultivariate decoding, face vs car:\n');
if exist(dec_file, 'file') == 2
    dec = load(dec_file);

    % scipy saves logicals as uint8; cast before using as a mask.
    sig_dec = logical(dec.sig);

    fprintf('  Subjects analysed: %d (%d-fold CV, C=%.2f, %d Hz)\n', ...
        numel(dec.subjects), dec.n_folds, dec.c_param, dec.sfreq_decode);
    fprintf('  Peak AUC: %.3f at %.0f ms\n', dec.peak_auc, dec.peak_ms);
    n_clust = sum(dec.cluster_pv < 0.05);
    fprintf('  Significant clusters: %d (%d / %d time points)\n', ...
        n_clust, sum(sig_dec), numel(sig_dec));
    if any(sig_dec)
        fprintf('  Significant time range: %.0f to %.0f ms\n', ...
            min(dec.times(sig_dec)), max(dec.times(sig_dec)));
        fprintf('  Smallest cluster p: %.4f\n', min(dec.cluster_pv));
        fprintf('  Decoding onset: %.0f ms\n', dec.onset);
    else
        fprintf('  No cluster above chance after permutation testing.\n');
    end

    [~, p_dec, ~, st_dec] = ttest(dec.n170, 0.5);
    fprintf('  N170 window: mean AUC=%.3f, t(%d)=%.2f, p=%.4f\n', ...
        mean(dec.n170), st_dec.df, st_dec.tstat, p_dec);
    fprintf('  Figures: fig8_time_resolved_decoding.png, fig9_temporal_generalisation.png\n');
else
    fprintf('  Not available. Run code/step6_decoding.py first.\n');
end

fprintf('\nGenerated outputs are in:\n  %s\n', results_dir);
fprintf('\nDone.\n');
