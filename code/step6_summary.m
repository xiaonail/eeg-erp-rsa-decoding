%% step6_summary.m - Summary of ERP, RSA, and TFR results.

project_root = fileparts(fileparts(mfilename('fullpath')));
load(fullfile(project_root, 'pipeline_setup.mat'));
addpath(functions_dir);

erp = load(fullfile(results_dir, 'erp_results.mat'));
rsa = load(fullfile(results_dir, 'rsa_results.mat'));
tfr_file = fullfile(results_dir, 'tfr_results.mat');
has_tfr = exist(tfr_file, 'file') == 2;
if has_tfr
    tfr = load(tfr_file, 'freqs', 'tfr_toi', 'plot_chan', 'channel_used');
end

fprintf('\n===== ERP CORE N170 tutorial summary (N=%d) =====\n', n_subjects);

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
if has_tfr
    fprintf('  Figure available: fig4_group_TFR.png\n');
    fprintf('  Plot channel: %s\n', tfr.plot_chan);
    fprintf('  Frequency range: %.0f-%.0f Hz\n', min(tfr.freqs), max(tfr.freqs));
else
    fprintf('  Not available. FieldTrip step was skipped or has not been run.\n');
end

fprintf('\nGenerated outputs are in:\n  %s\n', results_dir);
fprintf('\nRemember: these outputs are tutorial demonstrations unless you run the full dataset and report all analysis choices.\n');
