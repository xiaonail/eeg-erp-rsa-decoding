%% step0_setup.m - Configure paths, load toolboxes, verify inputs.

clear; clc;

%% Locate repository root and add project functions
project_root = fileparts(fileparts(mfilename('fullpath')));
functions_dir = fullfile(project_root, 'functions');
addpath(functions_dir);

%% Default configuration
raw_dir     = fullfile(project_root, 'raw_data');
proc_dir    = fullfile(project_root, 'processed');
results_dir = fullfile(project_root, 'results');

eeglab_path    = '';
fieldtrip_path = '';
subjects = {'1', '2', '3', '4', '5', '6'};
filename_pattern = '%s_N170.set';

%% Load local configuration if available
config_file = fullfile(project_root, 'code', 'config.m');
if exist(config_file, 'file') == 2
    run(config_file);
else
    fprintf('No config.m found. Using default relative folders.\n');
    fprintf('For local toolbox/data paths, copy code/config_template.m to code/config.m.\n');
end

%% Create output folders
if exist(proc_dir, 'dir') ~= 7
    mkdir(proc_dir);
end
if exist(results_dir, 'dir') ~= 7
    mkdir(results_dir);
end

%% EEGLAB setup
if ~isempty(eeglab_path)
    if exist(eeglab_path, 'dir') ~= 7
        error('eeglab_path does not exist: %s', eeglab_path);
    end
    addpath(eeglab_path);
end

if isempty(which('eeglab'))
    error(['EEGLAB is not on the MATLAB path. Set eeglab_path in config.m, ' ...
           'or add EEGLAB to the MATLAB path before running the pipeline.']);
end

eeglab;
close;

% step2_preprocess.m requires these plugins.
missing_plugins = {};
if isempty(which('iclabel'))
    missing_plugins{end+1} = 'ICLabel'; %#ok<SAGROW>
end
if isempty(which('pop_clean_rawdata')) && isempty(which('clean_rawdata'))
    missing_plugins{end+1} = 'clean_rawdata'; %#ok<SAGROW>
end
if ~isempty(missing_plugins)
    warning(['Missing recommended/required EEGLAB plugins: %s. ' ...
             'Install them through EEGLAB: File > Manage EEGLAB Extensions.'], ...
             strjoin(missing_plugins, ', '));
end

%% FieldTrip setup
if ~isempty(fieldtrip_path)
    if exist(fieldtrip_path, 'dir') ~= 7
        error('fieldtrip_path does not exist: %s', fieldtrip_path);
    end
    addpath(fieldtrip_path);
end

if isempty(which('ft_defaults'))
    error(['FieldTrip is not on the MATLAB path. Set fieldtrip_path in ' ...
           'code/config.m, or add FieldTrip to the MATLAB path before ' ...
           'running the pipeline. Download it from https://www.fieldtriptoolbox.org.']);
end

% ft_defaults adds FieldTrip's subfolders.
ft_defaults;

%% Verify raw subject files
n_subjects = numel(subjects);
all_found = true;

fprintf('\nChecking raw ERP CORE N170 files in:\n  %s\n\n', raw_dir);
for s = 1:n_subjects
    filename = sprintf(filename_pattern, subjects{s});
    file_path = fullfile(raw_dir, filename);

    if exist(file_path, 'file') == 2
        fprintf('  [OK]      %s\n', filename);
    else
        fprintf('  [MISSING] %s\n', filename);
        all_found = false;
    end
end

if ~all_found
    error(['Some raw files are missing. Check raw_dir, subjects, and ' ...
           'filename_pattern in config.m. Raw data are not included in Git.']);
end

%% Save reproducible setup state for downstream scripts
save(fullfile(project_root, 'pipeline_setup.mat'), ...
    'project_root', 'functions_dir', 'raw_dir', 'proc_dir', 'results_dir', ...
    'subjects', 'filename_pattern', 'n_subjects', ...
    'eeglab_path', 'fieldtrip_path');

fprintf('\nSetup complete. Next: code/step1_discover.m\n');
