%% step0_setup.m - Configure paths, load toolboxes, verify inputs. Run first.
%
% Machine-specific paths live in config.m (git-ignored). On a new machine,
% copy config_template.m to config.m and edit toolbox/data paths.

clear; clc;

%% Locate repository root and add project functions
project_root = fileparts(fileparts(mfilename('fullpath')));
functions_dir = fullfile(project_root, 'functions');
addpath(functions_dir);

%% Default configuration
% These defaults make the repository self-contained. Users can override any of
% them in config.m without editing public analysis scripts.
raw_dir     = fullfile(project_root, 'raw_data');
proc_dir    = fullfile(project_root, 'processed');
results_dir = fullfile(project_root, 'results');

eeglab_path    = '';
fieldtrip_path = '';
subjects = {'1', '2', '3', '4', '5', '6'};
filename_pattern = '%s_N170.set';

%% Load local configuration if available
config_file = fullfile(project_root, 'config.m');
if exist(config_file, 'file') == 2
    run(config_file);
else
    fprintf('No config.m found. Using default relative folders.\n');
    fprintf('For local toolbox/data paths, copy config_template.m to config.m.\n');
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

% step2_preprocess.m requires these plugins. We warn here so users can fix the
% environment before a long preprocessing run starts.
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

%% FieldTrip setup, optional for time-frequency analysis only
fieldtrip_available = false;
if ~isempty(fieldtrip_path)
    if exist(fieldtrip_path, 'dir') == 7
        addpath(fieldtrip_path);
        ft_defaults;
        fieldtrip_available = true;
    else
        warning('fieldtrip_path does not exist: %s', fieldtrip_path);
    end
elseif ~isempty(which('ft_defaults'))
    ft_defaults;
    fieldtrip_available = true;
end

if ~fieldtrip_available
    warning('FieldTrip not found. step4_timefreq.m will be skipped by run_all.m.');
end

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
    'eeglab_path', 'fieldtrip_path', 'fieldtrip_available');

fprintf('\nSetup complete. Next: code/step1_discover.m\n');
