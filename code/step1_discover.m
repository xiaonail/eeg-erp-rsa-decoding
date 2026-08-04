%% step1_discover.m - Inspect one subject's event codes, channels, sampling rate.

project_root = fileparts(fileparts(mfilename('fullpath')));
load(fullfile(project_root, 'pipeline_setup.mat'));
addpath(functions_dir);

eeglab;
close;

%% Load first subject as a metadata example
filename = sprintf(filename_pattern, subjects{1});
EEG = pop_loadset('filename', filename, 'filepath', raw_dir);
EEG = eeg_checkset(EEG);

%% Basic recording structure
fprintf('\n=== Example file: %s ===\n', filename);
fprintf('Sampling rate: %g Hz\n', EEG.srate);
fprintf('Channels:      %d\n', EEG.nbchan);
fprintf('Samples:       %d\n', EEG.pnts);
fprintf('Duration:      %.1f s\n', EEG.pnts / EEG.srate);

fprintf('\n=== Channel labels ===\n');
channel_labels = {EEG.chanlocs.labels};
for ch = 1:numel(channel_labels)
    fprintf('  %3d: %-10s', ch, channel_labels{ch});
    if mod(ch, 5) == 0
        fprintf('\n');
    end
end
fprintf('\n');

%% Event-code audit
fprintf('\n=== Event-code counts ===\n');
if isempty(EEG.event)
    warning('No EEG.event entries found. Check whether the raw file is valid.');
else
    event_types = {EEG.event.type};

    % Normalise to strings for counting across numeric and char formats.
    event_labels = cell(size(event_types));
    for i = 1:numel(event_types)
        if isnumeric(event_types{i})
            event_labels{i} = num2str(event_types{i});
        elseif isstring(event_types{i}) || ischar(event_types{i})
            event_labels{i} = char(event_types{i});
        else
            event_labels{i} = '<unreadable>';
        end
    end

    unique_labels = unique(event_labels);
    fprintf('  Code         | Count\n');
    fprintf('  -------------|------\n');
    for c = 1:numel(unique_labels)
        fprintf('  %-12s | %5d\n', unique_labels{c}, sum(strcmp(event_labels, unique_labels{c})));
    end
end

fprintf('\nDiscovery complete. Next: code/step2_preprocess.m\n');
