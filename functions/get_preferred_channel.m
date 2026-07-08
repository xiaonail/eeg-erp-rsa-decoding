function [channel_index, channel_label] = get_preferred_channel(labels, preferred_label, fallback_labels)
%GET_PREFERRED_CHANNEL Return a preferred channel index with ordered fallbacks.
%
%   [IDX, LABEL] = GET_PREFERRED_CHANNEL(LABELS, PREFERRED_LABEL, FALLBACK_LABELS)
%   searches LABELS first for PREFERRED_LABEL, then for each fallback label.
%   If none are present, the first channel is returned with a warning.
%
%   This keeps analysis scripts explicit about region-of-interest choices while
%   still allowing tutorial data with slightly different channel montages to run.

if nargin < 3 || isempty(fallback_labels)
    fallback_labels = {};
end

candidate_labels = [{preferred_label}, fallback_labels(:)'];
channel_index = [];
channel_label = '';

for i = 1:numel(candidate_labels)
    this_label = candidate_labels{i};
    channel_index = find(strcmp(labels, this_label), 1, 'first');
    if ~isempty(channel_index)
        channel_label = this_label;
        return;
    end
end

channel_index = 1;
channel_label = labels{1};
warning('None of the requested channels were found. Using first channel: %s.', channel_label);
end
