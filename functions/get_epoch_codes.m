function codes = get_epoch_codes(EEG)
%GET_EPOCH_CODES Stimulus code per epoch (event nearest latency 0).

codes = zeros(1, EEG.trials);
for ep = 1:EEG.trials
    lats = EEG.epoch(ep).eventlatency;
    if iscell(lats); lats = cell2mat(lats); end
    [~, k] = min(abs(lats));
    val = EEG.epoch(ep).eventtype;
    if iscell(val); val = val{k}; else; val = val(k); end
    if ischar(val); codes(ep) = str2double(val); else; codes(ep) = val; end
end
end
