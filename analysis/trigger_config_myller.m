% Müller 2003

config.position = {
    'start',      '0',  2;
    'other',      '1',  2;};

config.cond = {
    '1+2',       '1100',  3:6;
    '1+3',       '1010',  3:6;
    '2+4',       '0101',  3:6
    '3+4',       '0011',  3:6;};


config.hand = {
    'left',     '0', 7;
    'right',    '1', 7;};


config.resp = {
    'resp',     '1',     8;
    'stim', 	'0',     8;};

% Similarly, add other categories...

save('trigger_config_myller.mat', 'config');
