% Müller 2003

config.position = {
    'start',      '11',  6:7;
    'calib',      '10',  6:7;
    'resp',       '00',  6:7;};

config.training = {
    'training',     '1',     5;
    'experiment', 	'0',     5;};

config.cond = {
    '1+2',       '1100',  8:11;
    '1+3',       '1010',  8:11;
    '2+4',       '0101',  8:11;
    '3+4',       '0011',  8:11;
    '1',         '0001',  8:11;
    '2',         '0010',  8:11;
    '3',         '0011',  8:11;
    '4',         '0101',  8:11;
    '5',         '0110',  8:11;
    '6',         '0111',  8:11;};



% config.resp = {
%     'resp',     '1',     11;
%     'stim', 	'0',     11;};
% 
% config.position = {
%     'start',      '0',  5;
%     'other',      '1',  5;};
% 
% config.cond = {
%     '1+2',       '1100',  6:9;
%     '1+3',       '1010',  6:9;
%     '2+4',       '0101',  6:9
%     '3+4',       '0011',  6:9;};
% 
% 
% config.hand = {
%     'left',     '0', 10;
%     'right',    '1', 10;};


% Similarly, add other categories...

save('trigger_config_myller.mat', 'config');
