% Müller 2003
% trigN = '1'+ str(isTraining) + trigdic['calib'] + trigdic[str(current_point)]
% trigN = '1'+ str(isTraining) + trigdic['start'] + trigdic[attendCond]

config.training = {
    'training',     '1', 2;
    'experiment', 	'0', 2;};

config.position = {
    'start',      '11',  3:4;
    'calib',      '10',  3:4;
    'resp',       '00',  3:4;};

config.cond = {
    '1+2',       '1100',  5:8;
    '1+3',       '1010',  5:8;
    '2+4',       '0101',  5:8;
    '3+4',       '0011',  5:8;
    '1',         '0001',  5:8;
    '2',         '0010',  5:8;
    '3',         '0011',  5:8;
    '4',         '0101',  5:8;
    '5',         '0110',  5:8;
    '6',         '0111',  5:8;};


% config.resp = {
%     'resp',     '1',     8;
%     'stim', 	'0',     8;};
% 
% config.position = {
%     'start',      '0',  2;
%     'other',      '1',  2;};
% 
% config.cond = {
%     '1+2',       '1100',  3:6;
%     '1+3',       '1010',  3:6;
%     '2+4',       '0101',  3:6
%     '3+4',       '0011',  3:6;};
% 
% 
% config.hand = {
%     'left',     '0', 7;
%     'right',    '1', 7;};

% Similarly, add other categories...

save('trigger_config_myller.mat', 'config');
