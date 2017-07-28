%% System

%% Path

% irt
run([pwd, '/irt/setup.m']);

% ute
addpath(genpath([pwd, '/ute']))

% util
addpath([pwd])


%% Inline Function

chk_not_win = @(x) ~strcmp(computer, 'PCWIN') && ~strcmp(computer, 'PCWIN64');

