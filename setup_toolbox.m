%% System

%% Path

% irt
run([pwd, '/irt/setup.m']);

% bwdistsc
addpath(genpath([pwd, '/bwdistsc']))

% ute
addpath(genpath([pwd, '/ute']))

% ssqsm
addpath(genpath([pwd, '/ssqsm_bilgic']))

% util
addpath([pwd])



%% Inline Function

chk_not_win = @(x) ~strcmp(computer, 'PCWIN') && ~strcmp(computer, 'PCWIN64');

