path = pwd;

addpath(genpath(pwd));
rmpath([path, '/archive/'])
rmpath([path, '/scenarios/archive/'])
rmpath([path, '/classes/Targets/archive/'])
rmpath([path, '/functions/archive/'])

clear path
clc