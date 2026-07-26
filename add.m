desktop = com.mathworks.mde.desk.MLDesktop.getInstance();
desktop.closeGroup('Editor');
restoredefaultpath
addpath(genpath(pwd))
rmpath(genpath('archive/'))
rmpath(genpath('functions/archive/'))
rmpath(genpath('live_scripts/archive/'))
rmpath(genpath('scenarios/archive/'))
clear
clc