% startup file
src = pwd;

% add path
restoredefaultpath;
addpath(genpath(src));
rmpath('trash\')
rmpath('trash\others')