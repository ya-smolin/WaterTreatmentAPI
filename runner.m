s=dbstatus;
save('myBreakpoints.mat', 's');
clear all
clear classes
close all
load('myBreakpoints.mat');
dbstop(s);
format shortG;
controller = Controller();