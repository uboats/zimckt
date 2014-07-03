function ckthelp()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% ckthelp: print available instructions
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n');
%disp('  Available Instructions:');
disp('  - help: display all available instructions');
disp('  - ls: list files in current folder');
disp('  - run: start an analysis');
disp('     * 1st: circuit file is required');
disp('            including the path');
disp('     * 2nd: if doing transient analysis,');
disp('            the TRAN type should be claimed.');
disp('            Otherwise, just type ''Enter'' for default');
disp('  - clc: clean the screen');
disp('  - clear: clear all the memory');
disp('  - quit: quit the program');
fprintf('\n');

end

