%————————————————————————————————————————————
%  ██████╗███████╗████████╗██╗   ██╗██████╗ 
% ██╔════╝██╔════╝╚══██╔══╝██║   ██║██╔══██╗
% ╚█████╗ ███████╗   ██║   ██║   ██║██████╔╝
%  ╚═══██╗██╔════╝   ██║   ██║   ██║██╔═══╝ 
% ██████╔╝███████╗   ██║   ╚██████╔╝██║     
% ╚═════╝ ╚══════╝   ╚═╝    ╚═════╝ ╚═╝     
%————————————————————————————————————————————

% Load the correct python environment for the neural network training and
% prediction
try
    pyenv('Version',"C:\Users\maxja\AppData\Local\Programs\Python\Python310\python.exe");
catch
    error('Please go into the all_sim_setup.m file to change the python executable location, Neural Network will only give errors')
end

% Set some standard bode options, enabling phase wrapping and setting the
% frequency unit to hertz
opt = bodeoptions('cstprefs');
opt.PhaseWrapping = 'on';
opt.FreqUnits = 'Hz';

lw = 1.5; % Define the wanted linewidth for the plots
bodeRange = logspace(-6,-1,100); % Define a frequency range that includes all dynamics of the system

% A piece of code that sets the standard interpreter for all plots to 'latex'
% Source: https://nl.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end


