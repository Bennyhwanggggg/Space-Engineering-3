%% User Inputs
% This function allows user to input their desired optimisation technique
%
% Author Kuan Chun Hwang

function [FDS, CDS, SOD, BFGS, inclination, RAAN_input, vg, bg, mer, nomer] = userInput()

% GUI title
dlg_title = 'User Input';

% Good or bad estimate
prompt1 = {'Good or bad estimate (Good/Bad)'};
% Use merit function or not
prompt2 = {'Merit function? (Yes/No)'};
% Choses the Hessian Update method
prompt3 = {'Second order differential Hessian Approximation (1) or Quasi Newton BFGS Hessian Approximation (2) Note: chosing 2 requires Merit function and Central Differencing Scheme'};
% Choses between Forward differencing or Central differencing
prompt4 = {'Forward differencing (1) or Central Differencing (2)'};
% Input desired inclination
prompt5 = {'Inclination angle (deg) Note: For non-zero inclination, it is only accurate with central differencing, merit function and BFGS'};
% Input desired RAAN
% prompt6 = {'RAAN (deg)'};


% Create GUI
% User_Choices = inputdlg([prompt1,prompt2,prompt3,prompt4,prompt5,prompt6],dlg_title);
User_Choices = inputdlg([prompt1,prompt2,prompt3,prompt4,prompt5],dlg_title);
% Sort out user choices and pass on the information
SOD = strcmp('1',User_Choices(3));
BFGS = strcmp('2',User_Choices(3));
Fd = strcmp('1',User_Choices(4));
Cd = strcmp('2',User_Choices(4));
good = strcmp('Good',User_Choices(1));
bad = strcmp('Bad',User_Choices(1));
Merit = strcmp('Yes',User_Choices(2));
NoMerit = strcmp('No',User_Choices(2));
if Fd == 1
    FDS = 1;
    CDS = 2;
elseif Cd == 1
    FDS = 2;
    CDS = 1;
end
if SOD == 1
    SOD = 1;
    BFGS = 2;
elseif BFGS == 1
    SOD = 2;
    BFGS = 1;
end
if good == 1
    vg = 1;
    bg = 2;
elseif bad == 1
    vg = 2;
    bg = 1;
end
if Merit == 1
    mer = 1;
    nomer = 2;
elseif NoMerit == 1
    mer = 2;
    nomer = 1;
end
% Extract inclination input
inclination = deg2rad(str2num(User_Choices{5}));
% Extract RAAN input
% RAAN_input = deg2rad(str2num(User_Choices{6}));
RAAN_input = 0;
end