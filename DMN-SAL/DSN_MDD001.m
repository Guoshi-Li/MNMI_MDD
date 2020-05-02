
%=======================================================================================
% Program to estimate effective connectivity (EC) in the Default mode-Salience network 
% This default m-script estimates EC for the MDD subject "MDD001". 
% To estimate a different subject, one needs to change index number for the
% variable "Subject" at the beginning of the program and also inside the
% function FCobj. For example, for the subject MDD010, we need to set
%           Subject = 'MDD010';    
% 
% Guoshi Li (April 30, 2020)
%=======================================================================================

tic

clc;
close all;
clear all;

Subject = 'MDD001';     % MDD subject index
filename = [Subject '_NET1.mat'];   % File name to store the output data

rng(66,'twister'); % Seeds the random number generator

load SC;  % Load structural connectivity

MAX_GEN = 128;   % Maximal number of generations (iterations) for the genetic algorithm
FUN_TOL = 1e-3;  % Functional tolerance

NR = 7;  % Number of regions

% Lower and Upper bound for the model parameters

% Recurrent excitation
MIN_Wee  = 2;
MAX_Wee  = 4;

% Recurrent inhibition
MIN_Wie  = 2;
MAX_Wie  = 4;

% Inter-regional EC
MIN_gc =  -2;
MAX_gc =   2;

% External input
MIN_Spi = 0.2;
MAX_Spi = 0.4;

% Creat lower and upper bounds in vector form
MIN_WEE  = MIN_Wee*ones(1, NR);
MAX_WEE  = MAX_Wee*ones(1, NR);

MIN_WIE  = MIN_Wie*ones(1, NR);
MAX_WIE  = MAX_Wie*ones(1, NR);

MIN_GC = MIN_gc*ones(1,N1);
MAX_GC = MAX_gc*ones(1,N1);

MIN_SPI = MIN_Spi;
MAX_SPI = MAX_Spi;

lb = [MIN_WEE MIN_WIE MIN_GC MIN_SPI]; 
ub = [MAX_WEE MAX_WIE MAX_GC MAX_SPI]; 

NP = length(lb);  % number of parameters to be estimated

% Set optimization parameters

% Use this option if run in a local computer
% opts = optimoptions('ga', 'MaxGenerations', MAX_GEN, 'FunctionTolerance', FUN_TOL, ...
%     'UseParallel',true, 'Display','iter','PlotFcn', @gaplotbestf);

% Use this option if run in a cluster 
opts = optimoptions('ga', 'MaxGenerations', MAX_GEN, 'FunctionTolerance', FUN_TOL, ...
    'UseParallel',true, 'Display','iter');

% Call the ga function to estimate the parameters
[x,fval,eflag,outpt, population, scores] = ga(@FCobj, NP,...
    [],[],[],[],lb,ub,[],opts);


disp('The optimal solution is:');
x
fval

disp('The maximal correlation is:');
f = -FCobj(x)

% Save the output data
save (filename, 'x', 'fval', 'eflag', 'outpt', 'population', 'scores');

poolobj = gcp('nocreate');
delete(poolobj);

disp('Program completed!');

toc

% exit;


%==========================================================================
%           Objective Function
%==========================================================================

function f = FCobj(x)

rng(66,'twister');

Subject = 'MDD001';   % Subject index number

%  Load functional connectivity
ss = ['FC/' Subject];     % 
load (ss);

EFC  = FC1;

% Vectorize the upper part of the FC matrix (since it is symmetric matrix)
UEFC = triu(EFC, 1);
VEFC = UEFC(:);
VEFC(VEFC==0)=[]; 

load SC;             % Structural connectivity
SCM = NET1_SC;       % SC of the DMN-SAL network
 
Flag_Noise = 1;      % Introduce noise to the neural model
Flag_Mean_BOLD = 0;  % If 1: Remove the mean of the neural activity before introducing to the Hemodynamic model

NR = 7;   % Number of ROIs

DT = 10e-3;   % 10 ms; Integration step 
ST = 200;     % Total simulation time in sec

TR  = 2;      % Scan interval (2 sec)
NTR = TR/DT; 

TFC = 20;     % Start time to calculate simulated FC
NFC = TFC/DT;

W_EI0 = 3.0;  % Fixed weight from excitaotry neural population to inhibitory neural population 

% Vectorize the parameters
Wei = W_EI0*ones(NR, 1);
Wee = x(1:NR);   
Wie = x(NR+1:2*NR);
Wgc = x(2*NR+1:2*NR+N1);
SPI = x(end);

GC = zeros(NR,NR);

m = 1;

for i=1:NR
    for j=1:NR
       
       if (MAP1(i,j)==1)
           GC(i,j) = Wgc(m);
           m = m+1;
       end      
        
    end  
end


% Call the function "Model_NEURAL" to calculate the neural activity
[t, X] = Model_NEURAL(NR, DT, ST, SCM, GC, Wee, Wei, Wie, SPI, Flag_Noise);

% Call the function "Model_HEMO" to compute the hemodynamic responses
[tt, BOLD] = Model_HEMO(NR, DT, ST, X, Flag_Mean_BOLD);


SBOLD = BOLD(NFC:NTR:end, :);  % Downsample the BOLD signals

SFC = corrcoef(SBOLD);    % Simulated FC

% Vectorize the upper part of the simulated FC matrix (since it is symmetric matrix)
USFC = triu(SFC, 1);
VSFC = USFC(:);
VSFC(VSFC==0)=[]; 

% Pearson's correlation between the empirical FC and simulated FC
CO = corrcoef(VEFC, VSFC);
COF = CO(1,2);

% Value of the objective function to be minimized
f = -COF;


end


