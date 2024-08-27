%% Step 3 Program :  Plug the Amplitudes and Polarity Directions from Step 2 : Find the Spin States  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THIS SAMPEL CODE EXPLORES THE INTRINSIC SPIN OF ELASTIC GUIDED WAVE          %%
%% WHEN Eigen Values are Known. Find Eigen Vector (Polarity) and Find Spin State of a MODE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
%% First Select Which Spin Characteristic we want 
% Select RLSH  index = 1 for Rayleigh-Lamb Modes 
% Select RLSH index = 2 for All RL and SH mode  
%% SpinStates
% Please note for RLSH = 2 SpinState can only take values between 3 and 10.
% as 1 and 2 is not relevant in RLSH SpinState. 

% SpinState = 1 ; % s_fi
% SpinState = 2 ; % s_zi

% SpinState = 3 ; % s_fiziud
% SpinState = 4 ; % s_fizidu
% SpinState = 5 ; % s_fiziuu
% SpinState = 6 ; % s_fizidd
% SpinState = 7 ; % s_ziz3ud
% SpinState = 8 ; % s_ziz3du
% SpinState = 9 ; % s_ziz3uu
% SpinState = 10 ; % s_ziz3dd
%%
RLSH = 2; 
SpinState = 10;
%% Plot along x1 or x2 

plotindex = 3; 

xl=1;
xd=0;

x_min = 0;
x_max = 0.25; % total length along the plate axis 250 mm
num_points_x = 5000; % number of points in the x direction
num_points_y = 1000; % number of points in the y direction

del_x=(x_max-x_min)/num_points_x; % in m
Length=x_max; % in m

%% WAVE - AMPLITUDES : Get from Step 2
%Amplitude of P-waves
A1 = 0.1866; % amplitude of L-Wave 1
A2 = 0.1866; % amplitude of L-Wave 2
% Amplitude of S-waves
B1 = 0.6820; % amplitude of T-Wave 1
B2 = -0.6820; % amplitude of T-Wave 2

C1 = 0.7071; % amplitude of T-Wave 1
C2 = 0.7071; % amplitude of T-Wave 2
%% Wave POLARITY in deg : 
% Its Redundant as could be found inside the GWUTSpinFun_RLSH or GWUTSpinFun
tp=0;
tsv=50.76;
tsh=52.36;
%% Wave number and Frequency from Step 1.3 after Plotting Dispersion
k=1348.05;         % selected from k-w solution  It has 2*pi factor in it 
kh=1301.35;       % selected from k-w solution  It has 2*pi factor in it
w=1050000*2*pi;% As found from dispersion plot 
omg = w/(2*pi);

%% Material Properties and DENSITY
                        % %% Define Geometry of problem
                        % 
                        % h=1.0e-3;       % average height
                        % d=2*h;
                        % D=1*d;        % Length of period
                        % 
                        % kmax = pi/D;
                        % freq_start = 1e3;
                        % delfreq = 2e3;
                        % freq_end = 2e6;
                        % Corrugate_Coeff=0;
                        % e=Corrugate_Coeff*h;
                        % 
                        % %% material properties of problem
                        % E=69e9;
                        % nu=1/3;
                        % rho=2700;
                        % 
                        % lam=E*nu/((1-2*nu)*(1+nu));
                        % mu=E/(2*(1+nu));
                        % Cp=sqrt((lam+2*mu)/rho);
                        % Cs=sqrt(mu/rho);
                        % S=lam/mu;
                        % T=(4*pi^2*e/D^2);


%% OR LOAD THE Variables from a Prerun Data File from Step 1.1 
load ('wk_PSV_Disp.mat','h','d','Cp','Cs','E','nu','rho');
rho;    % in Kg/m^3
E;       % Pa
Cp;
Cs;


%% What Spin State we want to explore

%% SpinStates
% Please note for RLSH = 2 SpinState can only take values between 3 and 10.
% as 1 and 2 is not relevant in RLSH SpinState. 

% SpinState = 1 ; % s_fi
% SpinState = 2 ; % s_zi

% SpinState = 3 ; % s_fiziud
% SpinState = 4 ; % s_fizidu
% SpinState = 5 ; % s_fiziuu
% SpinState = 6 ; % s_fizidd
% SpinState = 7 ; % s_ziz3ud
% SpinState = 8 ; % s_ziz3du
% SpinState = 9 ; % s_ziz3uu
% SpinState = 10 ; % s_ziz3dd


if (RLSH ==1)
    % % RL mode only 
    GWUTSpinFun(k,w,A1,A2,B1,B2,SpinState, ...
                                           E, nu, rho,...
                                           d, xl, xd, x_min, x_max, ...
                                           num_points_x, num_points_y,...
                                           plotindex);
elseif (RLSH ==2)
    % RL & SH mode only 
    GWUTSpinFun_RLSH(k,kh,w,A1,A2,B1,B2,C1,C2,SpinState, ...
                                           E, nu, rho,...
                                           d, xl, xd, x_min, x_max, ...
                                           num_points_x, num_points_y,...
                                           plotindex);
end