%% Step 3 Program :  Plug the Amplitudes and Polarity Directions from Step 2 : Find the Spin States  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THIS SAMPEL CODE EXPLORES THE INTRINSIC SPIN OF ELASTIC GUIDED WAVE          %%
%% WHEN Eigen Values are Known. Find Eigen Vector (Polarity) and Find Spin State of a MODE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
%% First Select Which Spin Characteristic we want 
% Select RLSH  index = 1 for Rayleigh-Lamb Modes 
% Select RLSH index = 2 for All RL and SH mode  
%% SpinStates
% Please note for RLSH = 2 SpinState can only take values between 3 and 10.
% as 1 and 2 is not relevant in RLSH SpinState. 

%% P-SV  for RL modes
% SpinState = 1 ; % s_fiud
% SpinState = 2 ; % s_zi3ud
% SpinState = 3 ; % s_fizi3ud
% SpinState = 4 ; % s_fizi3du
% SpinState = 5 ; % s_fizi3uu
% SpinState = 6 ; % s_fizi3dd

%% RL-SH Hybrid mode
% SpinState = 7 ; % s_fizi1ud
% SpinState = 8 ; % s_fizi1du
% SpinState = 9 ; % s_fizi1uu
% SpinState = 10 ; % s_fizi1dd
% SpinState = 11 ; % s_zi1z3ud
% SpinState = 12 ; % s_zi1z3du
% SpinState = 13 ; % s_zi1z3uu
% SpinState = 14 ; % s_zi1z3dd

%% Plot along x1 or x2 

plotindex = 1; 



x_min = 0;
x_max = 1; % total length along the plate axis 500 mm
num_points_x = 5000; % number of points in the x direction
num_points_y = 1000; % number of points in the y direction

del_x=(x_max-x_min)/num_points_x; % in m
Length=x_max; % in m

%% WAVE - AMPLITUDES : Get from Step 2
%Amplitude of P-waves
A1 = 0.6601; % amplitude of L-Wave 1
A2 = 0.6601; % amplitude of L-Wave 2
% Amplitude of S-waves
B1 = 0.2535; % amplitude of T-Wave 1
B2 = -0.2535; % amplitude of T-Wave 2

C1 = 0.7071; % amplitude of T-Wave 1
C2 = -0.7071; % amplitude of T-Wave 2
%% Wave POLARITY in deg : 
% Its Redundant as could be found inside the GWUTSpinFun_RLSH or GWUTSpinFun
% tp=0;
% tsv=50.76;
% tsh=52.36;
%% Wave number and Frequency from Step 1.3 after Plotting Dispersion
k=295.25;         % selected from k-w solution  It has 2*pi factor in it 
kh=509.443;       % selected from k-w solution  It has 2*pi factor in it
w=251000*2*pi;% As found from dispersion plot 
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

xl=1;
xd=+h;
%% What Spin State we want to explore

%% SpinStates
% Please note for RLSH = 2 SpinState can only take values between 3 and 10.
% as 1 and 2 is not relevant in RLSH SpinState. 

% SpinState = 1 ; % s_fi
% SpinState = 2 ; % s_zi

%% SpinStates
% SpinState = 1 ; % s_fiud
% SpinState = 2 ; % s_zi3ud
% SpinState = 3 ; % s_fizi3ud
% SpinState = 4 ; % s_fizi3du
% SpinState = 5 ; % s_fizi3uu
% SpinState = 6 ; % s_fizi3dd
%% 
% SpinState = 7 ; % s_fizi1ud
% SpinState = 8 ; % s_fizi1du
% SpinState = 9 ; % s_fizi1uu
% SpinState = 10 ; % s_fizi1dd
% SpinState = 11 ; % s_zi1z3ud
% SpinState = 12 ; % s_zi1z3du
% SpinState = 13 ; % s_zi1z3uu
% SpinState = 14 ; % s_zi1z3dd

%%
Total_SAM = cell(14,1);
TotalNormSAM = zeros(num_points_x+1,3);
for wave = 1:2

    RLSH = wave;
    %% 
    if (RLSH ==1) 
        %%
        for SpinState = 1:6
            % % RL mode only 
            Spin_imag = GWUTSpinFun_SHM(k,w,A1,A2,B1,B2,SpinState, ...
                                                   E, nu, rho,...
                                                   d, xl, xd, x_min, x_max, ...
                                                   num_points_x, num_points_y,...
                                                   plotindex);

            Total_SAM{SpinState}=  Spin_imag;
            TotalNormSAM = TotalNormSAM + Spin_imag;
        end
    
    elseif (RLSH ==2)
        % RL & SH mode only 
        %%
        for SpinState = 7:14

            Spin_imag = GWUTSpinFun_RLSH_SHM(k,kh,w,A1,A2,B1,B2,C1,C2,SpinState, ...
                                                   E, nu, rho,...
                                                   d, xl, xd, x_min, x_max, ...
                                                   num_points_x, num_points_y,...
                                                   plotindex);
            Total_SAM{SpinState} = Spin_imag;
            TotalNormSAM = TotalNormSAM + Spin_imag;
        end
    end


end

% for i =1:length(TotalNormSAM(:,1))
% 
%     SAMMag = sqrt(TotalNormSAM(i,1)^2+TotalNormSAM(i,2)^2+TotalNormSAM(i,3)^2) ;
%     TotalSAM(i,1)=TotalNormSAM(i,1)/SAMMag;
%     TotalSAM(i,2)=TotalNormSAM(i,2)/SAMMag;
%     TotalSAM(i,3)=TotalNormSAM(i,3)/SAMMag;
% end
 %% Plotting of Normalized Final SAM
del_x=(x_max-x_min)/num_points_x; % in m
Length=x_max; % in m
x1=0:del_x:Length;   % discretize the x1 axis 

Spin_1 = TotalNormSAM(:,1);
Spin_2 = TotalNormSAM(:,2);
Spin_3 = TotalNormSAM(:,3);

close all;
try
figure; %plot(Spin_L_real,'k','LineWidth',4); hold on
    plot(x1,Spin_1,'LineWidth',4);
    xlim([-(0.1)*max(x1), max(x1)+(0.1)*max(x1)])
    
    title('Similar Spin State Between P-P waves')
    % Get the current colormap (jet)
    cmap = colormap('jet');
    %scaled_spin = (Spin_imag - min(Spin_imag)) / (max(Spin_imag) - min(Spin_imag));
    % Compute the color indices based on y values
    num_colors = size(cmap, 1);
    color_indices = round(interp1(linspace(min(Spin_1), max(Spin_1), num_colors), 1:num_colors, Spin_1));
    % Apply colors to the line plot
    hold on;
    for i = 1:length(x1)-1
        plot(x1(i:i+1), Spin_1(i:i+1), 'Color', cmap(color_indices(i), :), 'LineWidth', 6);
    end
    hold off;
    % Add colorbar for reference (optional)
    colorbar;
    clim([min(Spin_1) max(Spin_1)]);
catch ME
    display(ME.message);
end

try
figure;%plot(Spin_L_real,'k','LineWidth',4); hold on
    plot(x1,Spin_2,'LineWidth',4);
    xlim([-(0.1)*max(x1), max(x1)+(0.1)*max(x1)])
    
    title('Similar Spin State Between S-S waves')
        cmap = colormap('jet');
    % Compute the color indices based on y values
    num_colors = size(cmap, 1);
    color_indices = round(interp1(linspace(min(Spin_2), max(Spin_2), num_colors), 1:num_colors, Spin_2));
    % Apply colors to the line plot
    hold on;
    for i = 1:length(x1)-1
        plot(x1(i:i+1), Spin_2(i:i+1), 'Color', cmap(color_indices(i), :), 'LineWidth', 6);
    end
    hold off;
    % Add colorbar for reference (optional)
    colorbar;
    clim([min(Spin_2) max(Spin_2)]);
catch ME
    display(ME.message)
end


 %% Mathematical Derivation Says Spin_3 is zero i.e. x_3 component of the spin. 
try
figure;
    plot(x1,Spin_3,'LineWidth',4);
    xlim([-(0.1)*max(x1), max(x1)+(0.1)*max(x1)])
    
    title('Hybrid Spin State Between P-S waves')
        cmap = colormap('jet');
    % Compute the color indices based on y values
    num_colors = size(cmap, 1);
    color_indices = round(interp1(linspace(min(Spin_3), max(Spin_3), num_colors), 1:num_colors, Spin_3));
    % Apply colors to the line plot
    hold on;
    for i = 1:length(x1)-1
        plot(x1(i:i+1), Spin_3(i:i+1), 'Color', cmap(color_indices(i), :), 'LineWidth', 6);
    end
    hold off;

    % Add colorbar for reference (optional)
    colorbar;
    clim([min(Spin_3) max(Spin_3)]);
         
catch ME
    display(ME.message);
end