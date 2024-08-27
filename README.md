### QuantumSpinGuidedWave

# Citation Requirements

Citation:
If you use this code in your research or any other work, please cite the following paper and repository:
•	Paper Citation:
[Author(s)], "[Title of Your Paper]," [Journal/Conference Name], [Year], [DOI or URL if available].
•	Code Citation:
[Author(s)], "[Title of Your GitHub Repository]," GitHub, [Year]. Available at: [URL of the GitHub repository].


Three-Step Code used for publication in Ultrasonics and ...

Codes could be used freely (with or without modification) subjected to proper citation of the work by the author and the paper listed above.

The codes are divided into three steps. 

Step 1: It helps find the dispersion curves and wavemodes in a wave guided with user-selected material properties. The code finds the Wave numbers at given frequencies. 
Step 2: Modal eigenstates of Guided waves are calculated. From Step 1 the eigenvalues are known. For any k some frequency eigenvalues are found in Step 1. Select any mode and its eigenvalue, to plug in Step 2 code, to find the amplitudes of the wave potentials. Use the amplitude to calculate the wave displacement field. 
Step 3: In step 3 First Select Which Spin Characteristic we want to explore. 
% Select RLSH  index = 1 for Rayleigh-Lamb Modes 
% Select RLSH index = 2 for All RL and SH mode 
%% SpinStates
% Please note for RLSH = 2 SpinState can only take values between 3 and 10.
% as 1 and 2 are not relevant in RLSH SpinState. 

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

Plug the Wave potential amplitudes from Step 2 and proceed to explore the spin state. 
Call 
    % % RL mode only For Rayleigh-Lamb waves only 
    GWUTSpinFun(k,w,A1,A2,B1,B2,SpinState, ...
                                           E, nu, rho,...
                                           d, xl, xd, x_min, x_max, ...
                                           num_points_x, num_points_y,...
                                           plotindex);
                                           
                                           For Rayleigh-Lamb waves only 
Call 
    % RL & SH mode only : 
    GWUTSpinFun_RLSH(k,kh,w,A1,A2,B1,B2,C1,C2,SpinState, ...
                                           E, nu, rho,...
                                           d, xl, xd, x_min, x_max, ...
                                           num_points_x, num_points_y,...
                                           plotindex);

                                           Explore Hybrid Spin States when RL and SH are both present in a 3D waveguide. 
end
