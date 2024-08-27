


%% THIS SAMPEL CODE EXPLORES THE INTRINSIC SPIN OF ELASTIC WAVE          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFER CHAPTER 3 Eq. 3.61.1 and 3.62.2 in Book Metamaterial in Topological Acoustics, CRC Press, 2023
% Code written by Prof. Sourav Banerjee Date 06/24/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Explore SpinState
SpinState = 3; % 1 for P-P Interaction, 2 for S-S interaction, 3 for P-S Hybrid Interaction
SpinPlotAlong =0; % Angle in degree
%t1=44.59;
t2=44.59;

t1=69.14;
%t2=-69.14;


%% WAVE - AMPLITUDES 
%Amplitude of P-waves
A1 = 0.1866; % amplitude of L-Wave 1
A2 = 0.1866; % amplitude of L-Wave 2
% Amplitude of S-waves
B1 = 0.6820; % amplitude of T-Wave 1
B2 = -0.6820; % amplitude of T-Wave 2

% C1 = 0.7071; % amplitude of T-Wave 1
% C2 = 0.7071; % amplitude of T-Wave 2

% % Amplitude of P-waves
% A1 = 0.5583; % amplitude of L-Wave 1
% A2 = -0.5583; % amplitude of L-Wave 2
% 
% % Amplitude of S-waves
% B1 = 0.4340; % amplitude of T-Wave 1
% B2 = 0.4340; % amplitude of T-Wave 2
%% WAVE - VELOCITIES
% Velocities of P-wave
Cp1 = 6191.4; % m/sec wave velocity of P-wave 1
Cp2 = 6191.4; % m/sec wave velocity of P-wave 2
% Velocities of S-wave
Cs1 = 3095.7; % m/sec wave velocity of S-wave 1
Cs2 = 3095.7; % m/sec wave velocity of S-wave 2
%% FREQUENCY & DENSITY
w=1003000;  % in Hz frequency
rho=2700;    % in Kg/m^3
SF=rho*w/2; % a Spin Factor 
%% WAVENUMBER
% Wavenumber of P-waves - Two wave interaction
kp1 = w/Cp1; % wavenumber of P-wave 1
kp2 = w/Cp2;  % wavenumber of P-wave 2
%spin fiziud of S0 mode overwrite 
% kp2=1339.45/(2*pi);
% kp1=1339.45/(2*pi);
% Wavenumber of S-wave 
ks1 = w/Cs1; % wavenumber of wave 1
ks2 = w/Cs2;  % wavenumber of wave 2

%% 1-D SPATIAL DISCREATIZATION 
D=0.002;
x_min = 0;
x_max = 0.25; % total length along the plate axis 250 mm
num_points_x = 1000; % number of points in the x direction
del_x=(x_max-x_min)/num_points_x; % in m
Length=x_max; % in m
%% WAVE PROPAGATION DIRETION OF WAVE 1 and WAVE 2 & WAVEVECTOR
theta_1 = t1-SpinPlotAlong ; % wave propagation diraction of P-wave 1
theta_2 = t2-SpinPlotAlong; % wave propagation diraction of P-wave 2

% P-wave Vecotors for wave 1 and wave 2
kp1_1=kp1*cosd(theta_1);  % k1 wave number along x1 for P-wave 1
kp1_2=kp1*sind(theta_1);  % k2 wave number along x2 for P-wave 1

kp2_1=kp2*cosd(theta_2);  % k1 wave number along x1 for P-wave 2
kp2_2=kp2*sind(theta_2);  % k2 wave number along x2 for P-wave 2

% S-wave Vecotors for wave 1 and wave 2
ks1_1=ks1*cosd(theta_1);  % k1 wave number along x1 for S-wave 1
ks1_2=ks1*sind(theta_1);  % k2 wave number along x2 for S-wave 1

ks2_1=ks2*cosd(theta_2);  % k1 wave number along x1 for S-wave 2
ks2_2=ks2*sind(theta_2);  % k2 wave number along x2 for S-wave 2

%% WAVE POLARIZATION 
% Wave polarization direction of the P-wave
P_1 = theta_1 ; % wave polarization diraction of P-wave 1 along Theta 1 
P_2 = theta_2; % wave polarization diraction of P-wave 2 along Theta 2

% Wave polarization direction of the S-wave 
S_1=90+theta_1; %wave polarization diraction of S-wave 1 Orthogonal to Theta 1
S_2=-(90+theta_2); %wave polarization diraction of S-wave 2 Orthogonal to Theta 2

%% SPATIAL DISCRETIZATION
%i=sqrt(-1);

x1=0:del_x:Length; % discretize the x1 axis 
x2=0:del_x:Length; % discretize the x2 axis 

x1_plot=(zeros(length(x2),1))';
x2_plot=(zeros(length(x1),1))';

x3=(zeros(length(x1),1))';

%% WAVE POTENTIALS 
% P-wave potential 
fi1_exp = A1.*exp(1i*((kp1_1.*x1)+(kp1_2.*x2)+(kp1.*cosd(90).*x3)));
fi2_exp = A2.*exp(1i*((kp2_1.*x1)+(kp2_2.*x2)+(kp2.*cosd(90).*x3)));
% SV-wave potential 
Zi1_exp = B1.*exp(1i*((ks1_1.*x1)+(ks1_2.*x2)+(ks1.*cosd(90).*x3)));
Zi2_exp = B2.*exp(1i*((ks2_1.*x1)+(ks2_2.*x2)+(ks2.*cosd(90).*x3)));

%% WAVE DISPLACEMENT AMPLITUDE
% P-wave

u1_fi1 = 1i*kp1_1*fi1_exp; % u1 = dfi/dx1 & u2 = dfi/dx2
u2_fi1 = 1i*kp1_2*fi1_exp; % No minus sign V.Imp

u1_fi2 = 1i*kp2_1*fi2_exp; % u1 = dfi/dx1 & u2 = dfi/dx2
u2_fi2 = 1i*kp2_2*fi2_exp; % No minus sign V.Imp

u1_Zi1 = 1i*ks1_2*Zi1_exp; % See the reverse: u1 = dZi/dx2 & u2 = -dZi/dx1
u2_Zi1 = -1i*ks1_1*Zi1_exp; % Please mind the minus sign V.Imp

u1_Zi2 = 1i*ks2_2*Zi2_exp; % See the reverse: u1 = dZi/dx2 & u2 = -dZi/dx1
u2_Zi2 = -1i*ks2_1*Zi2_exp; % Please mind the minus sign V.Imp


%% SUPERPOSED WAVE DISPLACEMENT ALONG X1 and X2 DIRECTIONS
if (SpinState ==1)
       %% 2-P-wave Superposition
        u1=u1_fi1+u1_fi2;  % P-wave 1+2 disp Along x1
        u2=u2_fi1+u2_fi2;  % P-wave 1+2 disp Along x2
elseif (SpinState==2)
        %% 2-S-wave Superposition 
        u1=u1_Zi1+u1_Zi2; % S-wave 1+2 disp Along x1
        u2=u2_Zi1+u2_Zi2; % S-wave 1+2 disp Along x2
elseif(SpinState==3)
        %% 1P-1S-Hybrid wave P along theta1 so use fi1, S along theta2 so use Zi2
        u1=u1_Zi1+u1_fi2;  
        u2=u2_Zi1+u2_fi2;  
        % u1=u1_fi1+u1_Zi2;  
        % u2=u2_fi1+u2_Zi2;  


end
%%

%% EXPLORING THE INTRINSIC SPIN ONLY X1-X2 PLANE
%% 2-P-wave
sample = size(x1,2);                % At every position along X
Spin_imag=zeros(sample,1);
for i = 1:sample
    ULcross = u1(i)*conj(u2(i))-u2(i)*conj(u1(i));
    %Spin_real(i)=SF*real(ULcross);
    Spin_imag(i)=SF*imag(ULcross);
    %Spin_abs(i)=SF*abs(ULcross);
end

save('BloodyMatchIntrinsic.mat','u1_fi1','u2_fi1','u1_fi2','u2_fi2');
%% VISUALIZATION OF Polarization and SPIN

figure(1);
    quiver(x1,x2_plot,real(u1),real(u2),'k-'); hold on
    
figure(1);
    quiver(x1,x2_plot,imag(u1),imag(u2),'r-');
%legend('real u','imag u')
axis normal


if (SpinState==1)
    figure; %plot(Spin_L_real,'k','LineWidth',4); hold on

    plot(x1,Spin_imag,'LineWidth',4);
    xlim([-(0.1)*max(x1), max(x1)+(0.1)*max(x1)])
    title('Similar Spin State Between P-P waves')
    % Get the current colormap (jet)
    cmap = colormap('jet');
    %scaled_spin = (Spin_imag - min(Spin_imag)) / (max(Spin_imag) - min(Spin_imag));
    % Compute the color indices based on y values
    num_colors = size(cmap, 1);
    color_indices = round(interp1(linspace(min(Spin_imag), max(Spin_imag), num_colors), 1:num_colors, Spin_imag));
    % Apply colors to the line plot
    hold on;
    for i = 1:length(x1)-1
        plot(x1(i:i+1), Spin_imag(i:i+1), 'Color', cmap(color_indices(i), :), 'LineWidth', 14);
    end
    hold off;
    
    % Add colorbar for reference (optional)
    colorbar;
    clim([min(Spin_imag) max(Spin_imag)]);

elseif(SpinState==2)
    figure;%plot(Spin_L_real,'k','LineWidth',4); hold on
    plot(x1,Spin_imag,'LineWidth',4);
    xlim([-(0.1)*max(x1), max(x1)+(0.1)*max(x1)])
    title('Similar Spin State Between S-S waves')
        cmap = colormap('jet');
    % Compute the color indices based on y values
    num_colors = size(cmap, 1);
    color_indices = round(interp1(linspace(min(Spin_imag), max(Spin_imag), num_colors), 1:num_colors, Spin_imag));
    % Apply colors to the line plot
    hold on;
    for i = 1:length(x1)-1
        plot(x1(i:i+1), Spin_imag(i:i+1), 'Color', cmap(color_indices(i), :), 'LineWidth', 8);
    end
    hold off;
    
    % Add colorbar for reference (optional)
    colorbar;
    clim([min(Spin_imag) max(Spin_imag)]);
elseif (SpinState==3)
    figure;%plot(Spin_L_real,'k','LineWidth',4); hold on
    plot(x1,Spin_imag,'LineWidth',4);
    xlim([-(0.1)*max(x1), max(x1)+(0.1)*max(x1)])
    title('Hybrid Spin State Between P-S waves')
        cmap = colormap('jet');
    % Compute the color indices based on y values
    num_colors = size(cmap, 1);
    color_indices = round(interp1(linspace(min(Spin_imag), max(Spin_imag), num_colors), 1:num_colors, Spin_imag));
    % Apply colors to the line plot
    hold on;
    for i = 1:length(x1)-1
        plot(x1(i:i+1), Spin_imag(i:i+1), 'Color', cmap(color_indices(i), :), 'LineWidth', 8);
    end
    hold off;
    
    % Add colorbar for reference (optional)
    colorbar;
    clim([min(Spin_imag) max(Spin_imag)]);
end         
%% END
