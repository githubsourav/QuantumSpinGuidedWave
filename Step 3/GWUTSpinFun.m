


%% THIS SAMPEL CODE EXPLORES THE INTRINSIC SPIN OF ELASTIC WAVE          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% REFER CHAPTER 3 Eq. 3.61.1 and 3.62.2
% Code written by Prof. Sourav Banerjee Date 06/24/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GWUTSpinFun(k,w,Au,Ad,Bu,Bd,SpinState, ...
                                       E, nu, rho,...
                                       d, xl, xd, x_min, x_max, ...
                                       num_points_x, ~,...
                                       plotindex)
SpinPlotAlong = 0; % Angle in degree
%% WAVE - AMPLITUDES 
% Amplitude of P-waves
Ad=real(Ad); % amplitude of L-Wave Down Going
Au=real(Au); % amplitude of L-Wave Up Going
% % Amplitude of S-waves
Bd=real(Bd); % amplitude of T-Wave 1
Bu=real(Bu); % amplitude of T-Wave 2

y_max =+d/2;
y_min = -d/2;
D=d;
%xd = at depth along x2 where spin state is shout
%% Wave Vector 
% k;
%% SpinStates
% SpinState = 1 ; % s_fi
% SpinState = 2 ; % s_zi
% SpinState = 3 ; % s_fiziud
% SpinState = 4 ; % s_fizidu
% SpinState = 5 ; % s_fiziuu
% SpinState = 6 ; % s_fizidd

%% MATERIAL Propoerties
% E;
% nu;
% rho;
%% WAVE - VELOCITIES
lam=E*nu/((1-2*nu)*(1+nu));
mu=E/(2*(1+nu));
Cp=sqrt((lam+2*mu)/rho);
Cs=sqrt(mu/rho);
%% Spin Factor with FREQUENCY & DENSITY

SF=rho*w/2; % a Spin Factor 
%% WAVENUMBER
% Wavenumber of P-waves
kp = w/Cp; % wavenumber of P-wave 1
% Wavenumber of S-wave 
ks = w/Cs; % wavenumber of wave 1
%% WAVE PROPAGATION DIRETION OF WAVE 1 and WAVE 2 & WAVEVECTOR
theta_pu = real(acosd(k/kp)) ; % wave propagation diraction of P-wave 1
theta_pd = -real(acosd(k/kp)) ; % wave propagation diraction of P-wave 2


theta_su = real(acosd(k/ks)); % wave propagation diraction of S-wave 2
theta_sd = -real(acosd(k/ks)); % wave propagation diraction of S-wave 2



%%
%% WAVE POLARIZATION DIRETION OF WAVE 1 and WAVE 2 & WAVEVECTOR
theta_p1 = theta_pu-SpinPlotAlong ; % wave propagation diraction of P-wave 1
theta_p2 = theta_pd-SpinPlotAlong; % wave propagation diraction of P-wave 2

theta_s1 = theta_su-SpinPlotAlong ; % wave propagation diraction of P-wave 1
theta_s2 = theta_sd-SpinPlotAlong; % wave propagation diraction of P-wave 2


% P-wave Vecotors for wave 1 and wave 2
kp1_1=kp*cosd(theta_p1);  % k1 wave number along x1 for P-wave 1
kp1_2=kp*sind(theta_p1);  % k2 wave number along x2 for P-wave 1

kp2_1=kp*cosd(theta_p2);  % k1 wave number along x1 for P-wave 2
kp2_2=kp*sind(theta_p2);  % k2 wave number along x2 for P-wave 2

% S-wave Vecotors for wave 1 and wave 2
ks1_1=ks*cosd(theta_s1);  % k1 wave number along x1 for S-wave 1
ks1_2=ks*sind(theta_s1);  % k2 wave number along x2 for S-wave 1

ks2_1=ks*cosd(theta_s2);  % k1 wave number along x1 for S-wave 2
ks2_2=ks*sind(theta_s2);  % k2 wave number along x2 for S-wave 2

% Please Note k1_up = k1_dp = k1_us = k1_ds = k
% for GWUT kpd = kpu ; ksd=ksu ; 
%% 1-D SPATIAL DISCREATIZATION 
% Sp_samprate_x=(x_max-x_min)/num_points_x;
% Sp_samprate_y=(y_max-y_min)/num_points_y;
% 
% X_limit=x_max;

%num_points_x = 1000; % number of points in the x direction
%num_points_y = 500;
del_x=(x_max-x_min)/num_points_x; % in m
Length=x_max; % in m

%% WAVE POLARIZATION 
% Wave polarization direction of the P-wave
L_theta_d = -theta_pu ; % wave polarization diraction of P-wave 1
L_theta_u = theta_pu; % wave polarization diraction of P-wave 2
% Wave polarization direction of the S-wave 
T_theta_d=-(90+theta_su); %wave polarization diraction of S-wave 1
T_theta_u=90+theta_su; %wave polarization diraction of S-wave 2

%% SPATIAL DISCRETIZATION
i=sqrt(-1);
% 
if (plotindex ==1)
    x1=0:del_x:Length;   % discretize the x1 axis 
    x2 = xd*ones(length(x1),1)';
    x2_plot=(zeros(length(x1),1))';
elseif(plotindex==2)
    x2=y_min:Sp_samprate_y:y_max; % discretize the x1 axis 
    x1 = xl*ones(length(x2),1)';
    x1_plot=(zeros(length(x2),1))';
end
%x2 =zeros(length(x1))
%x2=(0:Sp_samprate:X_limit)+xd; % discretize the x2 axis 

x3=(zeros(length(x1),1))';
%% WAVE POTENTIALS 

% P-wave potential 
fi1_exp = Au.*exp(1i*((kp1_1.*x1)+(kp1_2.*x2)+(kp.*cosd(90).*x3)));
fi2_exp = Ad.*exp(1i*((kp2_1.*x1)+(kp2_2.*x2)+(kp.*cosd(90).*x3)));
% S-wave potential 
Zi1_exp = Bu.*exp(1i*((ks1_1.*x1)+(ks1_2.*x2)+(ks.*cosd(90).*x3)));
Zi2_exp = Bd.*exp(1i*((ks2_1.*x1)+(ks2_2.*x2)+(ks.*cosd(90).*x3)));
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
        u1=u1_fi1+u1_Zi2;  
        u2=u2_fi1+u2_Zi2;
elseif(SpinState==4)
        u1=u1_Zi1+u1_fi2;  
        u2=u2_Zi1+u2_fi2;  
elseif(SpinState==5)
        u1=u1_Zi1+u1_fi1;  
        u2=u2_Zi1+u2_fi1;  
elseif(SpinState==6)
        u1=u1_Zi2+u1_fi2;  
        u2=u2_Zi2+u2_fi2; 
end
save('ShouldMatchSpinFun.mat','u1_fi1','u2_fi1','u1_fi2','u2_fi2');

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
        plot(x1(i:i+1), Spin_imag(i:i+1), 'Color', cmap(color_indices(i), :), 'LineWidth', 8);
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
        plot(x1(i:i+1), Spin_imag(i:i+1), 'Color', cmap(color_indices(i), :), 'LineWidth', 4);
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
        plot(x1(i:i+1), Spin_imag(i:i+1), 'Color', cmap(color_indices(i), :), 'LineWidth', 4);
    end
    hold off;
    
    % Add colorbar for reference (optional)
    colorbar;
    clim([min(Spin_imag) max(Spin_imag)]);
         
elseif (SpinState==4)
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
        plot(x1(i:i+1), Spin_imag(i:i+1), 'Color', cmap(color_indices(i), :), 'LineWidth', 4);
    end
    hold off;
    
    % Add colorbar for reference (optional)
    colorbar;
    clim([min(Spin_imag) max(Spin_imag)]);
     
elseif (SpinState==5)
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
        plot(x1(i:i+1), Spin_imag(i:i+1), 'Color', cmap(color_indices(i), :), 'LineWidth', 4);
    end
    hold off;
    
    % Add colorbar for reference (optional)
    colorbar;
    clim([min(Spin_imag) max(Spin_imag)]);
         
elseif (SpinState==6)
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
        plot(x1(i:i+1), Spin_imag(i:i+1), 'Color', cmap(color_indices(i), :), 'LineWidth', 4);
    end
    hold off;
    
    % Add colorbar for reference (optional)
    colorbar;
    clim([min(Spin_imag) max(Spin_imag)]);
end     
%% END