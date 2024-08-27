


%% THIS SAMPEL CODE EXPLORES THE INTRINSIC SPIN OF ELASTIC WAVE          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% REFER CHAPTER 3 Eq. 3.61.1 and 3.62.2
% Code written by Prof. Sourav Banerjee Date 06/24/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GWUTSpinFun_RLSH(k,kh,w,Au,Ad,Bu,Bd,Cu,Cd,SpinState, ...
                                       E, nu, rho,...
                                       d, xl, xd, x_min, x_max, ...
                                       num_points_x, num_points_y,...
                                       plotindex)
SpinPlotAlong = 0; % Angle in degree
%% WAVE - AMPLITUDES 
% Amplitude of P-waves
Ad=real(Ad); % amplitude of L-Wave Down Going
Au=real(Au); % amplitude of L-Wave Up Going
% % Amplitude of SV-waves
Bd=real(Bd); % amplitude of T-Wave 1
Bu=real(Bu); % amplitude of T-Wave 2
% % Amplitude of SV-waves
Cd=real(Cd); % amplitude of T-Wave 1
Cu=real(Cu); % amplitude of T-Wave 2

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


theta_svu = real(acosd(k/ks)); % wave propagation diraction of S-wave 2
theta_svd = -real(acosd(k/ks)); % wave propagation diraction of S-wave 2

theta_shu = real(acosd(kh/ks)); % wave propagation diraction of S-wave 2
theta_shd = -real(acosd(kh/ks)); % wave propagation diraction of S-wave 2


%%
%% WAVE POLARIZATION DIRETION OF WAVE 1 and WAVE 2 & WAVEVECTOR
theta_p1 = theta_pu-SpinPlotAlong ; % wave propagation diraction of P-wave 1
theta_p2 = theta_pd-SpinPlotAlong; % wave propagation diraction of P-wave 2

theta_sv1 = theta_svu-SpinPlotAlong ; % wave propagation diraction of P-wave 1
theta_sv2 = theta_svd-SpinPlotAlong; % wave propagation diraction of P-wave 2

theta_sh1 = theta_shu-SpinPlotAlong ; % wave propagation diraction of P-wave 1
theta_sh2 = theta_shd-SpinPlotAlong; % wave propagation diraction of P-wave 2

% P-wave Vecotors for wave 1 and wave 2
kp1_1=kp*cosd(theta_p1);  % k1 wave number along x1 for P-wave 1
kp1_2=kp*sind(theta_p1);  % k2 wave number along x2 for P-wave 1

kp2_1=kp*cosd(theta_p2);  % k1 wave number along x1 for P-wave 2
kp2_2=kp*sind(theta_p2);  % k2 wave number along x2 for P-wave 2

% SV-wave Vecotors for wave 1 and wave 2
ksv1_1=ks*cosd(theta_sv1);  % k1 wave number along x1 for S-wave 1
ksv1_2=ks*sind(theta_sv1);  % k2 wave number along x2 for S-wave 1

ksv2_1=ks*cosd(theta_sv2);  % k1 wave number along x1 for S-wave 2
ksv2_2=ks*sind(theta_sv2);  % k2 wave number along x2 for S-wave 2

% SH-wave Vecotors for wave 1 and wave 2
ksh1_1=ks*cosd(theta_sh1);  % k1 wave number along x1 for S-wave 1
ksh1_2=ks*sind(theta_sh1);  % k2 wave number along x2 for S-wave 1

ksh2_1=ks*cosd(theta_sh2);  % k1 wave number along x1 for S-wave 2
ksh2_2=ks*sind(theta_sh2);  % k2 wave number along x2 for S-wave 2

% Please Note k1_up = k1_dp = k1_us = k1_ds = k
% for GWUT kpd = kpu ; ksd=ksu ; 
%% 1-D SPATIAL DISCREATIZATION 
% Sp_samprate_x=(x_max-x_min)/num_points_x;
% Sp_samprate_y=(y_max-y_min)/num_points_y;
% 
% X_limit=x_max;
x_min;
x_max; % total length along the plate axis D=0.002 m 
%num_points_x = 1000; % number of points in the x direction
%num_points_y = 500;
del_x=(x_max-x_min)/num_points_x; % in m
Length=x_max; % in m

%% WAVE POLARIZATION 
% Wave polarization direction of the P-wave
L_theta_d = -theta_pu ; % wave polarization diraction of P-wave 1
L_theta_u = theta_pu; % wave polarization diraction of P-wave 2
% Wave polarization direction of the SV-wave 
T_theta_d=-(90+theta_svu); %wave polarization diraction of S-wave 1
T_theta_u=90+theta_svu; %wave polarization diraction of S-wave 2

% Wave polarization direction of the SV-wave
% Orthogonal to the x1-x2 plane 
%% SPATIAL DISCRETIZATION
i=sqrt(-1);
% 
x1=0:del_x:Length;   % discretize the x1 axis 
x2 = xd*ones(length(x1),1)';
x2_plot=(zeros(length(x1),1))';

%x2 =zeros(length(x1))
%x2=(0:Sp_samprate:X_limit)+xd; % discretize the x2 axis 

x3=(zeros(length(x1),1))';
%% WAVE POTENTIALS 

% P-wave potential 
fi1_exp = Au.*exp(1i*((kp1_1.*x1)+(kp1_2.*x2)+(kp.*cosd(90).*x3)));
fi2_exp = Ad.*exp(1i*((kp2_1.*x1)+(kp2_2.*x2)+(kp.*cosd(90).*x3)));
% SV-wave potential 
Zi1_exp = Bu.*exp(1i*((ksv1_1.*x1)+(ksv1_2.*x2)+(ks.*cosd(90).*x3)));
Zi2_exp = Bd.*exp(1i*((ksv2_1.*x1)+(ksv2_2.*x2)+(ks.*cosd(90).*x3)));
% SH-wave potential 
Z31_exp = Cu.*exp(1i*((ksh1_1.*x1)+(ksh1_2.*x2)+(ks.*cosd(90).*x3)));
Z32_exp = Cd.*exp(1i*((ksh2_1.*x1)+(ksh2_2.*x2)+(ks.*cosd(90).*x3)));

%% WAVE DISPLACEMENT AMPLITUDE
% P-wave

u1_fiu = 1i*kp1_1*fi1_exp; % u1 = dfi/dx1 & u2 = dfi/dx2
u2_fiu = 1i*kp1_2*fi1_exp; % No minus sign V.Imp

u1_fid = 1i*kp2_1*fi2_exp; % u1 = dfi/dx1 & u2 = dfi/dx2
u2_fid = 1i*kp2_2*fi2_exp; % No minus sign V.Imp

% SV wave
u1_Ziu = 1i*ksv1_2*Zi1_exp; % See the reverse: u1 = dZi/dx2 & u2 = -dZi/dx1
u2_Ziu = -1i*ksv1_1*Zi1_exp; % Please mind the minus sign V.Imp

u1_Zid = 1i*ksv2_2*Zi2_exp; % See the reverse: u1 = dZi/dx2 & u2 = -dZi/dx1
u2_Zid = -1i*ksv2_1*Zi2_exp; % Please mind the minus sign V.Imp

% SH wave
u3_Z3u=Z31_exp;
u3_Z3d=Z32_exp;  % No Minus sign V.Imp. 

% % Finding Conjugates if necessary - Done for each Spin State Later 
% u1_fiu_c = conj(u1_fiu);
% u2_fiu_c = conj(u2_fiu);
% 
% u1_fid_c = conj(u1_fid);
% u2_fid_c = conj(u2_fid);
% 
% u1_Ziu_c = conj(u1_Ziu);
% u2_Ziu_c = conj(u2_Ziu);
% 
% u1_Zid_c = conj(u1_Zid);
% u2_Zid_c = conj(u2_Zid);
% 
% u3_Z3u_c = conj(u3_Z3u);
% u3_Z3d_c = conj(u3_Z3d);




%% Spin Matrices Already has the Imaginary Part
S1=[0 0 0 ; 0 0 1i ; 0 -1i 0];
S2=[0 0 1i ; 0 0 0 ; -1i 0 0];
S3=[0 -1i 0 ; 1i 0 0 ; 0 0 0];

%% SUPERPOSED WAVE DISPLACEMENT ALONG X1 and X2 DIRECTIONS

if(SpinState==3) % SfiZ3uu p1 = fi, p2=Z3
        u1=u1_fiu;  
        u2=u2_fiu;
        u3=u3_Z3u;
        

elseif(SpinState==4) %SfiZ3ud
        u1=u1_fiu;  
        u2=u2_fiu;
        u3=u3_Z3d;

elseif(SpinState==5) %SfiZ3du
        u1=u1_fid;  
        u2=u2_fid;  
        u3=u3_Z3u;

elseif(SpinState==6) %SfiZ3dd
        u1=u1_fid;  
        u2=u2_fid; 
        u3=u3_Z3d;

elseif(SpinState==7) %SZiZ3uu
        u1=u1_Ziu;  
        u2=u2_Ziu;
        u3=u3_Z3u;

elseif(SpinState==8) %SZiZ3ud
        u1=u1_Ziu;  
        u2=u2_Ziu;
        u3=u3_Z3d;

elseif(SpinState==9) %SZiZ3du
        u1=u1_Zid;  
        u2=u2_Zid;
        u3=u3_Z3u;

elseif(SpinState==10) %SZiZ3dd
        u1=u1_Zid;  
        u2=u2_Zid;
        u3=u3_Z3d;

end
%% Find Conjugate of the Displacements for Respective Spin State 
        u1_c=conj(u1);
        u2_c=conj(u2);
        u3_c=conj(u3);

        

%% EXPLORING THE INTRINSIC SPIN ONLY X1-X2 PLANE
%% 2-P-wave
sample = size(x1,2);                % At every position along X
Spin_1=zeros(sample,1);
Spin_2=zeros(sample,1);
Spin_3=zeros(sample,1);

for i = 1:sample
    
        % Mind the nomenclature : First Potential is p1, second ptential is p2
        u_p1 = [u1(i) u2(i) 0];
        u_p2 = [0 0 u3(i)];
        u_p1_c = [u1_c(i) u2_c(i) 0];
        u_p2_c = [0 0 u3_c(i)];

    ULcross_1 = (u_p1_c*S1*u_p2.')+(u_p2_c*S1*u_p1.');
    ULcross_2 = (u_p1_c*S2*u_p2.')+(u_p2_c*S2*u_p1.');
    ULcross_3 = (u_p1_c*S3*u_p2.')+(u_p2_c*S3*u_p1.');
    
    Spin_1(i)=SF*(ULcross_1);
    Spin_2(i)=SF*(ULcross_2);
    Spin_3(i)=SF*(ULcross_3);
    
    %Spin_abs(i)=SF*abs(ULcross);
end


%% VISUALIZATION OF Polarization and SPIN

% figure(1);
%     quiver(x1,x2_plot,real(u1),real(u2),'k-'); hold on
% 
% figure(1);
%     quiver(x1,x2_plot,real(u1),real(u3),'r-');
% %legend('real u','imag u')
% axis normal


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
        plot(x1(i:i+1), Spin_1(i:i+1), 'Color', cmap(color_indices(i), :), 'LineWidth', 16);
    end
    hold off;
    % Add colorbar for reference (optional)
    colorbar;
    clim([min(Spin_1) max(Spin_1)]);


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
        plot(x1(i:i+1), Spin_2(i:i+1), 'Color', cmap(color_indices(i), :), 'LineWidth', 16);
    end
    hold off;
    % Add colorbar for reference (optional)
    colorbar;
    clim([min(Spin_2) max(Spin_2)]);
 %% Mathematical Derivation Says Spin_3 is zero i.e. x_3 component of the spin. 
% figure;
%     plot(x1,Spin_3,'LineWidth',4);
%     xlim([-(0.1)*max(x1), max(x1)+(0.1)*max(x1)])
%     title('Hybrid Spin State Between P-S waves')
%         cmap = colormap('jet');
%     % Compute the color indices based on y values
%     num_colors = size(cmap, 1);
%     color_indices = round(interp1(linspace(min(Spin_3), max(Spin_3), num_colors), 1:num_colors, Spin_3));
%     % Apply colors to the line plot
%     hold on;
%     for i = 1:length(x1)-1
%         plot(x1(i:i+1), Spin_3(i:i+1), 'Color', cmap(color_indices(i), :), 'LineWidth', 4);
%     end
%     hold off;
% 
%     % Add colorbar for reference (optional)
%     colorbar;
%     clim([min(Spin_3) max(Spin_3)]);
         

%% END