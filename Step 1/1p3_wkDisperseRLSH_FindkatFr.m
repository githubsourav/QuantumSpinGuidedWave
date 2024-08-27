%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1.3 Program for RLSH :  Plot the Dispersion Curves.                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THIS SAMPEL CODE EXPLORES THE INTRINSIC SPIN OF ELASTIC GUIDED WAVE          %%
%% WHEN Eigen Values are Known. Find Eigen Vector (Polarity) and Find Spin State of a MODE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc
% clear


%data=open('D:\PhD Research\Corrugated Plate\SH_wave_present_1.mat');
% data=open('D:\PhD Research\Corrugated Plate\Different Coruugation with different D\D_2d\e_0.3_D_2d.mat');
%data=open('spinstateGWUT.mat');

%% Plotting RL Guided wave modes Composed of P-SV 
data_RL=open('wk_PSV_Disp.mat');
data_SH=open('wk_SH_Disp.mat');

h=data_RL.h;
D=data_RL.D;
sol_whole_RL=data_RL.sol_whole;
sol_whole_SH=data_SH.sol_whole;

sol_rl_sz_RL=data_RL.sol_rl_sz;% How many real or Imag sol we have 
sol_rl_sz_SH=data_SH.sol_rl_sz;% How many real or Imag sol we have 

Cp=data_RL.Cp;
Cs=data_RL.Cs;
L_sol_whole=length(sol_whole_RL);
w = sol_whole_RL(:,4); % w has 2*pi factor

kp=w./(Cp);
ks=w./(Cs);

% figure 
% scatter(sol_whole(1:sol_rl_sz(1),3)/1000,sol_whole(1:sol_rl_sz(1),4)./(2*pi*1000),'.')
% set(gca,'FontWeight','bold','FontSize',30, 'LineWidth',2.5); 
% xlabel('Wavenumber [1/mm]','FontSize',30,'FontWeight','bold')
% ylabel('Frequency [KHz]','FontSize',30,'FontWeight','bold')
% title({'Flat Plate'},'FontSize',24,'FontWeight','bold');
% 
% figure 
% scatter(sol_whole(1:sol_rl_sz(1),4)./(2*pi*1000),sol_whole(1:sol_rl_sz(1),3)/1000,'.')
% set(gca,'FontWeight','bold','FontSize',30, 'LineWidth',2.5); 
% xlabel(' Frequency [kHz]','FontSize',30,'FontWeight','bold')
% ylabel('Wavenumber [1/mm]','FontSize',30,'FontWeight','bold')
% title({'corrugated plate with 5 [mm] thickness \epsilon 0.25h'},'FontSize',24,'FontWeight','bold');


% figure 
% scatter(sol_whole(1:sol_rl_sz(1),3)*D,sol_whole(1:sol_rl_sz(1),4)/(2000*pi),'.')
% set(gca,'FontWeight','bold','FontSize',30, 'LineWidth',2.5); 
% hold on
% scatter(-sol_whole(sol_rl_sz(1)+1:L_sol_whole,3)*D,sol_whole(sol_rl_sz(1)+1:L_sol_whole,4)/(2000*pi),'.')
% 
% 
% xlabel('k.D ','FontSize',32,'FontWeight','bold')
% ylabel('Frequncy [KHz]','FontSize',32,'FontWeight','bold')
% title({'D=2d \epsilon=0.1h'},'FontSize',24,'FontWeight','bold');
% % title({'Flat plate'},'FontSize',24,'FontWeight','bold');
% hold off
% % axis([-200 200 0 40])
% box on

% figure
% scatter(sol_whole(1:sol_rl_sz(1),4),(sol_whole(1:sol_rl_sz(1),4))./sol_whole(1:sol_rl_sz(1),3),'r.')
% % hold on
% % scatter(-sol_whole(sol_rl_sz(1)+1:L_sol_whole,4)/(2*pi*1000),sol_whole(sol_rl_sz(1)+1:L_sol_whole,4)/(2*pi*1000)./sol_whole(sol_rl_sz(1)+1:L_sol_whole,3),'.')
% xlabel('Frequency [KHz] ','FontSize',30,'FontWeight','bold')
% ylabel(' Cp [Km/s]','FontSize',30,'FontWeight','bold')
% title({' pahse velocity D=1d \epsilon=0.5 '},'FontSize',24,'FontWeight','bold');
% set(gca,'FontWeight','bold','FontSize',30, 'LineWidth',2.5);
% box on


figure (1)
scatter(sol_whole_RL(1:sol_rl_sz_RL(1),4)/(2*pi)/1000,sol_whole_RL(1:sol_rl_sz_RL(1),3),'b.'); hold on
scatter(sol_whole_SH(1:sol_rl_sz_SH(1),4)/(2*pi)/1000,sol_whole_SH(1:sol_rl_sz_SH(1),3),'r.'); hold on


set(gca,'FontWeight','bold','FontSize',20, 'LineWidth',2.5); 
hold on
%scatter(-sol_whole(sol_rl_sz(1)+1:L_sol_whole,3)*2*h,sol_whole(sol_rl_sz(1)+1:L_sol_whole,4)*2*h/(2*pi*Cs),'.')
plot(w/(2*pi*1000), kp,'k','LineWidth',1.5); hold on % w divided by 2*pi gives the freq in Hz : Plot is in kHz thus divide by 1000
plot(w/(2*pi*1000), ks,'g','LineWidth',1.5); hold on % w divided by 2*pi gives the freq in Hz : Plot is in kHz thus divide by 1000
legend('k','k_h','k_p','k_s')
ylim([0 max(sol_whole_RL(:,3))]);
xlabel('Frequency [kHz]','FontSize',20,'FontWeight','bold')
ylabel('Wavenumber [1/m]','FontSize',20,'FontWeight','bold')
% title({'Obtained from literature [35]'},'FontSize',32,'FontWeight','bold');
title({'Wave Dispersion'},'FontSize',20,'FontWeight','bold');
hold off
box on

%% 

