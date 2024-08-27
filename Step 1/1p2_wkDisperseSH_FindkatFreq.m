%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1.2 Program for SH (Shear-Horizontal Modes):  Plug Material Properties and Geometry: Get Dispersion Eigen Values wave number 'k_h' for each 'w' Freq.                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THIS SAMPEL CODE EXPLORES THE INTRINSIC SPIN OF ELASTIC GUIDED WAVE          %%
%% WHEN Eigen Values are Known. Find Eigen Vector (Polarity) and Find Spin State of a MODE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all;
%% Define Geometry of problem

h=1.0e-3;       % average height
d=2*h;
D=1*d;        % Length of period

kmax = pi/D;
freq_start = 1e3;
delfreq = 2e3;
freq_end = 2e6;
Corrugate_Coeff=0;
e=Corrugate_Coeff*h;

%% material properties of problem
E=69e9;
nu=1/3;
rho=2700;

lam=E*nu/((1-2*nu)*(1+nu));
mu=E/(2*(1+nu));
Cp=sqrt((lam+2*mu)/rho);
Cs=sqrt(mu/rho);
S=lam/mu;
T=(4*pi^2*e/D^2);

%% coefficients
Wav_no=0.1:0.1:kmax;
Wav_no=Wav_no';
Freq=freq_start:delfreq:freq_end;
omega=2*pi*Freq;
omega=omega';

DET_COEF=zeros(length(Wav_no),1);
DET=[];
sol=[];
sol_whole=[];
sol_rl_sz=zeros(2,1);

%%
for nn=1:1
for kk=1:1:length(Freq)
        w=omega(kk,1);
        ks=w/Cs;
        kp=w/Cp;
    for j=1:length(Wav_no)
            
        
        k=(1i)^(nn-1)*Wav_no(j,1);
       
        beta=sqrt(ks^2-k^2);
        B1=exp(1i*h*beta);
        B_m1=exp(-1i*h*beta);
        
        C11=beta*besselj(0,e*beta)*B1;
        C12=beta*besselj(0,e*beta)*B_m1;
        C21=beta*besselj(0,e*beta)*B_m1;
        C22=beta*besselj(0,e*beta)*B1;
    
        
        C=[C11 C12; C21 C22;];
        DET_COEF(j,1)=(det(C));

        if j>=2
            Comp_1=imag(DET_COEF(j-1,1));
            REL_1=real(DET_COEF(j-1,1));
            Comp_2=imag(DET_COEF(j,1));
            REL_2=real(DET_COEF(j,1));
            
            if ( Comp_1==0 && Comp_2==0)
                if ( (REL_1*REL_2)<0 ) 
                sol=[sol;Wav_no(j-1),Wav_no(j),(Wav_no(j-1)+Wav_no(j))/2 ,w];
                DET=[DET; DET_COEF(j-1,1), DET_COEF(j,1)];
                end
                
            elseif ( Comp_1~=0 && Comp_2~=0)
                if ( (REL_1*REL_2)<0 &&(Comp_1*Comp_2)<0 ) 
                    sol=[sol;Wav_no(j-1),Wav_no(j), (Wav_no(j-1)+Wav_no(j))/2 ,w];
                    DET=[DET; DET_COEF(j-1,1), DET_COEF(j,1)];
                end
            end
            
        end

    end
    kk
end


DET_COEF=zeros(length(Wav_no),1);
sol_rl_sz(nn,1)=length(sol);
sol_whole=[sol_whole;sol];
sol=[];
end
save('wk_SH_Disp.mat');