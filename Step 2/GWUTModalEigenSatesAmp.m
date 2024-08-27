%% THIS SAMPEL CODE EXPLORES THE INTRINSIC SPIN OF ELASTIC GUIDED WAVE          %%
%% WHEN Eigen Values are Known. Find Eigen Vector (Polarity) and Find Spin State of a MODE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
%% Define Geometry of problem
load ('wk_PSV_Disp.mat','h','d','D','e','E','nu','rho');
%% Provide eigen Value Solution we are Interested in  
k=1348.05;%131;%57;        % selected from k-w solution  It has 2*pi factor in it 
w=1050000*2*pi; 
omg = w/(2*pi);
%125664;%25132.7e5; %selected from k-w solution, %2*pi*Freq;
% Discretization parameters for Spatial Visualization of the Field 
x_min = 0;
x_max = 10*D; % total length along the plate axis 
num_points_x = 1000; % number of points in the x direction
num_points_y = 500; 
%% Where we want to explore the SpinState: 
x1 = 1*D;
x2 = +0;
%% SpinState want to explore
%for sps =1:6
 SpinState = 1;%sps;
%end
% SpinState = 1 ; % s_fi
% SpinState = 2 ; % s_zi
% SpinState = 3 ; % s_fiziud
% SpinState = 4 ; % s_fizidu
% SpinState = 5 ; % s_fiziuu
% SpinState = 6 ; % s_fizidd
%% 1-D SPATIAL DISCREATIZATION 

%% Loaded Variables: Geometry
h;          % average height
d;
D;          % Length of period
e;          %corrugation coefficient
epsilon=e*h;    %corrugation
%% Loaded Variables: Material properties of problem
% E=69e9;
% nu=1/3;
% rho=2700;
%% Dependent Properties 
lam=E*nu/((1-2*nu)*(1+nu));
mu=E/(2*(1+nu));
Cp=sqrt((lam+2*mu)/rho);
Cs=sqrt(mu/rho);
S=lam/mu;
T=(4*pi^2*e/D^2);

%% Kp Ks Wave numbers 
ks=w/Cs;
kp=w/Cp;

%% WAVE PROPAGATION DIRETION OF WAVE 1 and WAVE 2 & WAVEVECTOR
theta_p = real(acosd(k/kp)) ; % wave propagation diraction of P-wave 1
theta_s = real(acosd(k/ks)); % wave propagation diraction of S-wave 2

eta=sqrt(kp^2-k^2);
beta=sqrt(ks^2-k^2);

E1=exp(1i*h*eta);
E_m1=exp(-1i*h*eta);
B1=exp(1i*h*beta);
B_m1=exp(-1i*h*beta);

eta=sqrt((w/Cp)^2-k^2);
beta=sqrt((w/Cs)^2-k^2);
        
C11=((-2*k^2-S*kp^2)*T*(besselj(1,e*eta)/eta)+(2*k^2-ks^2)*besselj(0,e*eta))*E1;
C12=((-2*k^2-S*kp^2)*T*(besselj(1,e*eta)/eta)+(2*k^2-ks^2)*besselj(0,e*eta))*E_m1;
C13=-2*k*beta*(T*(besselj(1,e*beta)/beta)-besselj(0,e*beta))*B1;
C14=2*k*beta*(T*(besselj(1,e*beta)/beta)-besselj(0,e*beta))*B_m1;
      
C21=((-2*k^2-S*kp^2)*T*(besselj(1,e*eta)/eta)+(2*k^2-ks^2)*besselj(0,e*eta))*E_m1;
C22=((-2*k^2-S*kp^2)*T*(besselj(1,e*eta)/eta)+(2*k^2-ks^2)*besselj(0,e*eta))*E1;
C23=-2*k*beta*(T*(besselj(1,e*beta)/beta)-besselj(0,e*beta))*B_m1;
C24=2*k*beta*(T*(besselj(1,e*beta)/beta)-besselj(0,e*beta))*B1;
        
C31=-2*k*eta*(T*(besselj(1,e*eta)/eta)-besselj(0,e*eta))*E1;
C32=2*k*eta*(T*(besselj(1,e*eta)/eta)-besselj(0,e*eta))*E_m1;
C33=(2*k^2-ks^2)*(T*(besselj(1,e*beta)/beta)-besselj(0,e*beta))*B1;
C34=(2*k^2-ks^2)*(T*(besselj(1,e*beta)/beta)-besselj(0,e*beta))*B_m1;
        
C41=-2*k*eta*(T*(besselj(1,e*eta)/eta)-besselj(0,e*eta))*E_m1;
C42=2*k*eta*(T*(besselj(1,e*eta)/eta)-besselj(0,e*eta))*E1;
C43=(2*k^2-ks^2)*(T*(besselj(1,e*beta)/beta)-besselj(0,e*beta))*B_m1;
C44=(2*k^2-ks^2)*(T*(besselj(1,e*beta)/beta)-besselj(0,e*beta))*B1;
        
C=[C11 C12 C13 C14; C21 C22 C23 C24; C31 C32 C33 C34; C41 C42 C43 C44];
clear C11 C12 C13 C14 C22 C23 C24 C31 C32 C33 C34 C41 C42 C43 C44 
clear E1 E_m1 B1 B_m1
display(C)
%% Find Eigen States - Eigen Vectors 
%Use Singular Value Decomposition to find the null space of C
[U, S, V] = svd(C); % this uses the Singular Value DEcomposition which I knew from online. 

%[eigVec, eigVal] = eig(C); % No need as V has a right null space of the C
%matrix 

A = V(:,end); % Solution vector A (non-trivial solution if exists) Given Amplitude / Polarity
clear C 
%% 
% Display the solution
disp('Solution vector A:');
disp(A);
componentNames = {'Apu', 'Apd', 'Asvu', 'Asvd'};
%% Plotting the Field for Verification 

% Create the grid
x_values = linspace(x_min, x_max, num_points_x);
y_values = linspace(-h - epsilon, h + epsilon, num_points_y);

[X, Y] = meshgrid(x_values, y_values);


% Equations for the corrugated surfaces
y_upper_surface = h + epsilon * cos(2 * pi * X / D);
y_lower_surface = -h - epsilon * cos(2 * pi * X / D);

% Determine the points inside the corrugated plate
inside_plate = Y < y_upper_surface & Y > y_lower_surface;

% Determine the points exactly on the upper and lower surface of the corrugated plate
%on_surface_plate = abs(Y - y_upper_surface) <= .0006 | abs(Y - y_lower_surface) <= .01;

% Collect the coordinates of the points inside and on the surface
inside_points = [X(inside_plate), Y(inside_plate)];
%on_surface_points = [X(on_surface_plate), Y(on_surface_plate)];


% % Display the geometry lines and gridpoints - Verified 
% figure;
% hold on;
% plot(x_values, y_upper_surface(1,:), 'k', x_values, y_lower_surface(1,:), 'k'); % Geometry lines
% plot(inside_points(:,1), inside_points(:,2), '.r', 'MarkerSize', 6); % Inside points
% %plot(on_surface_points(:,1), on_surface_points(:,2), '.b', 'MarkerSize', 6); % On-surface points
% axis equal;
% xlabel('X');
% ylabel('Y');
% axis equal; 
% title('Corrugated Plate Geometry with Mesh Points');


%% Calculate φ and ψ for each grid point
Apu = A(1);
Apd = A(2);
Asvu = A(3);
Asvd = A(4);
%% Phi Zi 

phi_pu = @(x, y) Apu * exp(1i * (k * x + eta * y));
phi_pd = @(x, y) Apd * exp(1i * (k * x - eta * y));
psi_SVu = @(x, y) Asvu * exp(1i * (k * x + beta * y));
psi_SVd = @(x, y) Asvd * exp(1i * (k * x - beta * y));

% Displacements u_x and u_y
u_1 = @(x, y) 1i * k * (phi_pu(x, y) + phi_pd(x, y)) + 1i * beta * (psi_SVu(x, y) - psi_SVd(x, y));
u_2 = @(x, y) 1i * eta * (phi_pu(x, y) - phi_pd(x, y)) - 1i * k * (psi_SVu(x, y) + psi_SVd(x, y));

%% Evaluate u_x and u_y for all grid points and plot
U1 = zeros(size(X));
U2 = zeros(size(Y));

for idx = 1:numel(X)
    if inside_plate(idx)
        U1(idx) = u_1(X(idx), Y(idx));
        U2(idx) = u_2(X(idx), Y(idx));
    end
end
% Real parts of Ux and Uy for contour plots, setting values outside the geometry to NaN
Ux_real = real(U1);
Uy_real = real(U2);

Ux_real(~inside_plate) = NaN; % Set values outside the geometry to NaN
Uy_real(~inside_plate) = NaN; % Set values outside the geometry to NaN

Ux_imag = imag(U1);
Uy_imag = imag(U2);

Ux_imag(~inside_plate) = NaN; % Set values outside the geometry to NaN
Uy_imag(~inside_plate) = NaN; % Set values outside the geometry to NaN
%% Plottign the Fields 
X=X*1000;
Y=Y*1000;
%% u1 

% Figure for the contour plot of u_x
figure;
set(gcf, 'Color', 'w'); 
subplot (2,1,1)
contourf(X, Y, Ux_real, 20, 'LineStyle', 'none'); 
colormap("jet"); 
colorbar; 
title('Contour Plot of the Real Part of u_1');
xlabel('x_1 [mm]','FontSize',12.5,'FontWeight','bold');
ylabel('x_2 [mm]','FontSize',12.5,'FontWeight','bold');
axis equal; 
set(gca, 'Color', 'none'); 
set(gca, 'box', 'on'); 
subplot (2,1,2)
contourf(X, Y, Ux_imag, 20, 'LineStyle', 'none'); 
colormap("jet"); 
colorbar; 
title('Contour Plot of the Imag Part of u_1');
xlabel('x_1 [mm]','FontSize',12.5,'FontWeight','bold');
ylabel('x_2 [mm]','FontSize',12.5,'FontWeight','bold');
axis equal; 
set(gca, 'Color', 'none'); 
set(gca, 'box', 'on'); 

%% u2
% Figure for the contour plot of u_x
figure;
set(gcf, 'Color', 'w'); 
subplot (2,1,1)
contourf(X, Y, Uy_real, 20, 'LineStyle', 'none'); 
colormap("jet"); 
colorbar; 
title('Contour Plot of the Real Part of u_2');
xlabel('x_1 [mm]','FontSize',12.5,'FontWeight','bold');
ylabel('x_2 mm]','FontSize',12.5,'FontWeight','bold');
axis equal; 
set(gca, 'Color', 'none'); 
set(gca, 'box', 'on'); 
subplot (2,1,2)
contourf(X, Y, Uy_imag, 20, 'LineStyle', 'none'); 
colormap("jet"); 
colorbar; 
title('Contour Plot of the Imag Part of u_2');
xlabel('x_1 [mm]','FontSize',12.5,'FontWeight','bold');
ylabel('x_2 [mm]','FontSize',12.5,'FontWeight','bold');
axis equal; 
set(gca, 'Color', 'none'); 
set(gca, 'box', 'on'); 

%% Displacement Conjugate 
% Compute the conjugates of Ux and Uy
U1_conj = conj(U1);
U2_conj = conj(U2);

% U_conj = [Ux_conj Uy_conj];
% U = [Ux ; Uy];
% S = [0 -1i ; 1i 0];

% Perform the required multiplication and subtraction
product = U1_conj .* U2 - U2_conj .* U1;
%product_Smat = U_conj.*S.*U;

% Take the imaginary part of the result
imaginary_part = imag(product);
%imaginary_part_Smat = imag(product_Smat);
% Set values outside the geometry to NaN to ensure they are not plotted
%imaginary_part(~inside_plate) = NaN;

%% Spin Field 
%Create the contour plot
figure;
set(gcf, 'Color', 'w'); 
contourf(X, Y, imaginary_part, 20, 'LineStyle', 'none'); 
colormap(jet); 
colorbar; 
title('SAM density');
xlabel('x_1 [mm]','FontSize',12.5,'FontWeight','bold');
ylabel('x_2 [mm]','FontSize',12.5,'FontWeight','bold');
axis equal; 
set(gca, 'Color', 'none');
set(gca, 'box', 'on'); 

% Create the contour plot
% figure;
% set(gcf, 'Color', 'w'); 
% contourf(X, Y, imaginary_part_Smat, 20, 'LineStyle', 'none'); 
% colormap(parula); 
% colorbar; 
% title('Contour Plot of Imag(u_x* . u_y - u_y* . u_x)');
% xlabel('X');
% ylabel('Y');
% axis equal; 
% set(gca, 'Color', 'none');
% set(gca, 'box', 'on');

% Calculate the magnitude of the vector (Ux, Uy)
magnitude = sqrt(real(U1).^2 + real(U2).^2);

%% Real Magnitude Field 
% Create the contour plot
figure;
set(gcf, 'Color', 'w'); 
contourf(X, Y, magnitude, 20, 'LineStyle', 'none'); 
colormap(parula); 
colorbar; 
title('Contour Plot of |sqrt(Ux^2 + Uy^2)|');
xlabel('X [mm]','FontSize',12.5,'FontWeight','bold');
ylabel('Y [mm]', 'FontSize',12.5,'FontWeight','bold');
axis equal; 
set(gca, 'Color', 'none');
set(gca, 'box', 'on'); 

%% Function Call to Get Spin State 
% SpinState = 1 ; % s_fi
% SpinState = 2 ; % s_zi
% SpinState = 3 ; % s_fiziud
% SpinState = 4 ; % s_fizidu
% SpinState = 5 ; % s_fiziuu
% SpinState = 6 ; % s_fizidd
%d=2;
% plotindex=1;
% GWUTSpinFun(k,w, Apu, Apd, Asvu, Asvd, SpinState,...
%                          E, nu, rho, ...
%                          d, x1, x2, x_min, x_max, ...
%                          num_points_x,num_points_y, ...
%                          plotindex)
% 
% %                          (k,w,Au,Ad,Bu,Bd,SpinState, ...
% %                           E, nu, rho,...
% %                           d, xl, xd, x_min, x_max, ...
% %                           num_points_x, num_points_y,...
% %                           plotindex)
%% END 
