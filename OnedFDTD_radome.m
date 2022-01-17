% ONE-DIMENSIONAL FINITE DIFFERENCE TIME DOMAIN ALGORITHM

close all;
clc;
clear all;

c0 = 299792458;
fmax = 2.4*(10^9); % max freq
nmax = 1.86; % highest refractive index
nsrc = 1; % source refractive index
lambda0 = c0/fmax; % free space wavelength
lambda_min = c0/fmax*nmax; % shortest wavelength
NRes = 25;
% NDRes = 10;
% dmin = 1;
nbc = 1; % Boundary condition refractive index

% Compute Grid Resolution
dz1 = lambda_min/NRes; % using minimum lambda and largest refractive index in the grid
% dz2 = dmin/NDRes; % using smallest feature size in the grid
% dz = min(dz1,dz2); 
dz = dz1;

% Snap Grid to Critical Dimensions
dc = lambda0/(4*nmax); %critical dimension (lattice constant/grating period/layer thickness)
N = ceil(dc/dz);
dz = dc/N;
Nz = 180;
% Nz = 200; % Number of grid cells in 1-D

ER = ones(1,Nz);  % free space permittivity
UR = ones(1,Nz);  % free space permeability
nz1=47;
nz2=nz1+37-1;
er_antenna=12;
ur_antenna=1;
er_radome = 3.46;
ur_radome = 1;
ER(nz1:nz2)=er_antenna;
UR(nz1:nz2)=ur_antenna;

nzr1 = 46;
nzr2 = nz2+1;

ER(nzr1:nz1)=er_radome;
ER(nz2:nzr2)=er_radome;
UR(nzr1:nz1)=ur_radome;
UR(nz2:nzr2)=ur_radome;

dt = (nbc*dz)/(2*c0); % wave travels 1 grid cell in 2 time steps

tau = 0.5/fmax;
t0 = 6*tau; % delay for source injection

MEy = (c0*dt)./ER; % update coefficient for Ey
MHx = (c0*dt)./UR; % update coefficient for Hx

tprop = (nmax*Nz*dz)/(c0); % time taken for a wave to propagate across the grid
T = (12*tau)+(5*tprop); % total simulation time 
STEPS = ceil(T/dt); % Number of iterations

Hx = zeros(1,Nz); % Magnetic field vector
Ey = zeros(1,Nz); % Electric field vector


% Source calculations:
nz_src = ceil(Nz/12); % position for source injection
t = [0:STEPS-1]*dt; % time array
delt = ((nsrc*dz)/(2*c0)) + dt/2; % delay between E and H
A = -sqrt(ER(nz_src)/UR(nz_src)); % amplitude factor of H
Esrc = exp(-((t-t0)/tau).^2); % Source Electric field
Hsrc = A*exp(-((t-t0+delt)/tau).^2); % Source Magnetic Field

% Initialize Fourier Transform
NFREQ = 24001;
FREQ = linspace(0, fmax, NFREQ);
K = exp(-1i*2*pi*dt*FREQ); % Kernels for each frequency
REF = zeros(1,NFREQ); % array for Reflectance spectrum
TRN = zeros(1,NFREQ); % array for transmittance spectrum
SRC = zeros(1,NFREQ); % array for Source spectrum

z = linspace(0,Nz*dz,Nz); % 1-D position array

E3=0;E2=0;E1=0;H3=0;H2=0;H1=0;


for T = 1:STEPS


    % Update H from E using Dirichlet Boundary condition
    for nz = 1:Nz-1
        Hx(nz) = Hx(nz)+MHx(nz)*(Ey(nz+1)-Ey(nz))/dz;
    end
    Hx(Nz) = Hx(Nz)+MHx(Nz)*(E3-Ey(Nz))/dz;
    
    Hx(nz_src-1) = Hx(nz_src-1) - MHx(nz_src-1)*Esrc(T)/dz; % H-field source inject
    
    H3 = H2; H2 = H1; H1 = Hx(1); % Record H at Perfect Absorbing Boundary
    
    % Update E from H using Dirichlet Boundary condition
    Ey(1) = Ey(1)+MEy(1)*(Hx(1)-H3)/dz;
    
    for nz = 2:Nz
        Ey(nz) = Ey(nz)+MEy(nz)*(Hx(nz)-Hx(nz-1))/dz;
    end
    
    Ey(nz_src) = Ey(nz_src) - MEy(nz_src)*Hsrc(T)/dz; % E-field source inject
    E3 = E2; E2 = E1; E1 = Ey(Nz); % Record E at Perfect Absorbing Boundary
    
    % Compute Fourier Transforms
    for nf=1:NFREQ
        REF(nf) = REF(nf)+(K(nf)^T)*Ey(1);
        TRN(nf) = TRN(nf)+(K(nf)^T)*Ey(Nz);
        SRC(nf) = SRC(nf)+(K(nf)^T)*Esrc(T);
    end
   
   figure(1)
   set(gcf,'doublebuffer','on');
   if ~mod(T,5)
   x = [nz1*dz nz2*dz nz2*dz nz1*dz nz1*dz];
   y = [-3 -3 3 3 -3];
   fill(x,y,'g');
   hold on
   x = [nzr1*dz nz1*dz nz1*dz nzr1*dz nzr1*dz];
   y = [-3 -3 3 3 -3];
   fill(x,y,'k');
   hold on
   x = [nz2*dz nzr2*dz nzr2*dz nz2*dz nz2*dz];
   y = [-3 -3 3 3 -3];
   fill(x,y,'k');
   hold on
   grid on;
   axis([0 Nz*dz -3 3]);
   plot(z,Ey,'b');
   hold on
   plot(z,Hx,'r');
   hold off
   axis([0 Nz*dz -3 3]);
   title(['Iteration: ',num2str(T)]);
   xlabel('z (m)');
   ylabel('Field Amplitude');
   end
   
end

REF = REF*dt;
TRN = TRN*dt;
SRC = SRC*dt;

REF = abs(REF./SRC).^2;
TRN = abs(TRN./SRC).^2;

figure(2);
plot(FREQ,REF,'-b');
hold on;
plot(FREQ,TRN,'-r');
plot(FREQ,REF+TRN,'-g');
xlim([FREQ(1) FREQ(NFREQ)]);
ylim([0 2]);
xlabel('Frequency (Hz)');
title('REFLECTANCE AND TRANSMITANCE');
hold off;

