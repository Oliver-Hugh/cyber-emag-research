% CODE FOR SIMULATION OF WAVES IN INDOOR ENVIRONMENT 
% VERSION: SEPT8, 2023
% ENVIRONMENT: CORNER REFLECTOR (2 WALLS) INCLUDING MULTIPLE SCATTERING

% _______________________

clear 

% SYSTEM

% consider a PEC corner reflector. The Region of Interest (ROI) is 
% within a square of sides L , and we can normalize L = 1 for 
% simplicity. The other dimensions, like wavelength, are then given as 
% a fraction of the room size. Source strength = 1 , no loss of generality.
% We seek to compute the field (denoted s_i) in that environment, 
% for the given transmitter position. 

% NOTE: 
% PHYSICIST VERSUS ENGINEERING NOTATION FOR FOURIER TRANSFORM: we use 
% physicist's notation, so exp(ikr) is wave traveling in positive r sense.

% _________________________

% FIXED VARIABLES 

% wavelength lambda:

lambda=0.2; % notice that the ROI has dimensions 1 times 1, so this 
% lambda is relative to the ROI length (L);

% wavenumber of the field = 2*pi/lambda;

k=2*pi/lambda;

% receiver position: X_r, Y_r; this is variable; 
% transmitter position: X_t, Y_t; we will fix this in this code, per run;

X_t=1.5

Y_t=1.5

% NOTE: Transmitter positions can be outside the ROI. The ROI is 
% for the receiver positions.

% Now we can add the positions of the images: 


%                          image 4     !     * original source 
%                                      !
%                                      !
%------------------------------------------------------------------
%
%                          image 3          image 2

X_t2=X_t;
Y_t2=-Y_t;
X_t3=-X_t;
Y_t3=-Y_t;
X_t4=-X_t;
Y_t4=Y_t;

Reverb=0 % simulates different reflectivities at the wall, between 0-1. 

for ii=1:100
    for jj=1:100

        X_r=ii/100;
        Y_r=jj/100;

        XX(ii,jj)=X_r;
        YY(ii,jj)=Y_r;

% line of sight contribution:
D=sqrt((X_r-X_t)^2+(Y_r-X_t)^2);
s_i1=-exp(i*k*D)/(4*pi*D);

% image 2 contribution: it is 'negative' image;
D=sqrt((X_r-X_t2)^2+(Y_r-Y_t2)^2);
s_i2=exp(i*k*D)/(4*pi*D);

% image 3 contribution: it is 'positive' image;
D=sqrt((X_r-X_t3)^2+(Y_r-Y_t3)^2);
s_i3=-exp(i*k*D)/(4*pi*D);

% image 4 contribution: is it 'negative' image;
D=sqrt((X_r-X_t4)^2+(Y_r-Y_t4)^2);
s_i4=exp(i*k*D)/(4*pi*D);

% total: 

s_i(ii,jj)=s_i1+Reverb*(s_i2+s_i3+s_i4);

%%%% This completes the calculation of the incident field, without
% scatterers in the ROI. This is the incident field only. 

    end
end

% RECOMMENDED PLOTS: 

% Plot mesh(XX,YY,real(s_i)), mesh(XX,YY,imag(s_i)),
% contour(XX,YY,real(s_i),50), etc. for different values of Reverb.
% no Reverb means the incident field is just a clean spherical /circular
% wavefront; if Reverb becomes significant there are interesting ripples.



%% FOR the future... 




%%%%%%%%%%%%%%%%%% adding the scatterers


%% single scatterer

% scattering strength tau1, position X_s1, Y_s1.

tau1=1;
% this is the scattering strength;

% positions: we can use the same positions used in the prior grid, ROI
% but we have to make sure that in the calculations we do not include 
% the same 'ii,jj' when computing the scattered field, since it would 
% give a singularity; we can use the neighbor value in that case; 

II=5;
JJ=5;
X_s1=II/100;
Y_s1=JJ/100;

%X_s1=0.035;

%Y_s1=0.055;

distancex=2*X_s1;
distancey=2*Y_s1;
distancesq=sqrt(distancex^2+distancey^2);

g02x=-exp(i*k*distancex)/(4*pi*distancex);
g02y=-exp(i*k*distancey)/(4*pi*distancey);
g02sq=-exp(i*k*distancesq)/(4*pi*distancesq);

% incident field at the target:

field_incident=s_i(II,JJ);

% total field: 

field_total=field_incident/(1+tau1*g02x+tau1*g02y-tau1*g02sq);

scattered_source=tau1*field_total;

% informal note: looks like a variant is to use ROI for the scatterers,
% but to let the receivers be randomly positioned, some inside ROI and 
% some outside. 

% RECEIVER POSITIONS: 
















