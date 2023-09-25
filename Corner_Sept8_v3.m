% CODE FOR SIMULATION OF WAVES IN INDOOR ENVIRONMENT 
% VERSION: SEPT8, 2023
% ENVIRONMENT: CORNER REFLECTOR (2 WALLS) INCLUDING MULTIPLE SCATTERING

% NOTE: This code focuses on the scattering response matrix for 
% multi-static sensing array outside ROI. 

%%%% Edwin Marengo Fuentes

% _______________________

clear 

% SYSTEM

% consider a PEC corner reflector. The Region of Interest (ROI) is 
% within a square of sides L , and we can normalize L = 1 for 
% simplicity. The other dimensions, like wavelength, are then given as 
% a fraction of the room size. Source strength = 1 , no loss of generality.
% Target is in the ROI. Multi-static array, transmit-receive is outside
% ROI. This code focuses on the computation of the response matrix K.

%% ROI is where the scatterer is located. 

% Schematic: 

%                                       rrrrrr receive array
%
%                               !                        t       transmit
%                               !                         t         array
%                               !                          t
%                               !                           t
%                               !  *                         t
%                               !  target
%_______________________________!__________________________________

% System: we have a transmit array, at an inclination angle, and 
% a receive array, at an inclination angle, located outside ROI 
% where the target is located. Target is a scatterer, it re-emits.
% We have 2 walls, a corner reflector, as the indoor environment in this 
% code. 

% REVERBERATION PARAMETER Reverb: can have Reverb any complex number, 
% to simulate a realistic wall. PEC, perfect electric conductor has 
% Reverb=1, perfect magnetic conductor is Reverb=-1, etc. but we can 
% have cases in between, where magnitude of Reverb <= 1. We can assume
% for simplicity Reverb is a real value between 0 and 1. If magnitude 
% of Reverb is below 1, we have dissipation or loss of energy at the walls
% in the sense that energy is propagated into the interior of the wall,
% and not reflected back into the ROI. 

Reverb=1;

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

x0_transmit=0.5;
y0_transmit=2;
alpha=pi/6; % radians, inclination of the transmit array;
number_transmitters=30;
transmitter_length=2;
transmitter_horizontal_separation=transmitter_length*(cos(alpha));
transmitter_vertical_separation=transmitter_length*(sin(alpha));
x0_receive=1.5;
y0_receive=0.5;
beta=pi/2+pi/6; % inclination of the receive array;
number_receivers=30;
receiver_length=1;
receiver_horizontal_separation=receiver_length*(cos(beta));
receiver_vertical_separation=receiver_length*(sin(beta));

% an epsilon added next, to avoid 'singularities' in the special case 
% of a single transmitter, etc. 

epsilon=1e-6;

for uu=1:number_transmitters % transmitters
    for vv=1:number_receivers % receivers 

        XXtransmit(uu,vv)=x0_transmit+(uu-1)*transmitter_horizontal_separation/(number_transmitters-1+epsilon);
        YYtransmit(uu,vv)=y0_transmit+(uu-1)*transmitter_vertical_separation/(number_transmitters-1+epsilon);
        XXreceive(uu,vv)=x0_receive+(vv-1)*receiver_horizontal_separation/(number_receivers-1+epsilon);
        YYreceive(uu,vv)=y0_receive+(vv-1)*receiver_vertical_separation/(number_receivers-1+epsilon);

X_t=XXtransmit(uu,vv);
Y_t=YYtransmit(uu,vv);
X_r=XXreceive(uu,vv);
Y_r=YYreceive(uu,vv);

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

% in this code we add different degrees of reverberation or reflectivity
% at the walls, via a new parameter Reverb varying between 0 and 1,
% where 0 means no reverberation or reflections, while 1 is the usual 
% PEC conditions. 


%for ii=1:100
%    for jj=1:100

%        X_r=ii/100;
%        Y_r=jj/100;

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

%Reverb=1

s_i(uu,vv)=s_i1+Reverb*(s_i2+s_i3+s_i4); % this is the background signal.

%%%% This completes the calculation of the incident field, without
% scatterers in the ROI. This is the incident field only. 

%%%%%%%%%%%%%%%%%% adding the scatterers.

%% single scatterer

% scattering strength tau1, position X_s1, Y_s1.

tau1=1;
% this is the scattering strength;

% positions: we can use the same positions used in the prior grid, ROI
% but we have to make sure that in the calculations we do not include 
% the same 'ii,jj' when computing the scattered field, since it would 
% give a singularity; we can use the neighbor value in that case; 

%II=ii;
%JJ=jj;
%X_s1=II/100;
%Y_s1=JJ/100;

X_s1=0.2;

Y_s1=0.5;

distancex=2*X_s1;
distancey=2*Y_s1;
distancesq=sqrt(distancex^2+distancey^2);

g02x=-exp(i*k*distancex)/(4*pi*distancex);
g02y=-exp(i*k*distancey)/(4*pi*distancey);
g02sq=-exp(i*k*distancesq)/(4*pi*distancesq);

% incident field at the target:

%%%%%%%%%%%%%%%%% 

X_r=X_s1;
Y_r=Y_s1; % probing field coming from the transmit array; arrives at target;

% line of sight contribution:
D=sqrt((X_r-X_t)^2+(Y_r-X_t)^2);
sTarget_i1=-exp(i*k*D)/(4*pi*D);

% image 2 contribution: it is 'negative' image;
D=sqrt((X_r-X_t2)^2+(Y_r-Y_t2)^2);
sTarget_i2=exp(i*k*D)/(4*pi*D);

% image 3 contribution: it is 'positive' image;
D=sqrt((X_r-X_t3)^2+(Y_r-Y_t3)^2);
sTarget_i3=-exp(i*k*D)/(4*pi*D);

% image 4 contribution: is it 'negative' image;
D=sqrt((X_r-X_t4)^2+(Y_r-Y_t4)^2);
sTarget_i4=exp(i*k*D)/(4*pi*D);

% total: 

%Reverb=1

field_incident=sTarget_i1+Reverb*(sTarget_i2+sTarget_i3+sTarget_i4);


%%%%%%%%%%%%%%%%%

% total field: 

field_total=field_incident/(1+tau1*g02x+tau1*g02y-tau1*g02sq);

scattered_source=tau1*field_total; % this is the amplitude of the 
% scatterer's source, which radiates into the receiver array; 

%%% this radiates, into the receiver array: we evaluate the received field
%%% at the array: 

%%%%%%%%%%%%%%%%%%%%%%%%%%

X_t=X_s1;
Y_t=Y_s1; % target acts as emitter;
X_r=XXreceive(uu,vv);
Y_r=YYreceive(uu,vv); % wave signal sent to the receive array;

X_t2=X_t;
Y_t2=-Y_t;
X_t3=-X_t;
Y_t3=-Y_t;
X_t4=-X_t;
Y_t4=Y_t;

% line of sight contribution:
D=sqrt((X_r-X_t)^2+(Y_r-X_t)^2);
sReceive_i1=-exp(i*k*D)/(4*pi*D);

% image 2 contribution: it is 'negative' image;
D=sqrt((X_r-X_t2)^2+(Y_r-Y_t2)^2);
sReceive_i2=exp(i*k*D)/(4*pi*D);

% image 3 contribution: it is 'positive' image;
D=sqrt((X_r-X_t3)^2+(Y_r-Y_t3)^2);
sReceive_i3=-exp(i*k*D)/(4*pi*D);

% image 4 contribution: is it 'negative' image;
D=sqrt((X_r-X_t4)^2+(Y_r-Y_t4)^2);
sReceive_i4=exp(i*k*D)/(4*pi*D);

% total: 

%Reverb=1

sReceive_i=sReceive_i1+Reverb*(sReceive_i2+sReceive_i3+sReceive_i4);
% this is the Green's function or propagator from target to receivers;

K(uu,vv)=scattered_source*sReceive_i;

% matrix K: this is the 'response matrix' of the scattering system.

%%%%%%%%%%%%%%%%%%%%%%%%%%



    end
end

%%% For the signal processing, we may need the vectorized form of the 
%%% matrices.

s_iv=reshape(s_i,[],1);  % column vector; probing field signals/data.
Kv=reshape(K,[],1); % data perturbation, scattered field signals/data. 
% In the coherent change detection methods we use projections of the two.

% Recommended plots: 

% Geometry: 

plot(XXreceive',YYreceive','o');
hold
plot(XXtransmit,YYtransmit,'+');
plot(X_s1,Y_s1,'*');
axis('square')
axis([0 5 0 5]);
xlabel('x');
ylabel('y');


% plot real, imag parts of s_i, K, subarrays of those, etc. 

% FUTURE additions: add the signal processing for detection.















