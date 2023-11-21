clear
% Code to simulate scattering by a point target in a
% waveguide made of 2 parallel plates. Simulation is in 3D, 
% but sources, targets are in the same plane so only 2 spatial 
% variables are used to denote position in the simulations. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% single frequency variant, multiple positions for
%%%%%%%%%% transmit,receive.

%%%%%%%%%%%%%%%%%%%%%%% FOCUS OF THIS CODE APPEARS TO BE THE 
%%%%%%%%%%%%%%%%%%%%%% CALCULATION OF SNR FOR ENERGY VERSUS OT DETECTORS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    t t t t t t t t t



%            X target 



%    r     r    r     r     r 

%%%%%%%%%%%%%%%%%% single frequency 

% Code 1 - uses a single PEC plate ; 

% Foldy-Lax

% Xs = (0,d) position of the receiver; 
% tildeXs=(0,-d) position of the mirror image of the receiver; 

%%%%%%%% multiple positions 

%for kk=1:100

   % d=0.537+kk/10;

d=3; % meters % receiver vertical position

% target strength = tau; position Xt=(xt,yt);

tau=1;

% xt changes in the simulation; NO, xt is -4;

xt=-4;
%yt=0.125/4; 

yt=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y=18; % transmitter vertical position 

% speed of light = c 

c=3e8;

%%%%%%%%%%%%%%%%%%%%%% single frequency 

w=2*pi*2.4e9; %f=2.4 gigahertz;

%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% adding some noise---- for the effect of noise: 

% circularly symmetric: half variance  in real, in imag parts. -

SNR_energy=0.1;

Ntransmitters=20;
Nreceivers=200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nsamples=Nreceivers*Ntransmitters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for pp=1:Nreceivers % receivers ------

    xr=-10+20*pp/Nreceivers;
    Xrec(pp)=xr;

% r = (x,y) ; position of the transmitter; 

for ll=1:Ntransmitters % transmitters -------

    x=-1+2*ll/Ntransmitters; % transmitter positions 
    XX(ll)=x;
%x=20;


%for kk=1:100
%    % frequencies
%    w=(1+kk/25)*3e8*pi;
%    ww(kk)=w;
    dd1=sqrt((x-xr)^2+(y-d)^2);
    dd2=sqrt((x-xr)^2+(y+d)^2);
    G0(pp,ll)=-1/(4*pi)*1/dd1*exp(i*w/c*dd1)+1/(4*pi)*1/dd2*exp(i*w/c*dd2);
    G0_late(pp,ll)=1/(4*pi)*1/dd2*exp(i*w/c*dd2);

    G0_early(pp,ll)=-1/(4*pi)*1/dd1*exp(i*w/c*dd1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dd1h=sqrt((x-xt)^2+(y-yt)^2);
    dd2h=sqrt((x-xt)^2+(y+yt)^2); % for G0_inc;
    dd1out=sqrt((xt-xr)^2+(yt-d)^2);
    dd2out=sqrt((xt-xr)^2+(yt+d)^2); % for GO_out;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    g0(pp,ll)=G0(pp,ll);
    g0aux(pp,ll)=G0(pp,ll);
    
    G0_inc(pp,ll)=-1/(4*pi)*1/dd1h*exp(i*w/c*dd1h)+1/(4*pi)*1/dd2h*exp(i*w/c*dd2h);
    G0_out(pp,ll)=-1/(4*pi)*1/dd1out*exp(i*w/c*dd1out)+1/(4*pi)*1/dd2out*exp(i*w/c*dd2out);




%%%%%%%%%%

%%%% adding target 

% U_inc=G0;

% receiver in the vicinity of the target; transmitter far away at (x,y);

% U_inc - tau U (-1/4pi)(1/(2d)) exp(ik*2d)= U ; 

    U(pp,ll)=G0_inc(pp,ll)/(1-tau*(1/(4*pi))*(exp(i*w/c*2*yt))/(2*yt));
    U_s(pp,ll)=G0_out(pp,ll)*tau*U(pp,ll);
    u_s(pp,ll)=U_s(pp,ll);

    ot_detector_global(pp,ll)=(U_s(pp,ll))*conj(g0aux(pp,ll));
    ot_detector_local(pp,ll)=(U_s(pp,ll))*conj(g0aux(pp,ll));
    energy_detector_global(pp,ll)=(abs(U_s(pp,ll))).^2;

    signal_ot_detector(pp,ll)=(abs(ot_detector_local(pp,ll)))^2;

ot_detector_local_LATE(pp,ll)=(U_s(pp,ll))*conj(G0_late(pp,ll));

ot_detector_local_EARLY(pp,ll)=(U_s(pp,ll))*conj(G0_early(pp,ll));

end

ot_detector_local1(pp)=abs(sum(ot_detector_local(pp,:)));
ot_detector_local2(pp)=real(sum(ot_detector_local(pp,:)));
ot_detector_local3(pp)=imag(sum(ot_detector_local(pp,:)));
energy_detector_local1(pp)=sum(energy_detector_global(pp,:));

ot_detector_local1_LATE(pp)=abs(sum(ot_detector_local_LATE(pp,:)));
ot_detector_local1_EARLY(pp)=abs(sum(ot_detector_local_EARLY(pp,:)));

signal_ot_detector_local1(pp)=(abs(ot_detector_local1(pp)))^2;
signal_ot_detector_local2(pp)=(abs(ot_detector_local2(pp)))^2;
signal_ot_detector_local3(pp)=(abs(ot_detector_local3(pp)))^2;

signal_energy_detector_local1(pp)=(abs(energy_detector_local1(pp)));

signal_ot_detector_local1_LATE(pp)=(abs(ot_detector_local1_LATE(pp)))^2;
signal_ot_detector_local1_EARLY(pp)=(abs(ot_detector_local1_EARLY(pp)))^2;
% sept29,2023 updates: -- based on the theory i add the following: 

normalized_s_ot_det_1(pp)=signal_ot_detector_local1(pp)/(norm(g0aux(pp,:)))^2;

normalized_s_ot_det_1LATE(pp)=signal_ot_detector_local1_LATE(pp)/(norm(G0_late(pp,:)))^2;

normalized_s_ot_det_1EARLY(pp)=signal_ot_detector_local1_EARLY(pp)/(norm(G0_early(pp,:)))^2;

average_s_energy_det_1(pp)=signal_energy_detector_local1(pp)/Ntransmitters;
% average is over the transmitters here; this is the pre-filtering at
% sensor;

ratio_SNRs(pp)=normalized_s_ot_det_1(pp)/average_s_energy_det_1(pp);

ratio_SNRs_LATE(pp)=normalized_s_ot_det_1LATE(pp)/average_s_energy_det_1(pp);

ratio_SNRs_EARLY(pp)=normalized_s_ot_det_1EARLY(pp)/average_s_energy_det_1(pp);

% this ratio is the key QUANTITY HERE. 

end

ot_detector_global1=sum(ot_detector_local1);
ot_detector_global2=sum(ot_detector_local2);
ot_detector_global3=sum(ot_detector_local3);
energy_detector_global1=sum(energy_detector_local1);

signal_ot_detector_global1=(abs(ot_detector_global1))^2;
signal_ot_detector_global2=(abs(ot_detector_global2))^2;
signal_ot_detector_global3=(abs(ot_detector_global3))^2;
signal_energy_detector_global1=(abs(energy_detector_global1));



%%%%%%%%%%%%%%%%%%%%%%%%
% this is noise-free data. let us add noise now. 

var=(norm(U_s))^2/(Nsamples*SNR_energy); % this noise is measured 
% for the global data ; full transmit-receive

for pp=1:Nreceivers
    for ll=1:Ntransmitters
Noise_ot_detector(pp,ll)=var*(abs(g0aux(pp,ll)))^2;
%SNR_ot_detector_local(pp,ll)=signal_ot_detector(pp,ll)/Noise_ot_detector(pp,ll);

    end
    
    Noise_ot_detector1(pp)=sum(Noise_ot_detector(pp,:));


end

Noise_ot_detector_global=sum(Noise_ot_detector1);

SNR_ot_detector_local=signal_ot_detector./Noise_ot_detector;
SNR_ot_1=signal_ot_detector./Noise_ot_detector; % measures SNR for ot det.
% note: identical to prior quantity ;; - note added sept.29/2023
SNR_energy_detector=energy_detector_global./var;

% this is per pp,ll;

SNR_ot_detector_local1=signal_ot_detector_local1./Noise_ot_detector1;
SNR_ot_detector_global1=signal_ot_detector_global1./Noise_ot_detector_global;
SNR_energy_detector1=signal_energy_detector_local1./(Ntransmitters*var);
SNR_energy_detector_global=signal_energy_detector_global1./(Nreceivers*Ntransmitters*var);

% Pd , Pfa, using paper with Gruber, 2013. 

%%%% use sensor at position pp=1; 

SNR=SNR_energy;
NSR=1/SNR;
sigma=sqrt(NSR*signal_energy_detector_local1(160)/Ntransmitters);

for tn=1:10000
% tn = index for normalized threshold , relative to || pi|| , || incident
% field ||

tnorm(tn)=sigma*tn/1000;
refval5=normalized_s_ot_det_1(100)/sigma;

pd(tn)=marcumq(sqrt(2)*sqrt(normalized_s_ot_det_1(100))/sigma,sqrt(2)*tnorm(tn)/sigma,1);
pfa(tn)=exp(-(tnorm(tn))^2/sigma^2);
pfa2(tn)=marcumq(0,sqrt(2)*tnorm(tn)/sigma,1);

pd_EARLY(tn)=marcumq(sqrt(2)*sqrt(normalized_s_ot_det_1EARLY(100))/sigma,sqrt(2)*tnorm(tn)/sigma,1);
pd_LATE(tn)=marcumq(sqrt(2)*sqrt(normalized_s_ot_det_1LATE(100))/sigma,sqrt(2)*tnorm(tn)/sigma,1);


mu=2*signal_energy_detector_local1(160)/sigma^2;
tn2=mu*tn/1000;
pd_energy(tn)=1-ncx2cdf(tn2,2*Ntransmitters,mu);
pfa_energy(tn)=1-ncx2cdf(tn2,2*Ntransmitters,0);
end



M=5
for test=1:1000
    mu=2*M;
   tn3=test/100*mu;
    pd_energy1(test)=1-ncx2cdf(tn3,2*M,mu);
    pfa_energy1(test)=1-ncx2cdf(tn3,2*M,0);
    pfa_energy2(test)=marcumq(0,sqrt(tn3),M);
    pd_energy2(test)=marcumq(sqrt(mu),sqrt(tn3),M);
    pd3(test)=marcumq(M,sqrt(tn3),1);
    pfa3(test)=marcumq(0,sqrt(tn3),1);
    pd4(test)=marcumq(sqrt(M),sqrt(tn3),1);
    pd5(test)=marcumq(1,sqrt(tn3),1);
    pd6(test)=marcumq(sqrt(2*M),sqrt(tn3),1);
end

for uu=1:100
    xground(uu)=-12+24*uu/100;
end

plot(XX,y,'x');
hold
plot(xground,0,'+')
plot(Xrec,d,'*');
plot(xt,yt,'o');
axis([-15 15 -5 25]);