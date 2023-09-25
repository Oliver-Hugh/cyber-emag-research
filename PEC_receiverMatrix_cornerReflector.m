clear
% Code to simulate scattering by a point target in a
% waveguide made of 2 parallel plates. Simulation is in 3D,
% but sources, targets are in the same plane so only 2 spatial
% variables are used to denote position in the simulations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Receiver Matrix------- see notes


%%%%%%%%%%%% simulates fixed target position;

%%%% ask target position , tgtX, tgtY ; this is fixed for the simulation ;


%%% meaning of the experiment ---- we have a single transmitter;
%%% we have multiple receiver ---- we fix target;
%%% the output is a plot of each detector averaged over a range of
%%% frequencies

%%%% meaning of the radar experiment: twofold: we can detect locally,
%% at each 'transmitter', in this case via reciprocity they operate as
%% receivers ; the other meaning is they are transmitters and the
%% 'receiver' is the receiver; detection can be post- or pre-processed.
%% post-processed is if the receiver was the transmitter;
%% pre-processed is if the receiver is the receiver and we have used
%% time reversal ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%        |   R   R   Xt target (fixed)   R   R   R   R   R   R   R   R   R   R
%        |
%        P   R   R   R   R   R   R   R   R   R   R   R   R   R   R   R   R   R
%        E
%        C   R   R   R   R   R   R   R   R   R   R   R   R   R   R   R   R   R
%        |
%        |   R   R   R   R   R   R   R   R   R   R   R   R   R   R   R   R   R
% Mirror3|
%        |   R   R   R   R   R   R   R   R   R   R   R   R   R   R   R   R   R
%        |        t transmitter (fixed)
%        |   R   R   R   R   R   R   R   R   R   R   R   R   R   R   R   R   R
%        |
% -------------PEC plate -----------------------------------------
% Mirror2|   Mirror1

% DATA ARE MULTI FREQUENCY, AND FOR MULTIPLE TRANSMIT POSITIONS


%define constants
i = sqrt(-1);
c = 3e8;


% R = (recX, recY) positions of the receivers
num_receivers = 101;
recX_1pos = 1;
recX_lastpos = 21;
recX = linspace(recX_1pos, recX_lastpos, num_receivers);

recY_1pos = 1;
recY_lastpos = 21;
recY = linspace(recY_1pos, recY_lastpos, num_receivers);


% t = (transX, transY) position of the transmitter
transX = 3;
transY = 0.7;


% Xt = (tgtX, tgtY) position of the target
tgtX=6;
tgtY=5;

% target strength = tau
tau=0.2;

% frequency range
num_frequencies = 100;

% Creates array of angular frequencies (rad/s) ranging from 1.04c(pi) to 5c(pi)
w = zeros(1,num_frequencies);
for index = 1:length(w)
    w(index) = (1 + index/25) * c * pi;
end

% representative frequencies
% typical distance is the distance between PEC and target
dist = (tgtY+tgtX)/2;
% distance = 0.2 wavelength
wSpec(1) = 2 * pi * c * 0.2 / dist;
% distance = 0.5 wavelength
wSpec(2) = 2 * pi * c * 0.5 / dist;
% distance = 1 wavelength
wSpec(3) = 2 * pi * c / dist;
% distance = 2 wavelength
wSpec(4) = 2 * pi * c * 2 / dist;
% distance = 5 wavelength
wSpec(5) = 2 * pi * c * 5 / dist;

% X slice value
Xslice = tgtX;
[minDifference, Xslice_index] = min(abs(recX - Xslice));

% Y slice value
Yslice = tgtY;
[minDifference, Yslice_index] = min(abs(recY - Yslice));

% Create empty matrices to hold values
U_bknd = zeros(length(recX), length(recY), length(w));
U_out = zeros(length(recX), length(recY), length(w));
U_s = zeros(length(recX), length(recY), length(w));

U_inc = zeros(length(w));
U_tgt = zeros(length(w));

ot_detector = zeros(length(recX), length(recY), length(w));
ot_detector_avg = zeros(length(recX), length(recY));
ot_detector_Xslice = zeros(size(recY));
ot_detector_Yslice = zeros(size(recX));

energy_detector = zeros(length(recX), length(recY), length(w));
energy_detector_avg = zeros(length(recX), length(recY));
energy_detector_Xslice = zeros(size(recY));
energy_detector_Yslice = zeros(size(recX));

% for specific frequencies
U_bkndSpec = zeros(length(recX), length(recY), length(wSpec));
U_outSpec = zeros(length(recX), length(recY), length(wSpec));
U_sSpec = zeros(length(recX), length(recY), length(wSpec));

U_incSpec = zeros(length(wSpec));
U_tgtSpec = zeros(length(wSpec));

ot_detectorSpec = zeros(length(recX), length(recY), length(wSpec));
energy_detectorSpec = zeros(length(recX), length(recY), length(wSpec));

ot_detector_XsliceSpec = zeros(length(recY), length(wSpec));
ot_detector_YsliceSpec = zeros(length(recX), length(wSpec));
energy_detector_XsliceSpec = zeros(length(recY), length(wSpec));
energy_detector_YsliceSpec = zeros(length(recX), length(wSpec));


% Formula Calculations

greens = @(distance, freq) -1/(4*pi)*1/distance*exp(i*freq/c*distance);
tgt_ref = @(distance, freq) -(tau*(1/(4*pi))*(exp(i*freq/c*distance))/(distance));



% Traverse receiver matrix
for XX = 1:length(recX)
    for YY = 1:length(recY)
        
        % Distance between transmitter and receiver
        d_trans_rec = sqrt((transX-recX(XX))^2 + (transY-recY(YY))^2);
        % Distance between transmitter and receiver mirror images
        d_trans_recMirror(1) = sqrt((transX-recX(XX))^2 + (transY+recY(YY))^2);
        d_trans_recMirror(2) = sqrt((transX+recX(XX))^2 + (transY+recY(YY))^2);
        d_trans_recMirror(3) = sqrt((transX+recX(XX))^2 + (transY-recY(YY))^2);
        
        % Distance between transmitter and target
        d_trans_tgt = sqrt((transX-tgtX)^2 + (transY-tgtY)^2);
        % Distance between transmitter and target mirror images
        d_trans_tgtMirror(1) = sqrt((transX-tgtX)^2 + (transY+tgtY)^2);
        d_trans_tgtMirror(2) = sqrt((transX+tgtX)^2 + (transY+tgtY)^2);
        d_trans_tgtMirror(3) = sqrt((transX+tgtX)^2 + (transY-tgtY)^2);
    
        % Distance between target and receiver
        d_tgt_rec = sqrt((tgtX-recX(XX))^2 + (tgtY-recY(YY))^2);
        % Distance between target and receiver mirror images
        d_tgt_recMirror(1) = sqrt((tgtX-recX(XX))^2 + (tgtY+recY(XX))^2);
        d_tgt_recMirror(2) = sqrt((tgtX+recX(XX))^2 + (tgtY+recY(XX))^2);
        d_tgt_recMirror(3) = sqrt((tgtX+recX(XX))^2 + (tgtY-recY(XX))^2);
        
        %target self reflection distances
        d_tgt_mirror(1) = 2*tgtY;
        d_tgt_mirror(2) = 2*sqrt(tgtX^2 + tgtY^2);
        d_tgt_mirror(3) = 2*tgtX;
    
    
        % Traverse frequencies
        for WW = 1:length(w)
            
            % Background field transmitter to receiver
            U_bknd_comps(1) = greens(d_trans_rec, w(WW));
            U_bknd_comps(2) = -1*greens(d_trans_recMirror(1), w(WW));
            U_bknd_comps(3) = greens(d_trans_recMirror(2), w(WW));
            U_bknd_comps(4) = -1*greens(d_trans_recMirror(3), w(WW));
            U_bknd(XX,YY,WW) = sum(U_bknd_comps);
            
            % Incident field transmitter to target
            U_inc_comps(1) = greens(d_trans_tgt, w(WW));
            U_inc_comps(2) = -1*greens(d_trans_tgtMirror(1), w(WW));
            U_inc_comps(3) = greens(d_trans_tgtMirror(2), w(WW));
            U_inc_comps(4) = -1*greens(d_trans_tgtMirror(3), w(WW));
            U_inc(WW) = sum(U_inc_comps);
            
            
            % Target reflections
            tgt_ref_comps(1) = tgt_ref(d_tgt_mirror(1), w(WW));
            tgt_ref_comps(2) = -1 * tgt_ref(d_tgt_mirror(2), w(WW));
            tgt_ref_comps(3) = tgt_ref(d_tgt_mirror(3), w(WW));
            
            % Total field at target
            U_tgt(WW)=U_inc(WW)/(1 + sum(tgt_ref_comps));
            
            
            % Field at receiver not influenced by target field or tau
            U_out_comps(1) = greens(d_tgt_rec, w(WW));
            U_out_comps(2) = -1 * greens(d_tgt_recMirror(1), w(WW));
            U_out_comps(3) = greens(d_tgt_recMirror(2), w(WW));
            U_out_comps(4) = -1 * greens(d_tgt_recMirror(3), w(WW));
            U_out(XX,YY,WW) = sum(U_out_comps);
       
            % Field from target to receiver
            U_s(XX,YY,WW)=U_out(XX,YY,WW)*tau*U_tgt(WW);
          
            % Optical Theorem Detector value
            ot_detector(XX,YY,WW) = U_s(XX,YY,WW) * conj(U_bknd(XX,YY,WW));
            
            % Energy Detector value
            energy_detector(XX,YY,WW) = abs(U_s(XX,YY,WW)).^2;
        end
       
        % Average Optical Theorem Detector values over all frequencies
        ot_detector_avg(XX,YY)= sum(ot_detector(XX,YY,:))/length(w);
        
        % Average Energy Detector values over all frequencies
        energy_detector_avg(XX,YY) = sum(energy_detector(XX,YY,:))/length(w);
        
        % Average Optical Theorem Detector values over all frequencies
        ot_detector_Xslice(YY)= ot_detector_avg(Xslice_index,YY);
        
        % Average Energy Detector values over all frequencies
        energy_detector_Xslice(YY) = sum(energy_detector(Xslice_index,YY,:))/length(w);
        
        % Average Optical Theorem Detector values over all frequencies
        ot_detector_Yslice(XX)= ot_detector_avg(XX,Yslice_index);
        
        % Average Energy Detector values over all frequencies
        energy_detector_Yslice(XX) = sum(energy_detector(XX,Yslice_index,:))/length(w);
        
         % Traverse specific frequencies
        for WW = 1:length(wSpec)
            % Background field transmitter to receiver
            U_bknd_comps(1) = greens(d_trans_rec, w(WW));
            U_bknd_comps(2) = -1*greens(d_trans_recMirror(1), w(WW));
            U_bknd_comps(3) = greens(d_trans_recMirror(2), w(WW));
            U_bknd_comps(4) = -1*greens(d_trans_recMirror(3), w(WW));
            U_bkndSpec(XX,YY,WW) = sum(U_bknd_comps);
            
            % Incident field transmitter to target
            U_inc_comps(1) = greens(d_trans_tgt, w(WW));
            U_inc_comps(2) = -1*greens(d_trans_tgtMirror(1), w(WW));
            U_inc_comps(3) = greens(d_trans_tgtMirror(2), w(WW));
            U_inc_comps(4) = -1*greens(d_trans_tgtMirror(3), w(WW));
            U_incSpec(WW) = sum(U_inc_comps);
            
            
            % Target reflections
            tgt_ref_comps(1) = tgt_ref(d_tgt_mirror(1), w(WW));
            tgt_ref_comps(2) = -1 * tgt_ref(d_tgt_mirror(2), w(WW));
            tgt_ref_comps(3) = tgt_ref(d_tgt_mirror(3), w(WW));
            
            % Total field at target
            U_tgtSpec(WW)=U_incSpec(WW)/(1 + sum(tgt_ref_comps));
            
            
            % Field at receiver not influenced by target field or tau
            U_out_comps(1) = greens(d_tgt_rec, w(WW));
            U_out_comps(2) = -1 * greens(d_tgt_recMirror(1), w(WW));
            U_out_comps(3) = greens(d_tgt_recMirror(2), w(WW));
            U_out_comps(4) = -1 * greens(d_tgt_recMirror(3), w(WW));
            U_outSpec(XX,YY,WW) = sum(U_out_comps);
       
            % Field from target to receiver
            U_sSpec(XX,YY,WW)=U_outSpec(XX,YY,WW)*tau*U_tgtSpec(WW);
          
            % Optical Theorem Detector value
            ot_detectorSpec(XX,YY,WW) = U_sSpec(XX,YY,WW) * conj(U_bkndSpec(XX,YY,WW));
            
            % Energy Detector value
            energy_detectorSpec(XX,YY,WW) = abs(U_sSpec(XX,YY,WW)).^2;
        
            
            % Optical Theorem Detector X-slice
            ot_detector_XsliceSpec(YY,WW)= ot_detectorSpec(Xslice_index,YY,WW);
            
            % Energy Detector X-slice
            energy_detector_XsliceSpec(YY,WW) = energy_detectorSpec(Xslice_index,YY,WW);
            
            % Optical Theorem Detector Y-slice
            ot_detector_YsliceSpec(XX,WW)= ot_detectorSpec(XX,Yslice_index,WW);
            
            % Energy Detector Y-slice
            energy_detector_YsliceSpec(XX,WW) = energy_detectorSpec(XX,Yslice_index,WW);
        end
        
    end
end


% 3D mesh
figure('Position',[100 100 1500 1000]);
subplot(2,3,1)
surf(recX, recY, transpose(abs(ot_detector_avg)));
title('Optical Theorem Detector: f_{ot,m} values averaged over frequency');
xlabel('Receiver X position');
ylabel('Receiver Y position');
zlabel('SNR');

subplot(2,3,2)
surf(recX, recY, transpose(real(ot_detector_avg)));
title("Optical Theorem Detector: f_{ot,r} values averaged over frequency");
xlabel('Receiver X position');
ylabel('Receiver Y position');
zlabel('SNR');


subplot(2,3,3)
surf(recX, recY, transpose(imag(ot_detector_avg)));
title("Optical Theorem Detector: f_{ot,i} values averaged over frequency");
xlabel('Receiver X position');
ylabel('Receiver Y position');
zlabel('SNR');


subplot(2,3,4)
surf(recX, recY, transpose((energy_detector_avg)));
title("Energy Detector: f_e values averaged over frequency");
xlabel('Receiver X position');
ylabel('Receiver Y position');
zlabel('SNR');

% 2D X slice
% figure('Position',[100 100 1500 1000]);
% subplot(2,3,1)
% plot(recY, abs(ot_detector_Xslice));
% title(["Optical Theorem Detector: f_{ot,m} values averaged over frequency at value X=",num2str(Xslice)]);
% xlabel("Receiver Y position");
% ylabel('SNR');
% 
% subplot(2,3,2)
% plot(recY, real(ot_detector_Xslice));
% title(["Optical Theorem Detector: f_{ot,r} values averaged over frequency at value X=",num2str(Xslice)]);
% xlabel("Receiver Y position");
% ylabel('SNR');
% 
% subplot(2,3,3)
% plot(recY, imag(ot_detector_Xslice));
% title(["Optical Theorem Detector: f_{ot,i} values averaged over frequency at value X=",num2str(Xslice)]);
% xlabel("Receiver Y position");
% ylabel('SNR');
% 
% subplot(2,3,4)
% plot(recY, (energy_detector_Xslice));
% title(["Energy Detector: f_e values averaged over frequency at value X=",num2str(Xslice)]);
% xlabel("Receiver Y position");
% ylabel('SNR');

% 2D Y slice
% figure('Position',[100 100 1500 1000]);
% subplot(2,3,1)
% plot(recX, abs(ot_detector_Yslice));
% title(["Optical Theorem Detector: f_{ot,m} values averaged over frequency at value Y=",num2str(Yslice)]);
% xlabel("Receiver X position");
% ylabel('SNR');
% 
% subplot(2,3,2)
% plot(recX, real(ot_detector_Yslice));
% title(["Optical Theorem Detector: f_{ot,r} values averaged over frequency at value Y=",num2str(Yslice)]);
% xlabel("Receiver X position");
% ylabel('SNR');
% 
% subplot(2,3,3)
% plot(recX, imag(ot_detector_Yslice));
% title(["Optical Theorem Detector: f_{ot,i} values averaged over frequency at value Y=",num2str(Yslice)]);
% xlabel("Receiver X position");
% ylabel('SNR');
% 
% subplot(2,3,4)
% plot(recY, (energy_detector_Yslice));
% title(["Energy Detector: f_e values averaged over frequency at value X=",num2str(Yslice)]);
% xlabel("Receiver X position");
% ylabel('SNR');

% 3D mesh specific frequencies
figure('Position',[100 100 2000 1000]);
subplot(2,5,1)
surf(recX, recY, transpose(abs(ot_detectorSpec(:,:,1))));
title('Optical Theorem Detector: f_{ot,m} values at distance = 0.2\lambda');
xlabel('Receiver X position');
ylabel('Receiver Y position');
zlabel('SNR');

subplot(2,5,6)
surf(recX, recY, transpose(abs(energy_detectorSpec(:,:,1))));
title('Energy Detector: f_{e} values at distance = 0.2\lambda');
xlabel('Receiver X position');
ylabel('Receiver Y position');
zlabel('SNR');

subplot(2,5,2)
surf(recX, recY, transpose(abs(ot_detectorSpec(:,:,2))));
title('Optical Theorem Detector: f_{ot,m} values at distance = 0.5\lambda');
xlabel('Receiver X position');
ylabel('Receiver Y position');
zlabel('SNR');

subplot(2,5,7)
surf(recX, recY, transpose(abs(energy_detectorSpec(:,:,2))));
title('Energy Detector: f_{e} values at distance = 0.5\lambda');
xlabel('Receiver X position');
ylabel('Receiver Y position');
zlabel('SNR');

subplot(2,5,3)
surf(recX, recY, transpose(abs(ot_detectorSpec(:,:,3))));
title('Optical Theorem Detector: f_{ot,m} values at distance = 1\lambda');
xlabel('Receiver X position');
ylabel('Receiver Y position');
zlabel('SNR');

subplot(2,5,8)
surf(recX, recY, transpose(abs(energy_detectorSpec(:,:,3))));
title('Energy Detector: f_{e} values at distance = 1\lambda');
xlabel('Receiver X position');
ylabel('Receiver Y position');
zlabel('SNR');

subplot(2,5,4)
surf(recX, recY, transpose(abs(ot_detectorSpec(:,:,4))));
title('Optical Theorem Detector: f_{ot,m} values at distance = 2\lambda');
xlabel('Receiver X position');
ylabel('Receiver Y position');
zlabel('SNR');

subplot(2,5,9)
surf(recX, recY, transpose(abs(energy_detectorSpec(:,:,4))));
title('Energy Detector: f_{e} values at distance = 2\lambda');
xlabel('Receiver X position');
ylabel('Receiver Y position');
zlabel('SNR');

subplot(2,5,5)
surf(recX, recY, transpose(abs(ot_detectorSpec(:,:,5))));
title('Optical Theorem Detector: f_{ot,m} values at distance = 5\lambda');
xlabel('Receiver X position');
ylabel('Receiver Y position');
zlabel('SNR');

subplot(2,5,10)
surf(recX, recY, transpose(abs(energy_detectorSpec(:,:,5))));
title('Energy Detector: f_{e} values at distance = 5\lambda');
xlabel('Receiver X position');
ylabel('Receiver Y position');
zlabel('SNR');

% X slice specific frequencies
figure('Position',[100 100 2000 1000]);
subplot(2,5,1)
plot(recY, abs(ot_detector_XsliceSpec(:,1)));
title(["Optical Theorem Detector: f_{ot,m} values at distance = 0.2\lambda at value X=",num2str(Xslice)]);
xlabel("Receiver Y position");
ylabel('SNR');

subplot(2,5,6)
plot(recY, abs(energy_detector_XsliceSpec(:,1)));
title(["Energy Detector: f_{e} values at distance = 0.2\lambda at value X=",num2str(Xslice)]);
xlabel("Receiver Y position");
ylabel('SNR');

subplot(2,5,2)
plot(recY, abs(ot_detector_XsliceSpec(:,2)));
title(["Optical Theorem Detector: f_{ot,m} values at distance = 0.5\lambda at value X=",num2str(Xslice)]);
xlabel("Receiver Y position");
ylabel('SNR');

subplot(2,5,7)
plot(recY, abs(energy_detector_XsliceSpec(:,2)));
title(["Energy Detector: f_{e} values at distance = 0.5\lambda at value X=",num2str(Xslice)]);
xlabel("Receiver Y position");
ylabel('SNR');

subplot(2,5,3)
plot(recY, abs(ot_detector_XsliceSpec(:,3)));
title(["Optical Theorem Detector: f_{ot,m} values at distance = 1\lambda at value X=",num2str(Xslice)]);
xlabel("Receiver Y position");
ylabel('SNR');

subplot(2,5,8)
plot(recY, abs(energy_detector_XsliceSpec(:,3)));
title(["Energy Detector: f_{e} values at distance = 1\lambda at value X=",num2str(Xslice)]);
xlabel("Receiver Y position");
ylabel('SNR');

subplot(2,5,4)
plot(recY, abs(ot_detector_XsliceSpec(:,4)));
title(["Optical Theorem Detector: f_{ot,m} values at distance = 2\lambda at value X=",num2str(Xslice)]);
xlabel("Receiver Y position");
ylabel('SNR');

subplot(2,5,9)
plot(recY, abs(energy_detector_XsliceSpec(:,4)));
title(["Energy Detector: f_{e} values at distance = 2\lambda at value X=",num2str(Xslice)]);
xlabel("Receiver Y position");
ylabel('SNR');

subplot(2,5,5)
plot(recY, abs(ot_detector_XsliceSpec(:,5)));
title(["Optical Theorem Detector: f_{ot,m} values at distance = 5\lambda at value X=",num2str(Xslice)]);
xlabel("Receiver Y position");
ylabel('SNR');

subplot(2,5,10)
plot(recY, abs(energy_detector_XsliceSpec(:,5)));
title(["Energy Detector: f_{e} values at distance = 5\lambda at value X=",num2str(Xslice)]);
xlabel("Receiver Y position");
ylabel('SNR');


[maxValue, linearIndex] = max(ot_detector_avg(:));
[rowIndex, columnIndex] = ind2sub(size(ot_detector_avg), linearIndex);

disp(['Optical Theorem Predicted Target location: (', num2str(recX(rowIndex)), ', ', num2str(recY(columnIndex)), ')']);

[maxValue, linearIndex] = max(energy_detector_avg(:));
[rowIndex, columnIndex] = ind2sub(size(energy_detector_avg), linearIndex);

disp(['Energy Detector Predicted Target location: (', num2str(recX(rowIndex)), ', ', num2str(recY(columnIndex)), ')']);



