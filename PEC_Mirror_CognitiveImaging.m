%% Intoduction
% CODE FOR SIMULATION OF WAVES IN INDOOR ENVIRONMENT
% AUTHORS: ZACHARY SPIEGEL, SUSAN BADANAI
% ENVIRONMENT: PARALLEL PLATE REFLECTOR w/ COGNITIVE IMAGING

% Adapted from "Corner_Sept8_v3.m" and "PEC_target6modified.m" from Edwin A. Marengo,
% "PEC_receiverMatrix_cornerReflector.m" from Jordan Leong, and
% "PEC_Corner_MovingTarget.m"

close all
clear
clc

% Schematic
% TODO - Schematic Here

%% Constants and Initial Setup:

c = 3e8;    % Speed of Light 

tau = 1;    % Target strength     

[receiveCoordinates] = line_array(-10, 3, 201, 20, 0);
[transmitCoordinates] = line_array(-1, 18, 21, 2, 0);

% Target position as a function of time
timeArray = 0:0.01:1;

tgtX = -10:20 / (length(timeArray) - 1):10;
tgtY = sin(-10 + (20/max(timeArray) * timeArray)) + 1.5;

% Define Fixed frequency --> 2.4GHz
w = 2 * pi * 2.4e9;

% Define total number of samples (number of Transmitters * number of Receivers)
totalSamples = length(receiveCoordinates) * length(transmitCoordinates);

% Loop for each time step
for t = 1:length(timeArray)
    % Loop for each receiver
    for a = 1:length(receiveCoordinates)
        % Loop Over Transmitters
        for b = 1:length(transmitCoordinates)
            % Sets coordinates for Transmitters and Receivers
            recX = receiveCoordinates(a, 1);
            recY = receiveCoordinates(a, 2);  
            transX = transmitCoordinates(b, 1);
            transY = transmitCoordinates(b, 2);        
    
            % Find distance from transmitter to receiver (early wave)
            d_trans_rec_early = sqrt((transX - recX)^2 + (transY - recY)^2);
    
            % Find distance from transmitter to receiver (late wave)
            d_trans_rec_late = sqrt((transX - recX)^2 + (transY + recY)^2);
    
            % Find distance from transmitter to target (early wave)
            d_trans_tgt_early = sqrt((transX - tgtX(t))^2 + (transY - tgtY(t))^2);

            % Find distance from transmitter to target (late wave)
            d_trans_tgt_late = sqrt((transX - tgtX(t))^2 + (transY + tgtY(t))^2);

            % Find distance from target to receiver (early wave)
            d_tgt_rec_early = sqrt((tgtX(t) - recX)^2 + (tgtY(t) - recY)^2);

            % Find distance from target to receiver (late wave)
            d_tgt_rec_late = sqrt((tgtX(t) - recX)^2 + (tgtY(t) + recY)^2);
    
            % Calculate respective Green's Function (Impulse Response) for
            % target to receiver
            G_0_early(a, b, t) = greens(d_trans_rec_early, w, c);
            G_0_late(a, b, t) = -1 * greens(d_trans_rec_late, w, c);
    
            G_0(a, b, t) = G_0_late(a, b, t) + G_0_early(a, b, t);

            % Calculate respective Green's Function (Impulse Response) for
            % transmitter to target (incident)
            G_0_inc_early(a, b, t) = greens(d_trans_tgt_early, w, c);
            G_0_inc_late(a, b, t) = -1 * greens(d_trans_tgt_late, w, c);

            G_0_inc(a, b, t) = G_0_inc_early(a, b, t) + G_0_inc_late(a, b, t);

            % Calculate respective Green's Function (Impulse Response) for
            % target to receiver (scatter field)
            G_0_scatter_early(a, b, t) = greens(d_tgt_rec_early, w, c);
            G_0_scatter_late(a, b, t) = -1 * greens(d_tgt_rec_late, w, c);

            G_0_scatter(a, b, t) = G_0_scatter_early(a, b, t) + G_0_scatter_late(a, b, t);

            % Calculate Total Fields
            U_inc(a, b, t) = G_0_inc(a, b, t) / (1 - tau * (1 / (4*pi))) * (exp(1i * w / c * 2 * tgtY(t))) / (2 * tgtY(t));
            U_scatter(a, b, t) = G_0_scatter(a, b, t) * tau * U_inc(a, b, t);

            % Calculate Optical Detector Values
            ot_detector_early(a, b, t) = U_scatter(a, b, t) * conj(G_0_early(a, b, t));
            ot_detector_late(a, b, t) = U_scatter(a, b, t) * conj(G_0_late(a, b, t));
            ot_detector(a, b, t) = U_scatter(a, b, t) * conj(G_0(a, b, t));

            energy_detector(a, b, t) = abs(U_scatter(a, b, t)) ^ 2;
            signal_ot_detector(a, b, t) = abs(ot_detector(a, b, t)) ^ 2;
        end

    % More Code Here
    a = 1;

    % Pre-calculate detector values for each receiver and all transmitters
    ot_detector_early_abs(a, t) = abs(sum(ot_detector_early(a, :, t)));
    ot_detector_late_abs(a, t) = abs(sum(ot_detector_late(a, :, t)));

    ot_detector_abs(a, t) = abs(sum(ot_detector(a, :, t)));
    ot_detector_real(a, t) = real(sum(ot_detector(a, :, t)));
    ot_detector_imag(a, t) = imag(sum(ot_detector(a, :, t)));

    energy_detector_sum(a, t) = sum(energy_detector(a, :, t));

    % Calculat etotal Signal energy
    signal_ot_detector_early_abs(a, t) = abs(ot_detector_early_abs(a, t)) ^ 2;
    signal_ot_detector_late_abs(a, t) = abs(ot_detector_late_abs(a, t)) ^ 2;

    signal_ot_detector_abs(a, t) = abs(ot_detector_abs(a, t)) ^ 2;
    signal_ot_detector_real(a, t) = abs(ot_detector_real(a, t)) ^ 2;
    signal_ot_detector_imag(a, t) = abs(ot_detector_imag(a, t)) ^ 2;

    signal_energy_detector_sum(a, t) = abs(energy_detector_sum(a, t)) ^ 2;

    % Normalized Signal energy
    normalized_signal_ot_detector_early_abs(a, t) = (signal_ot_detector_early_abs(a, t) / norm(G_0_early(a, :, t))) ^ 2;
    normalized_signal_ot_detector_late_abs(a, t) = (signal_ot_detector_late_abs(a, t) / norm(G_0_late(a, :, t))) ^ 2;
    normalized_signal_ot_detector_abs(a, t) = (signal_ot_detector_abs(a, t) / norm(G_0(a, :, t))) ^ 2;

    average_signal_energy(a, t) = signal_energy_detector_sum(a, t) / length(transmitCoordinates);

    % Calculate SNR
    SNR_early(a, t) = normalized_signal_ot_detector_early_abs(a, t) / average_signal_energy(a, t);
    SNR_late(a, t) = normalized_signal_ot_detector_late_abs(a, t) / average_signal_energy(a, t);
    SNR(a, t) = normalized_signal_ot_detector_abs(a, t) / average_signal_energy(a, t);
    end
end

% Line 189
a = 1;


%% Plots

figure
plot(receiveCoordinates(:,1), receiveCoordinates(:,2), 'o');
hold;
plot(transmitCoordinates(:,1), transmitCoordinates(:,2), '+');
plot(tgtX, tgtY, 'v');
plot(-15:0.01:15, 0, 'x', "Color", 'k');
legend('Receivers', 'Transmitters', 'Target Positions', 'PEC');
axis('square')
axis([-15 15 -5 25]);
xlabel('x');
ylabel('y');

title("Simplified Model of Environment Situation");


%% FUNCTIONS

% Function to create linear array of transmitters or receivers
function [XYCoordinates] = line_array(xPos, yPos, numElements, length, angle)
    XXVals = zeros(numElements, 1);
    YYVals = zeros(numElements, 1);
    % an epsilon added next, to avoid 'singularities' in the special case 
    % of a single transmitter, etc. 
    
    epsilon=1e-6;

    horizontalSeparation = length * cos(angle);
    verticalSeparation = length * sin(angle);

    for i = 1:numElements
        XXVals(i) = xPos + (i - 1) * horizontalSeparation / (numElements - 1 + epsilon);
        YYVals(i) = yPos + (i - 1) * verticalSeparation / (numElements - 1 + epsilon);
    end

    XYCoordinates(:,1) = XXVals';
    XYCoordinates(:,2) = YYVals';
end


% Function to calculate Green's Function
function [result] = greens (distance, frequency, c)
    % Formula for Green's Function
    % NOTE 1: 'frequency / c' = k = (2 * pi) / lambda

    result = -1 / (4 * pi) * 1 / distance * exp(1i * (frequency / c) * distance);
end