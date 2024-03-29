%% Intoduction
% CODE FOR SIMULATION OF WAVES IN INDOOR ENVIRONMENT 
% AUTHORS: ZACHARY SPIEGEL, MATTHEW LOW, OLIVER HUGH, MICHAEL SEKENSKI
% ENVIRONMENT: CORNER REFLECTOR (2 WALLS) INCLUDING MULTIPLE SCATTERING

% Adapted from "Corner_Sept8_v3.m" from Edwin A. Marengo and
% "PEC_receiverMatrix_cornerReflector.m" from Jordan Leong

close all
clear
clc

% Schematic
% * = target
% R = receiver
% T = transmitter
%        |
%        |       T [TRANSMIT ARRAY]
%        |      T
%        P     T
%        E    T
%        C   T
%        | 
%        |              R
% Mirror3|               R
%        |                R     [RECEIVE ARRAY]
%        |  *              R
%        |*                 R
%        |     *
%        |  [MULTIPLE TARGET POSITIONS AS A FUNCTION OF TIME]
%        |         *
%        |
% -------------PEC plate -----------------------------------------
% Mirror2|   Mirror1

% Conditions to account for:
%   1) Field between any given transmitter and receiver
%   2) Field between any given transmitter and target
%   3) Field between and given receiver and target



%% Constants and Initial Setup:

c = 3e8;    % Speed of Light 


% REVERBERATION PARAMETER Reverb: can have Reverb any complex number, 
% to simulate a realistic wall. PEC, perfect electric conductor has 
% Reverb=1, perfect magnetic conductor is Reverb=-1. If magnitude 
% of Reverb is below 1, we have dissipation or loss of energy at the walls
% in the sense that energy is propagated into the interior of the wall,
% and not reflected back into the ROI.

Reverb = 1;

% Creates XY coordinate pairs for each transmitter and receiver in the
% situation. Definitions for "line_array" function are:
%   1) Starting X-Position
%   2) Starting Y-Position
%   3) Number of array elements
%   4) Array Length (meters)
%   5) Angle of array (radians)
[transmitCoordinates] = line_array(0.5, 2, 30, 2, pi/6);
[receiveCoordinates] = line_array(1.5, 0.5, 30, 1, (pi/2 + pi/6));

% Sets the time scale for the movement of the target
timeArray = 0:0.001:10;

% tgtX, tgtY: they indicate the targets' positions. 
% Also, we have the target strengths tau. 

tgtX = zeros(1, length(timeArray));
tgtY = zeros(1, length(timeArray));

tgtX = (0.1 .* timeArray) + 0.1;
tgtY = (0.01 .* timeArray .^ 2) + 0.1;
tau = 1;

% Sets Range of Frequencies
% num_frequencies = 225;

% Creates array of angular frequencies (rad/s) ranging from 1.04c(pi) [lambda ~= 1.92] to 10c(pi) [lambda = 0.2]
% w = zeros(1,num_frequencies);
% for index = 1:length(w)
%     w(index) = (1 + index/25) * c * pi;
% end

% Setting omega = 10*c*pi (lambda = 0.2)
w = 10*c*pi;

% Preallocate matrix sizes
U_bknd = zeros(length(transmitCoordinates), length(receiveCoordinates), length(timeArray));
U_inc = zeros(length(transmitCoordinates), length(receiveCoordinates), length(timeArray));
U_inc_total = zeros(length(transmitCoordinates), length(receiveCoordinates), length(timeArray));
U_scatt = zeros(length(transmitCoordinates), length(receiveCoordinates), length(timeArray));
U_out = zeros(length(transmitCoordinates), length(receiveCoordinates), length(timeArray));
K = zeros(length(transmitCoordinates), length(receiveCoordinates), length(timeArray));



%% CALCULATE TOTAL FIELD
for UU = 1:length(transmitCoordinates)
    for VV = 1:length(receiveCoordinates)
        % Sets coordinates for Transmitters and Receivers
        transX = transmitCoordinates(UU, 1);
        transY = transmitCoordinates(UU, 2);
        recX = receiveCoordinates(VV, 1);
        recY = receiveCoordinates(VV, 2);

        % Distance between transmitter and receiver
        d_trans_rec = sqrt((recX - transX)^2 + (recY - transY)^2);
        % Distance between transmitter and receiver mirror images
        d_trans_recMirror(1) = sqrt((recX - transX)^2 + (recY + transY)^2);
        d_trans_recMirror(2) = sqrt((recX + transX)^2 + (recY + transY)^2);
        d_trans_recMirror(3) = sqrt((recX + transX)^2 + (recY - transY)^2);

        % Background field transmitter to receiver
        U_bknd_comps(1) = greens(d_trans_rec, w, c);
        U_bknd_comps(2) = -1 * greens(d_trans_recMirror(1), w, c);
        U_bknd_comps(3) = greens(d_trans_recMirror(2), w, c);
        U_bknd_comps(4) = -1 * greens(d_trans_recMirror(3), w, c);
        U_bknd(UU, VV) = U_bknd_comps(1) + Reverb*(U_bknd_comps(2) + U_bknd_comps(3) + U_bknd_comps(4));    
        
        for TT = 1:length(timeArray)
            % Calculates Green's Function for target to determine total field due to scatterer
            d_tgtX = 2 * tgtX(TT);
            d_tgtY = 2 * tgtY(TT);
            d_tgtSQRT = sqrt(d_tgtX^2 + d_tgtY^2);

            U_tgtX = greens(d_tgtX, w, c);
            U_tgtY = greens(d_tgtY, w, c);
            U_tgtSQRT = greens(d_tgtSQRT, w, c);

            % Distance between transmitter and target
            d_trans_tgt = sqrt((tgtX(TT) - transX)^2 + (tgtY(TT) - transY)^2);
            % Distance between transmitter and target mirror images
            d_trans_tgtMirror(1) = sqrt((tgtX(TT) - transX)^2 + (tgtY(TT) + transY)^2);
            d_trans_tgtMirror(2) = sqrt((tgtX(TT) + transX)^2 + (tgtY(TT) + transY)^2);
            d_trans_tgtMirror(3) = sqrt((tgtX(TT) + transX)^2 + (tgtY(TT) - transY)^2); 

            % Distance between target and receiver
            d_tgt_rec = sqrt((recX - tgtX(TT))^2 + (recY - tgtY(TT))^2);
            % Distance between target and receiver mirror images (target acts as emitter)
            d_tgt_recMirror(1) = sqrt((recX - tgtX(TT))^2 + (recY + tgtY(TT))^2);
            d_tgt_recMirror(2) = sqrt((recX + tgtX(TT))^2 + (recY + tgtY(TT))^2);
            d_tgt_recMirror(3) = sqrt((recX + tgtX(TT))^2 + (recY - tgtY(TT))^2);    
            % Incident field transmitter to target, with calculation of
            % total scatterer field
            U_inc_comps(1) = greens(d_trans_tgt, w, c);
            U_inc_comps(2) = -1 * greens(d_trans_tgtMirror(1), w, c);
            U_inc_comps(3) = greens(d_trans_tgtMirror(2), w, c);
            U_inc_comps(4) = -1 * greens(d_trans_tgtMirror(3), w, c);
            U_inc(UU, VV, TT) = U_inc_comps(1) + Reverb*(U_inc_comps(2) + U_inc_comps(3) + U_inc_comps(4));  
            U_inc_total(UU, VV, TT) = U_inc(UU, VV, TT) / (1 + tau*U_tgtX + tau*U_tgtY - tau*U_tgtSQRT);
            U_scatt(UU, VV, TT) = tau * U_inc_total(UU, VV, TT);
         
            % Field at receiver not influenced by target or tau
            U_out_comps(1) = greens(d_tgt_rec, w, c);
            U_out_comps(2) = -1 * greens(d_tgt_recMirror(1), w, c);
            U_out_comps(3) = greens(d_tgt_recMirror(2), w, c);
            U_out_comps(4) = -1 * greens(d_tgt_recMirror(3), w, c);
            U_out(UU, VV, TT) = U_out_comps(1) + Reverb*(U_out_comps(2) + U_out_comps(3) + U_out_comps(4));

            % Total Field K
            K(UU, VV, TT) = U_out(UU, VV, TT) * U_scatt(UU, VV, TT);
        end
    end
end



%% PLOTS

% Plot for schematic of environment
% figure
% plot(receiveCoordinates(:,1), receiveCoordinates(:,2), 'o');
% hold;
% plot(transmitCoordinates(:,1), transmitCoordinates(:,2), '+');
% plot(tgtX, tgtY, 'v');
% axis('square')
% axis([0 4 0 4]);
% xlabel('x');
% ylabel('y');

% Plot for Environment and K
f = figure;

% Maps out X and Y coordinates of K
recXCoordinates = (receiveCoordinates(:, 1))';
recYCoordinates = (receiveCoordinates(:, 2))';

% Sets position and size of Figure
f.Position = [900 500 1300 550];

% Loops through each time step
% for i = 1:length(timeArray)

% Plot for schematic of environment at t = 0
subplot(1,2,1);
plot(receiveCoordinates(:, 1), receiveCoordinates(:, 2), 'o');
hold on;
plot(transmitCoordinates(:, 1), transmitCoordinates(:, 2), '+');
plot(tgtX(1), tgtY(1), 'v');
hold off;
axis('square');
axis([0 4 0 4]);
xlabel('X', FontWeight = "bold");
ylabel('Y', FontWeight = "bold");
title('Environment Schematic');
legend("Receiver", "Transmitter", "Target");

% Plot for real part of K
subplot(1,2,2);
surf(recXCoordinates, recYCoordinates, real(K(:,:,1)))
title(sprintf('3D Surface Plot of Re(K) at %d milliseconds', (1-1) * (timeArray(length(timeArray))-timeArray(1)) / (length(timeArray) - 1)*1000));
xlabel('X', FontWeight = "bold");
ylabel('Y', FontWeight = "bold");
zlabel('SNR', FontWeight = "bold");
axis([1 1.5 0.5 1.5 -0.2 0.2]);

% end
    
% Create a slider
slider = uicontrol('Style', 'slider', 'Min', 1, 'Max', length(timeArray), ...
                   'Value', 1, 'SliderStep', [1/(length(timeArray)-1), 1/(length(timeArray)-1)], ...
                   'Units', 'normalized', 'Position', [0.1, 0.02, 0.8, 0.05], ...
                   'Callback', @(src, event) updatePlot(src, event, K, receiveCoordinates, transmitCoordinates, tgtX, tgtY, timeArray));




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


% Function to update the plot based on the slider value
function updatePlot(source, ~, K, receiveCoordinates, transmitCoordinates, tgtX, tgtY, timeArray)
    timeSlice = round(get(source, 'Value'));
    % surf(real(K(:,:,timeSlice)));
    % xlabel('Transmitter #');
    % ylabel('Reciever #');
    % zlabel('Field');
    % title(['3D Surface Plot of Re(K) at Time Slice ', num2str(timeSlice)]);
    
    % Set fixed axis limits
    % axis([0 30 0 30 -0.2 0.2])

    % Plot for schematic of environment
    subplot(1,2,1);
    plot(receiveCoordinates(:, 1), receiveCoordinates(:, 2), 'o');
    hold on;
    plot(transmitCoordinates(:, 1), transmitCoordinates(:, 2), '+');
    plot(tgtX(timeSlice), tgtY(timeSlice), 'v');
    hold off;
    axis('square');
    axis([0 4 0 4]);
    xlabel('X', FontWeight = "bold");
    ylabel('Y', FontWeight = "bold");
    title('Environment Schematic');
    legend("Receiver", "Transmitter", "Target");
    
    % Plot for real part of K
    subplot(1,2,2);
    surf(receiveCoordinates(:, 1), receiveCoordinates(:, 2), real(K(:,:,timeSlice)))
    title(sprintf('3D Surface Plot of Re(K) at %d milliseconds', (timeSlice-1) * (timeArray(length(timeArray))-timeArray(1)) / (length(timeArray) - 1)*1000));
    xlabel('X', FontWeight = "bold");
    ylabel('Y', FontWeight = "bold");
    zlabel('SNR', FontWeight = "bold");
    axis([1 1.5 0.5 1.5 -0.2 0.2]);


end
