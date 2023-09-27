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
%        |       R [RECEIVE ARRAY]
%        |      R
%        P     R
%        E    R
%        C   R
%        | 
%        |              T
% Mirror3|               T
%        |                T     [TRANSMIT ARRAY]
%        |  *              T
%        |*                 T
%        |     *
%        |  [MULTIPLE TARGET POSITIONS AS A FUNCTION OF TIME]
%        |         *
%        |
% -------------PEC plate -----------------------------------------
% Mirror2|   Mirror1

% Conditions to account for:
%   1) Field between any given transmitter and receiver
%   2) Field between any given transmitter and target
%   3) Field between and given receiver and targer



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
%   4) Array Length
%   5) Angle of array
[transmitCoordinates] = line_array(0.5, 2, 30, 2, pi/6);
[receiveCoordinates] = line_array(2.5, 0.5, 30, 2, (pi/2 + pi/6));

% tgtX, tgtY: they indicate the targets' positions. 
% Also, we have the target strengths tau. 

tgtX = 0.2;     % Update as a function of time
tgtY = 0.8;     % Update as a function of time
tau = 1;

% frequency range
num_frequencies = 225;

% Creates array of angular frequencies (rad/s) ranging from 1.04c(pi) [lambda ~= 1.92] to 10c(pi) [lambda = 0.2]
w = zeros(1,num_frequencies);
for index = 1:length(w)
    w(index) = (1 + index/25) * c * pi;
end



%% CALCULATE TOTAL FIELD
for UU = 1:length(transmitCoordinates)
    for VV = 1:length(receiveCoordinates)
        transX = transmitCoordinates(UU, 1);
        transY = transmitCoordinates(UU, 2);
        recX = receiveCoordinates(VV, 1);
        recY = receiveCoordinates(VV, 2);

        % Distance between transmitter and receiver
        d_trans_rec = sqrt((transX - recX)^2 + (transY - recY)^2);
        % Distance between transmitter and receiver mirror images
        d_trans_recMirror(1) = sqrt((transX - recX)^2 + (transY + recY)^2);
        d_trans_recMirror(2) = sqrt((transX + recX)^2 + (transY + recY)^2);
        d_trans_recMirror(3) = sqrt((transX + recX)^2 + (transY - recY)^2);

        % for TT = number of time elements
            % Lines 166 - 193 from Corner_Sept8_V3?
            
            % % Distance between transmitter and target
            % d_trans_tgt(ZZ) = sqrt((transX - tgtX(ZZ))^2 + (transY - tgtY(ZZ))^2);
            % % Distance between transmitter and target mirror images
            % d_trans_tgtMirror1(ZZ) = sqrt((transX - tgtX(ZZ))^2 + (transY + tgtY(ZZ))^2);
            % d_trans_tgtMirror2(ZZ) = sqrt((transX + tgtX(ZZ))^2 + (transY + tgtY(ZZ))^2);
            % d_trans_tgtMirror3(ZZ) = sqrt((transX + tgtX(ZZ))^2 + (transY - tgtY(ZZ))^2); 
            % 
            % % Distance between target and receiver
            % d_tgt_rec(ZZ) = sqrt((recX - tgtX(ZZ))^2 + (recY - tgtY(ZZ))^2);
            % % Distance between target and receiver mirror images (target acts as emitter)
            % d_tgt_recMirror1(ZZ) = sqrt((recX - tgtX(ZZ))^2 + (recY + tgtY(ZZ))^2);
            % d_tgt_recMirror2(ZZ) = sqrt((recX + tgtX(ZZ))^2 + (recY + tgtY(ZZ))^2);
            % d_tgt_recMirror3(ZZ) = sqrt((recX + tgtX(ZZ))^2 + (recY - tgtY(ZZ))^2);    

            for WW = 1:length(w)
                % Background field transmitter to receiver
                U_bknd_comps(1) = greens(d_trans_rec, w(WW), c);
                U_bknd_comps(2) = -1 * greens(d_trans_recMirror(1), w(WW), c);
                U_bknd_comps(3) = greens(d_trans_recMirror(2), w(WW), c);
                U_bknd_comps(4) = -1 * greens(d_trans_recMirror(3), w(WW), c);
                U_bknd(UU,VV,WW) = U_bknd_comps(1) + Reverb*(U_bknd_comps(2) + U_bknd_comps(3) + U_bknd_comps(4));    

                % % Incident field transmitter to target
                % U_inc_comps(1) = greens(d_trans_tgt(ZZ), w(WW), c);
                % U_inc_comps(2) = -1 * greens(d_trans_tgtMirror1(ZZ), w(WW), c);
                % U_inc_comps(3) = greens(d_trans_tgtMirror2(ZZ), w(WW), c);
                % U_inc_comps(4) = -1 * greens(d_trans_tgtMirror3(ZZ), w(WW), c);
                % U_inc(UU, ZZ, WW) = U_inc_comps(1) + Reverb*(U_inc_comps(2) + U_inc_comps(3) + U_inc_comps(4));  
                % 
                % % Field at receiver not influenced by target or tau
                % U_out_comps(1) = greens(d_tgt_rec(ZZ), w(WW), c);
                % U_out_comps(2) = -1 * greens(d_tgt_recMirror1(ZZ), w(WW), c);
                % U_out_comps(3) = greens(d_tgt_recMirror2(ZZ), w(WW), c);
                % U_out_comps(4) = -1 * greens(d_tgt_recMirror3(ZZ), w(WW), c);
                % U_out(UU, ZZ, WW) = U_out_comps(1) + Reverb*(U_out_comps(2) + U_out_comps(3) + U_out_comps(4));      
            end
    % end
    end
end



%% PLOTS

% Plot for schematic of environment
plot(receiveCoordinates(:,1), receiveCoordinates(:,2), 'o');
hold;
plot(transmitCoordinates(:,1), transmitCoordinates(:,2), '+');
plot(tgtX, tgtY, 'v');
axis('square')
axis([0 4 0 4]);
xlabel('x');
ylabel('y');

a = 1;


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

    result = -1 / (4 * pi) * 1 / distance * exp(i * (frequency / c) * distance);
end
