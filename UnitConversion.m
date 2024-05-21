% Script to get people acquainted with FK coronae rescaling using the Fast
% Discrete Curvelet Transform
% James Atterholt Caltech 2024
clear

% Add path to the CurveLab Toolbox (put your own path in here). You can
% download this software package at www.curvelet.org
addpath("~/SeismoPrograms/CurveLab-2.1.3/fdct_wrapping_matlab/")

% Add path to your favorite seismic colormap. Here I use bluewhitered.
% Link: https://www.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered
addpath("~/SeismoPrograms/ColorMaps/")

% References:
% E. J. Candes, L. Demanet, D. L. Donoho, L. Ying, Fast Discrete Curvelet 
%   Transforms, 2005.
% E. J. Candes, D. L. Donoho, New Tight Frames of Curvelets and Optimal 
%   Representations of Objects with Smooth Singularities, 2002.

%% Please cite:
% Atterholt, J., Zhan, Z., Shen, Z., & Li, Z. (2021). A unified 
%   wavefield-partitioning approach for distributed acoustic sensing. 
%   Geophysical Journal International, 228(2), 1410–1418. 
%   https://doi.org/10.1093/gji/ggab407

%% Load in some strain rate data

MAT = load("EQ_strainrate.mat").strainrate_MAT;

%% Do some preprocessing

% Define a priori variables
dt = (1/100); % s (sampling rate)
ds = 8; % m (station spacing)

% Define bad and good stations for this array (hand picked)
BadSta = [1:24 83:90 143:147 256:271 455:473 663:670 843:872 1060:1067];
GoodSta = setdiff(1:size(MAT,2),BadSta);

% Get rid of the bad stations
MAT = MAT(:,GoodSta);

% Get the parameters for the filtering 
npts=size(MAT,1);
nch=size(MAT,2);

% Set high and low frequencies for the filter
freq1=1;
freq2=10;

% Create a butterworth filter and a mild taper
Fn=1/dt/2;
w1=freq1/Fn;
w2=freq2/Fn;
[B,A]=butter(2,[w1 w2],'bandpass');
taper=tukeywin(npts,0.05);

% Perform the filtering along with detrending. Tapers are important for
% edge effects. Filtering can be turned off here.
for ii=1:nch
    MAT(:,ii)=detrend(MAT(:,ii));
    MAT(:,ii)=filtfilt(B,A,MAT(:,ii).*taper); 
    MAT(:,ii)=MAT(:,ii).*taper;
end

% Apply a median filter (a,b,c selected empirically, see function for
% details). Median filter reduces spurious peaks / discontinuities that
% create artifacts when transforming.
MAT = MedianFilter(MAT,50,5,10);

% Apply a tapering to the left and right sides too
taper_horiz=tukeywin(length(GoodSta),0.05);
for i = 1:size(MAT,1)
    MAT(i,:) = MAT(i,:) .* taper_horiz';
end

% Get distance and time
t = (1:size(MAT,1))*dt;
d = (1:size(MAT,2))*ds;

%% Plot the preprocessed data
figure_1 = figure(1);
figure_1.Position(3:4) = [900 700];
imagesc(d,t,MAT);
ylim([2 18])
caxis([-1 1]*2.5e-7)
colormap(bluewhitered(256))
colorbar
xlabel("Distance (m)")
ylabel("Time (s)")
title("Strain-Rate (strain/s)")

%% Curvelet transform
% 1st: The matrix that you are transforming
% 2nd: Real or complex valued curvelets (0 = complex, 1 = real)
% 3rd: Objects at the finest scale, I choose curvelets here to keep the
% velocity compartmentalization at the finest scale.
% 4th and 5th: These define the structure of the polar tiling you are
% using (4th: number of scales, 5th: number of angles at the coursest scale). 
% We end up in the weeds here. I choose 8 and 16 because those are
% optimized sizes for the example and most DAS arrays. Number of angles may
% be increased for more finely compartmentalized velocity bounds.
C_i = fdct_wrapping(MAT,1,1,8,16);

f = 1/dt; % 1/s (frequency)
k = 1/ds; % 1/m (wavenumber)

% Allocate for the calculation of the mean for each wedge
Wedge_velocities = cell(size(C_i));
% Zero the 1st wedge (contains velocities 0-inf). It's super low frequency.
Wedge_velocities{1} = 0; 

% Now run through each cell and wedge and compute the velocity for
% each wedge. Velocities are those that bisect each wedge.
for s = 2:length(C_i)
    v_bounds = zeros((length(C_i{s})/4)+1,1);
    halfway = length(C_i{s})/8;
    v_bounds(halfway+1,1) = (f)/(k);
    for i = 1:(halfway)
        % Sorry about the parentheses here...
        % Velocities approaching the f/k value
        v_bounds(i,1) = ((((i-1)/halfway)+((1/2)*(1/halfway)))*f)/(k);
        % Velocities exceeding the f/k value
        v_bounds(halfway+i+1,1) = (f)/(((((halfway-i)/halfway+(1/2)*(1/halfway))))*k);
    end
    % Ensuring that the 8 half quadrants are accounted for. Correct sign
    % for the conversion is accounted for here. 
    wedgevels_use1 = (halfway+2):halfway*2+1; 
    wedgevels_use2 = halfway*2+1:-1:(halfway+2); 
    wedgevels_use3 = halfway:-1:1;
    wedgevels_use4 = 1:halfway;
    wedgevels_shaped = [v_bounds(wedgevels_use1) ; -v_bounds(wedgevels_use2) ; -v_bounds(wedgevels_use3) ; v_bounds(wedgevels_use4) ; v_bounds(wedgevels_use1) ; -v_bounds(wedgevels_use2) ; -v_bounds(wedgevels_use3) ; v_bounds(wedgevels_use4)];
    
    % Store the bisecting velocities in the cell shaped like the curvelet
    % transform version of your data.
    Wedge_velocities{s} = wedgevels_shaped;
end
%%
% Get a copy of C_i.
C_i_copy = C_i;

% Start your loop towards multiplying the velocities.
for s=1:length(C_i)
    for w = 1:length(C_i{s})
        
        % Same halfway as before (# of wedges from 0 to f/k).
        halfway = length(C_i{s})/8;
        
        % Classify wedges that contain infinite velocities.
        infinites = [halfway, halfway+1, halfway + length(C_i{s})/2, halfway + length(C_i{s})/2 + 1];
        
        % Perform the multiplication. Wedges with infinite velocities are
        % zeroed for stability.
        if ((s > 1) && ismember(w,infinites))
            C_i_copy{s}{w} = C_i{s}{w}(:,:)*0;
        else
            C_i_copy{s}{w} = C_i{s}{w}(:,:)*Wedge_velocities{s}(w);
        end
        
        
    end
end

% Perform the inverse curvelet transform
acceleration_MAT = real(ifdct_wrapping(C_i_copy,1,size(MAT,1),size(MAT,2)));

%% Plot up your acceleration seismogram
figure_2 = figure(2);
figure_2.Position(3:4) = [900 700];
imagesc(d,t,acceleration_MAT);
caxis([-1 1]*5e-4);
ylim([2 18])
colormap(bluewhitered(256))
colorbar
xlabel("Distance (m)")
ylabel("Time (s)")
title("Acceleration (m/s^2) (infinity wedges removed)")

%% Save to MAT file

save('EQ_acceleration.mat','acceleration_MAT')
