%% How to manually set LSPs in QuaDRiGa (Satellite Scenario)
%
% This tutorial explains, how to generate a time-series of channel
% coefficients with manual selection of LSPs in a satellite scenario.
%
% By default, QuaDRiGa automatically generates correlated LSPs based on
% statistics extracted from measurements. However, in some cases it might
% be preferable to fix these parameters, e.g. when the large-scale fading
% is provided by an external source such as a measured profile. This can be
% done by providing specific valued for the LSPs along with the terminal
% trajectory.

%% Setting general parameters
% We set up some basic parameters such as center frequency and sample density.

close all
clear all

set(0,'defaultTextFontSize', 14)    % Set default font size for the plots
set(0,'defaultAxesFontSize', 14)

s = simulation_parameters;          % Basic simulation parameters
s.center_frequency      = 2.185e9;  % Center Frequency
s.sample_density        = 4;        % 4 samples / half wave length


%% Defining Antenna Arrays
% We set up our antenna arrays for the transmitter at the satellite and the
% receiver. We use synthetic patch antennas for this case. Two elements
% are crossed by an angle of 90 degree. The signal is then split and fed
% with a 90 degree phase shift to both elements generating RHCP and 
% LHCP signals.

% Create a patch antenna with 120 degree opening
a = array('custom',120,120,0);

% Copy element 1 to element 2.
a.copy_element(1,2);

% Rotate the second pattern by 90 degree around the x-axis.
a.rotate_pattern(90,'x',2);

%%
% Set the coupling between the elements. The Tx-signal for the first
% element is shifted by -90 degree out of phase and put on the second element.
% The signal for the second element is shifted by +90 degree and copied to the
% first element. Both antennas thus radiate a LHCP and a RHCP wave.

a.coupling = 1/sqrt(2) * [1 1;1j -1j];  % LHCP and RHCP

% Create a copy of the array for the receiver.
b = a.copy_objects;

% Rotate the receive antenna array to face sky-wards.
b.rotate_pattern(-90,'y');

b.visualize;                    % Plot the pattern of the Rx-Antenna

%% Setting the satellite orbital position
% In this step, we set the antennas and the position of the
% satellite into a simulation layout. A layout object contains all the
% geometric information that are necessary to run the simulation. First,
% we define the position of the satellite. Since the model uses Cartesian
% coordinates, we have to transform the position of the satellite first.
%
% Choose a random satellite position (Astra 2, seen from Berlin).
% The distance only needs to be big enough to ensure insignificant changes
% in the reception angle on the ground.

sat_el      = 28.4;             % Elevation angle
sat_az      = 161.6;            % Azimuth angle (South = 180 degree)
rx_latitude = 51;               % Latitude of the Rx

% Approximate the satellite distance for GEO orbit
dist_x      = 35786 + rx_latitude/90 * 6384;    % [km]
dist_y      = (1-rx_latitude/90) * 6384;        % [km]
sat_dist    = sqrt(dist_x^2 + dist_y^2);        % [km]
sat_dist    = sat_dist*1e3;                     % [m]

% Transform angles to Cartesian coordinates
sat_x = sat_dist * cosd(sat_el) * cosd( -sat_az+90 );
sat_y = sat_dist * cosd(sat_el) * sind( -sat_az+90 );
sat_z = sat_dist * sind(sat_el);

% We also turn the antenna of the satellite so it points to the receiver.
a.rotate_pattern( sat_el , 'y' );
a.rotate_pattern( 270-sat_az , 'z' );

% Create a new layout
l = layout( s );            

% Set the satellite position in the layout
l.tx_position = [ sat_x ; sat_y ; sat_z ];

% Set the transmit antena in the layout
l.tx_array = a;

l.randomize_rx_positions( 500,0,0,0 );          % Random Rx position


%% Setting up a user track and assigning parameters
% QuaDRiGa needs the positions of transmitter and receiver, e.g. for
% calculating the polarization or the arrival and departure angels. The
% positioning information of the Tx and Rx is essential also when the LSPs
% are not calculated. The following generates a linear track with 20 m
% length having a random direction. The track is further split into 4 segments of
% 5 m length.
%
% The splitting is done by calling the method 'split_segment' of the track
% object. The first two arguments of that function are the minimum
% length of a segment (1 m) and the maximum length of the segment (6 m).
% Each existing segment that is longer then the maximum length is split
% into subsegments. The length of those segments is random where the third
% and fourth parameter determine the mean and the STD of the length of new
% subsegment. Hence, 't.split_segment( 1,6,5,0 )' splits all segment longer
% than 6 m into subsegments of 5 m length.

t = track('linear',20);                         % Linear track, 20 m length
t.interpolate_positions( s.samples_per_meter);  % Interpolate to sample density
t.split_segment( 1,6,5,0 )                      % Split in 4 segments

t.initial_position = l.rx_position;

%%
% Each segment gets assigned a scenario. This is also essential since many
% parameters (such as the number of clusters, the XPR etc.) are
% scenario-specific. Hence, they are the same for the entire scenario.
% Here, we set the first the segments to NLOS, the third to LOS and the
% last to NLOS.

Sn = 'MIMOSA_10-45_NLOS';
Sl = 'MIMOSA_10-45_LOS';
t.scenario = {Sn,Sn,Sl,Sn};                     % Set scenarios

%%
% Parameters are assigned to the track by using the "par" field. There are
% 8 parameters that can be set:
% 
% 	Delay Spread (ds) in units of [sec], 
% 	Ricean K-Factor (kf) in [dB], 
% 	Path gain (pg) in [dB], 
% 	Azimuth spread of Departure (asD) in [degree], 
% 	Azimuth spread of Arrival (asA) in [degree], 
% 	Elevation spread of Departure (asD) in [degree], and
% 	Elevation spread of Arrival (asA) in [degree]
%   Cross-Polarization Ratio (XPR) in [DB]
% 
% Delay Spread (ds), the angular spreads (asD, asA, esD and esA) and the
% XPR are given once for each segment of the track. K-Factor (kf) and Path
% Gain (pg) are given for each snapshot. QuaDRiga can fill the fields
% automatically using the map-based parameter generation method.   

t.generate_parameters

%% 
% Now, we want to adjust the delay spread manually. We do this by editing
% the par-field of the track. Note that missing parameters are automatically 
% generated by the channel model.

t.par.ds = [ 0.45 0.33 0.12 0.60 ]*1e-6;

%%
% Next, we manually set the path gain. We create a time-series for the PG
% using a cubic interpolation method. This is also set in the parameter
% object.  

pg = [ -102 , -97 , -82 , -99 ];                % Path gain per segment
pg = [ pg ; pg ];

ind = [ t.segment_index(2:end)-50 ; t.segment_index(2:end)+50 ];
ind = [ 1 ; ind(:) ; t.no_snapshots ];          % Segment positions

pgi = pchip( ind , pg(:) , 1:t.no_snapshots );  % Cubic interpoaltion

t.par.pg = pgi;


%% Calculate channel coefficients
% Now we calculate the coefficients and the received power along the path.
% The following command calculate the channel coefficients. We then check
% the number of clusters that have been produced for each segment.

l.track = t;                % Set track
l.rx_array = b;             % Set Rx-Antenna

c = l.get_channels;         % Calculate the channel coefficients


%%
% We plot the power along the path. You can see the power levels of
% around -102, -97, -82 and -99 dBm which have been set earlier. Note that
% there are additional antenna gains.

power = sum(abs(c.coeff).^2,3);
power = 10*log10(power);
power = reshape( power, [] , size(power,4) ).';

figure
[~,dist] = t.get_length;
plot(dist,t.par.pg,'--k','Linewidth',2)
hold on
plot(dist,power)
hold off
title('Simulated Power')
xlabel('Distance from start point [m]')
ylabel('Received Power [dBm]')
axis([0,20,-115,-65])
grid on
legend('Setting','LHCP-LHCP','LHCP-RHCP','RHCP-LHCP','RHCP-RHCP',2)


%%
% The last plot shows the DS along the path. The results reflect the
% settings of 0.45, 0.33, 0.12 and 0.60 quiet well. As for the power, there
% is an overlap in between the segments. For example, in between 7.5 and
% 10m the DS drops from 0.33 to 0.12 microseconds. Additional fluctuations
% are caused by small scale fading.

coeff = c.coeff;
delay = c.delay;

pow_tap = squeeze( mean(mean(abs( coeff ).^2,1),2) );
pow_sum = sum(pow_tap);
mean_delay = sum( pow_tap.*delay) ./ pow_sum;
ds = sqrt( sum( pow_tap.*delay.^2)./ pow_sum - mean_delay.^2 );

dss = zeros( 1,t.no_snapshots);
ind = [ t.segment_index t.no_snapshots ];
for n = 1:t.no_segments
    dss(ind(n) : ind(n+1)) = t.par.ds(n);
end

figure
plot(dist,dss*1e6,'--k','Linewidth',2)
hold on
plot(dist,ds*1e6)
hold off
title('Simulated Delay Spread')
xlabel('Distance from start point [m]')
ylabel('RMS DS [\mus]')
axis([0,20,0,1])
grid on
legend('Setting','Actual DS',2)

%%
close all
disp(['QuaDRiGa Version: ',simulation_parameters.version])

