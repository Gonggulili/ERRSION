% First, we set up the simulation environment.
close all
clear all

set(0,'defaultTextFontSize', 14)    % Set default font size for the plots
set(0,'defaultAxesFontSize', 14)

s = simulation_parameters;          % Basic simulation parameters
s.center_frequency = 2.53e9;        % Center Frequency
s.sample_density = 4;               % 4 samples / half wave length

% Setting up the antenna arrays
l=layout(s);                        % create a new layout
l.tx_array.generate('3gpp-3d',2,2)  % create 3gpp-3d 4*4 antenna array for Tx 
l.rx_array.generate('dipole')       % create dipole antenna for Rx

% Setting up the Tx and Rx position and track
l.no_tx = 2;                        % number of Tx
l.no_rx = 2;                        % number of Rx
l.tx_position(:,1) = [10;20;25];    % position of Tx01
l.tx_position(:,2) = [80;90;30];    % position of Tx02

l.track(1).generate('linear', 60 ,3*pi/4 ); % track of Rx01
l.track(1).initial_position = [50;20;0];           % initial position of Rx01
l.track(1).scenario = {'WINNER_UMa_C2_NLOS';'WINNER_Indoor_A1_LOS'};   % Two Scenarios of Rx01 for two Tx
l.track(2).generate('linear', 40,0 );     % track of Rx02
l.track(2).initial_position = [ 100;50;0 ];        % initial position of Rx02
l.track(2).scenario = {'WINNER_UMa_C2_LOS';'WINNER_UMa_C2_NLOS'};   % Two Scenarios of Rx02 for two Tx

l.track.interpolate_positions( s.samples_per_meter );

l.visualize;                                    % Plot all tracks
view(-33,45)

% Generate parameter sets and channel coefficients
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));
p = l.create_parameter_sets;         % Create parameter sets
c = p.get_channels;             % generate channel coefficients