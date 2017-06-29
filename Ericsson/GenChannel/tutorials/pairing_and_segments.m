%% Tutorial - Setting up pairing and segments
%
% This quick tutorial shows how to set up scenarios with several
% transmitters and receiveres and the use of scenarios. 
% First, we set up a basic simulation with two transmitters.

clear all

s = simulation_parameters;
s.center_frequency = 2.53e9;
s.drifting_precision = 0;

l = layout(s);

l.no_tx = 2;
l.tx_position(:,1) = [ -142 ; 355 ; 64 ];   % Outdoor TX
l.tx_position(:,2) = [ 5 ; 0; 10 ];         % Indoor Tx

l.no_rx = 2;

%%
% Second, we set two different receivers.
% Track 1 is indoors. The Link to Tx1 is in scenario Un
% The link to Tx2 is in A1l. This track has no segments. The rows in 
% track.scenario indicate the scenario for each Tx. If there is only one
% row, then all Tx get the same scenario.

l.track(1).generate('linear', 0.2 );
l.track(1).scenario = {'WINNER_UMa_C2_NLOS';'WINNER_Indoor_A1_LOS'};   % Two Scenarios

%%
% The second track is outdoors, far away from the Indoor transmitter. 
% The first part of track 2 is in LOS, the second is in NLOS. The columns
% of track.scenario indocate the segments. Here, all Tx get the same
% scenarios.

l.track(2).generate('linear', 0.2 );
l.rx_position(:,2) = [ 100;50;0 ];
l.track.interpolate_positions( s.samples_per_meter );

l.track(2).segment_index = [1 3];
l.track(2).scenario = {'WINNER_UMa_C2_LOS','WINNER_UMa_C2_NLOS'};

l.visualize;
view(-33, 60);
%%
% We calculate the channel coefficients and plot the list of created
% segments.

p = l.create_parameter_sets;
c = p.get_channels;

strvcat( c.name )


%% 
% As we can see, 6 segments were generated. However, the channel Tx2_Rx2
% will most likely not be needed because of the large distance. 
% We thus remove the link from the pairing matrix and recompute the
% channels.

l.pairing = [1 2 1 ; 1 1 2 ];

p = l.create_parameter_sets;
c = p.get_channels;

strvcat( c.name )

%% 
% At last, we can combine the segments and generate the final chanels.

cn = c.merge;
strvcat( cn.name )

%%
close all
disp(['QuaDRiGa Version: ',simulation_parameters.version])
