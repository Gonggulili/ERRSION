%% Network Setup and Parameter Generation
%
% The tutorial demonstrates how to setup a simple Layout with multiple
% receivers, how to adjust parameters manually, generate channel
% coefficients, and how to calculate simple parameters from the data.
% 
% The channel model class 'parameter_set' generates correlated values for
% the LSPs. The channel builder then uses those values to create coefficients 
% that have the specific properties defined in 'parameter_set'. One important 
% question is therefore: Can the same properties which are defined in 
% 'parameter_set' also be found in the generated coefficients? This is an 
% important test to verify, if all components of the channel builder work 
% correctly.

%% Channel model setup and coefficient generation
% We first set up the basic parameters. We do not need drifting here, since
% no time varying channels are generated.

close all
clear all

set(0,'defaultTextFontSize', 14)
set(0,'defaultAxesFontSize', 14)

s = simulation_parameters;
s.center_frequency = 2.53e9;
s.sample_density = 2;
s.use_absolute_delays = 1;
s.drifting_precision = 0;

%%
% We have one transmitter and 250 receiver positions. Each receiver gets a
% specific channel. However, the receivers LSPs will be correlated. We use
% omni directional antennas at all terminals.

l = layout( s );
l.no_rx = 250;
l.randomize_rx_positions( 200 , 1.5 , 1.5 , 1 );  % 200 m radius, 1.5 m Rx height
l.track.set_scenario('V2V_Freeway_LOS');

l.tx_position(3) = 25;                            % 25 m tx height
l.tx_array.generate( 'omni' );
l.rx_array = l.tx_array;

l.visualize([],[],0);
view(-33, 60);

% Figure: Distribution of the users in the scenario.

%%
% We set up the scenario such that there is no XPR. I.e. all vertical
% polarized paths will remain vertical after a reflection. The same result
% would be achieved with a perfectly X-polarized array at the receiver and
% summing up the power over all elements. We further increase the KF to
% have a wider spread. This allows us to study the parameters at a wider
% range when evaluating the results.

p = l.create_parameter_sets(0,1);
p.plpar = [];
p.scenpar.xpr_mu    = 100;              % Disable XPR
p.scenpar.xpr_sigma = 0;
p.scenpar.KF_mu     = 5;                % Increase KF-Range
p.scenpar.KF_sigma  = 15;
p.scenpar.DS_mu     = log10(0.6e-6);    % Median DS = 600 ns
p.scenpar.DS_sigma  = 0.3;              % 300-1200 ns range
p.update_parameters;

c = p.get_channels;

coeff = squeeze( cat( 1, c.coeff ) );
delay = permute( cat(3,c.delay) , [3,1,2] );


%% Results and discussion
% In the following four plots, we extract parameters from the generated
% coefficients and compare them with the initial ones which were generated
% by the 'parameter_set' object (P). The values in (P) can be seen as a
% request to the channel builder and the values in the generated
% coefficients (C) as a delivery.
%
% We first calculate the SF from the channel data by summing up the power
% over all 20 taps. We see, that the values are almost identical.

sf = sum(mean( abs(coeff).^2 ,3),2);

figure
plot(-35:35,-35:35,'k')
hold on
plot([-35:35]+3,-35:35,'--k')
plot([-35:35]-3,-35:35,'--k')
plot( 10*log10(p.sf') , 10*log10(sf) , '.'  )
hold off
axis([ -25 , 25 , -25, 25 ])

legend('Equal','+/- 3dB','Location','SouthEast');
xlabel('SF_P [dB]');
ylabel('SF_C [dB]');
title('Shadow Fading - Requested vs. generated value');


%%
% Next, we repeat the same calculation for the K-Factor. Again, we see that
% the values are almost identical.

p_nlos = sum(mean( abs(coeff(:,2:end,:)).^2 ,3),2);
p_los  = mean( abs(coeff(:,1,:)).^2 ,3);
kf = p_los./p_nlos;

figure
plot(-35:35,-35:35,'k')
hold on
plot([-35:35]+3,-35:35,'--k')
plot([-35:35]-3,-35:35,'--k')
plot( 10*log10(p.kf') , 10*log10(kf) , '.'  )
hold off
axis([ -30 , 30 , -30, 30 ])

legend('Equal','+/- 3dB','Location','SouthEast');
xlabel('KF_P [dB]');
ylabel('KF_C [dB]');
title('K-Factor - Requested vs. generated value');

%%
% Now we repeat the calculation for the RMS delays spread.

pow_tap = abs(coeff).^2;
pow_sum = sum(pow_tap,2);
mean_delay = sum( pow_tap.*delay,2) ./ pow_sum;
ds = sqrt( sum( pow_tap.*delay.^2 ,2)./ pow_sum - mean_delay.^2 );
ds = mean(ds,3);

figure
plot([0:0.1:2],[0:0.1:2],'k')
hold on
plot([0:0.1:2]*1.1,[0:0.1:2],'--k')
plot([0:0.1:2],[0:0.1:2]*1.1,'--k')
plot( p.ds'*1e6 , (ds')*1e6 , '.'  )
hold off
axis([ 0,1.5,0,1.5 ])

legend('Equal','+/- 10% Error','Location','SouthEast');
xlabel('DS_P [\mus]');
ylabel('DS_C [\mus]');
title('Delay Spread - Requested vs. generated value');

%%
% The following plot shows the RMSDS of the requested and generated values
% (in dB) vs. the K-factor. A value of +3 means, that the RMSDS of the 
% generated coefficients is twice a high as in the parameter_set (P). We see,
% that for a K-Factor of up to 30 dB, the DS difference is small (less than
% 3 dB).

figure
plot([-35,35],[0,0],'k')
hold on
plot([-35,35],[-3,-3],'--k')
plot([-35,35],[3,3],'--k')
plot( 10*log10(p.kf'), 10*log10(ds./p.ds') , '.' )
hold off
axis([ -30 , 30 , -6 6 ])

legend('Equal','+/- 3dB','Location','SouthWest');
xlabel('KF_P [dB]');
ylabel('DS_P - DS_C [dB]');
title('Delay Spread difference vs. K-factor');

%%
%close all
disp(['QuaDRiGa Version: ',simulation_parameters.version])
