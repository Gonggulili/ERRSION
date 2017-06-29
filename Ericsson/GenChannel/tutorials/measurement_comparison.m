%% Simulating a Measured Scenario
%
% This more complex tutorial shows how to manually define a state sequence,
% manipulate antennas, run the PPG and obtain a time series of channel
% coefficients.
%
% This script recreates a measured drive test from the Park Inn Hotel at
% Berlin Alexanderplatz. The transmitter was at the rooftop of the hotel
% while the mobile receiver was moving south on Grunerstraße. A simplified
% version of the scenario is recreated in the simulation where the
% scenarios along the track were classified by hand.    


%% Channel model setup and coefficient generation
% First, we set up the channel model. 

set(0,'defaultTextFontSize', 14)
set(0,'defaultAxesFontSize', 14)
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));

close all
clear all

s = simulation_parameters;                      % Basic simulation parameters
s.center_frequency = 2.185e9;
s.sample_density = 2;
s.use_absolute_delays = 1;

t = track('linear',500,-135*pi/180);            % Track with 500m length, direction SE
t.initial_position = [120;-120;0];              % Start position
t.interpolate_positions( 1 );                   % Interpolate to 1 sample per meter

t.segment_index = [1    45   97   108  110  160  190  215  235  245 ...
    280  295  304  330  400  430 ];             % Set segments (states)

Sl = 'MIMOSA_10-45_LOS';
Sn = 'MIMOSA_10-45_NLOS';
t.scenario      = {Sn,Sl,Sn,Sl,Sn,Sn,Sn,Sl,Sn,Sl,Sn,Sl,Sn,Sn,Sn,Sn};
t.interpolate_positions( 3 );

l = layout( s ); 
l.tx_position = [0;0;125];                      % Set the position of the Tx
l.track = t;                                    % Set the rx-track

l.tx_array = array('lhcp-rhcp-dipole');         % Generate Tx antenna
l.tx_array.rotate_pattern(30,'y');              % 30 deg Tilt
l.tx_array.rotate_pattern(-90,'z');             % point southwards

l.rx_array = array('lhcp-rhcp-dipole');         % Rx-Antenna
l.rx_array.rotate_pattern(-90,'y');             % point skywards

l.visualize;
view(-33, 45);

lnk = [ l.tx_position, ...
    l.track.positions(:,l.track.segment_index(2))+l.track.initial_position ];

hold on
plot3( lnk(1,:),lnk(2,:),lnk(3,:) , '--' )
hold off

%% Generate channel coefficients
% Next, we calculate the channel coefficients.

[p,cb] = l.create_parameter_sets(0);
p(2).scenpar.NumClusters = 14;
p.update_parameters;

c = cb.get_channels;
cn = c.merge(0.2);


%% Results
% First, we plot the PDP vs distance from the start point. 

h = cn(1).fr( 20e6,256 );
pdp = squeeze(sum(sum( abs(ifft(h,[],3)).^2 , 1),2));
pdp = 10*log10(pdp.');

figure
imagesc(pdp(end:-1:1,1:192));

cm = colormap('hot');
colormap(cm(end:-1:1,:));

caxis([ max(max(pdp))-60 max(max(pdp))-5 ]);
colorbar

title('Time variant power delay profile');

set(gca,'XTick',1:32:192);
set(gca,'XTickLabel',(0:32:192)/20e6*1e6);
xlabel('Delay [\mus]');

ind = sort( cn.no_snap : -cn(1).no_snap/10 : 1);
set(gca,'YTick', ind );
set(gca,'YTickLabel', round(sort(500-ind / 3,'descend')) );
ylabel('Distance [m]');


%%
% The next plot shows the total received power along the path. Green shaded
% ares are LOS. The rest is NLOS.

dist = (1:cn.no_snap)*l.track.get_length/cn.no_snap;
ind  = find(strcmp(l.track.scenario,Sl));
los  = [];
for n = 1:numel(ind)
    los = [los l.track.segment_index(ind(n)) : l.track.segment_index(ind(n)+1)];
end

power = 10*log10( sum( reshape( abs(cn.coeff).^2 , [] , cn.no_snap ) ,1)/4 );
ar = zeros(1,cn.no_snap);
ar(los) = -200;

figure
a = area(dist,ar);
set(a(1),'FaceColor',[0.7 0.9 0.7]);
set(a,'LineStyle','none') 

hold on
plot(dist,power)
hold off

title('Position dependent power')
xlabel('Track [m]');
ylabel('Power [dB]');
axis([0 500 min(power)-5 max(power)+5])
legend('LOS','P_{total}','location','SouthEast');
grid on


%%
% The following plot shows the distribution (PDF) of the received power for
% both, the LOS and NLOS segments.

bins   = -150:2:-80;
p_los  = hist(power(los),bins)/cn.no_snap*100;
p_nlos = hist(power(setdiff(1:cn.no_snap,los)),bins)/cn.no_snap*100;

figure
bar(bins,[p_los;p_nlos]')
axis([-124.5,-83,0,ceil(max([p_los,p_nlos]))])
grid on
colormap('Cool')

title('Empirical PDF of the LOS and NLOS power')
xlabel('P_{total} [dB]');
ylabel('Probability [%]');
legend('LOS','NLOS','location','NorthEast')


%%
% The next plot shows the RMS delay spread along the path. Again, shaded
% ares are for the LOS segments.

pow_tap = squeeze(sum(sum(abs(cn.coeff).^2,1),2));
pow_sum = sum( pow_tap,1 );
mean_delay = sum( pow_tap.*cn.delay ,1) ./ pow_sum;
ds = sqrt( sum( pow_tap.*cn.delay.^2 ,1)./ pow_sum - mean_delay.^2 );
ar = zeros(1,cn.no_snap);
ar(los) = 10;

figure
a = area(dist,ar);
set(a(1),'FaceColor',[0.7 0.9 0.7]);
set(a,'LineStyle','none') 

hold on
plot( dist , ds*1e6  )
hold off

ma = 1e6*( max(ds)+0.1*max(ds) );

axis([0 500 0 ma])
title('Position dependant delay spread');
xlabel('Track [m]');
ylabel('Delay Spread [dB]');
legend('LOS','\sigma_\tau','location','NorthEast')

grid on


%%
% The following plot shows the distribution (PDF) of the RMS delay spread for
% both, the LOS and NLOS segments.

bins = 0:0.03:3;
ds_los  = hist(ds(los)*1e6,bins)/cn.no_snap*100;
ds_nlos = hist(ds(setdiff(1:cn.no_snap,los))*1e6,bins)/cn.no_snap*100;

figure
bar(bins,[ds_los;ds_nlos]')
axis([0,1.5,0,ceil(max([ds_los,ds_nlos]))])
grid on
colormap('Cool')

title('Empirical PDF of the LOS and NLOS RMSDS')
xlabel('\sigma_\tau [\mus]');
ylabel('Probability [%]');
legend('LOS','NLOS','location','NorthEast')


%%
% The Quadriga model calculates specular components for each cluster. Here
% we visualize these components as seen by the receiver.  

p(2).no_positions = 1;
p(2).rx_track = p(2).rx_track(1);
p(2).rx_array = p(2).rx_array(1);
p(2).rx_track.interpolate_positions( s.samples_per_meter );

[ c_detailed ,cb ] = p(2).get_channels;


%% 
% The following plots are for the first LOS segment on the track and show
% the dependency of the received power on the Delay, AoA and EaA at the
% receiver position.     

delays = c_detailed.delay(:,1).'*1e6;
pow = 10*log10( squeeze( sum(mean( abs(c_detailed.coeff(:,:,:,1)).^2 ,1),2) )' );
az = cb.AoA(1,:) * 180/pi;
el = cb.EoA(1,:) * 180/pi;

figure
plot(delays,pow,'*' ,'Markersize',10)
xlabel('Delay [\mus]')
ylabel('Power [dB]')
title('Power versus Delay Spread (LOS segment 1)')
grid on
a = axis;
axis([0,3,a(3:4)])

figure
plot(pow,el,'*' ,'Markersize',10)
ylabel('Elevation angle [deg]')
xlabel('Power [dB]')
title('Power versus Elevation (LOS segment 1)')
grid on
a = axis;
axis([a(1:2),-90,90])


figure
plot(az,pow,'*','Markersize',10)
hold on
plot(az,10*log10(abs(squeeze(c_detailed.coeff(1,1,:,1))).^2)','or' ,'Markersize',10)
plot(az,10*log10(abs(squeeze(c_detailed.coeff(1,2,:,1))).^2)','og','Markersize',10 )
hold off
xlabel('Azimuth angle [deg]')
ylabel('Power [dB]')
title('Power versus Azimuth (LOS segment 1)')
legend('Sum Power','P(H11)','P(H12)','location','NorthWest');
grid on
a = axis;
axis([-180,180,a(3:4)])


pow_total = 10*log10(squeeze(sum( mean(mean(abs( c_detailed.coeff ).^2,1),2) ,3)));
pow_direct = 10*log10(squeeze( mean(mean(abs( c_detailed.coeff(:,:,1,:) ).^2,1),2) ));

figure
plot( 1:c_detailed.no_snap , pow_total, 'r' );
hold on
plot( 1:c_detailed.no_snap , pow_direct ,'b');
hold off
grid on

title('Power versus Sanpshot index (LOS segment 1)')
xlabel('Snapshot Index')
ylabel('Power [dB]')
legend('Total Power','Direct Component')

axis([ 1, c_detailed.no_snap , min(pow_total)-5 ,max(pow_total)+5   ])


pd = squeeze(abs(c_detailed.coeff(:,:,1,:)).^2);
xpd_d = squeeze(pd(1,1,:)+pd(2,2,:)) ./ squeeze(pd(1,2,:)+pd(2,1,:));

pt = squeeze(sum( abs( c_detailed.coeff ).^2 ,3));
xpd_t = squeeze(pt(1,1,:)+pt(2,2,:)) ./ squeeze(pt(1,2,:)+pt(2,1,:));

%%
%close all
disp(['QuaDRiGa Version: ',simulation_parameters.version])
