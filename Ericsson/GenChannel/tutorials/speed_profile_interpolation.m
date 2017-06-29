%% Applying Varying Speeds (Channel Interpolation)
%
% This tutorial shows how to adjust the speed of the terminal, e.g. when
% breaking or accelerating. First, a simple scenario defined. Channel
% coefficients are calculated at a constant speed and then interpolated to
% match the varying speed of the terminal.
%
% One new feature that makes the simulations more realistic is the function
% to apply arbitrary speed- and movement profiles, e.g. accelerating,
% breaking or moving at any chosen speed. These profiles are defined in the
% track class. The profiles are then converted in to effective sampling
% points which aid the interpolation of the channel coefficients.    

%% Channel model setup
% First, we set up the simulation parameters. Note the sample density of 2.5
% which enables very fast simulations even with drifting.

close all
clear all

set(0,'defaultTextFontSize', 14)
set(0,'defaultAxesFontSize', 14)

s = simulation_parameters;
s.center_frequency = 2.53e9;
s.sample_density = 2.5;
s.use_absolute_delays = 1;

%%
% Second, we define a track. It has a length of 20 m, starts at 10 m east
% of the transmitter and consists of three segments (LOS, NLOS, LOS). The
% positions are interpolated to match the sample density defined above. The
% track is then plugged into a network layout with one transmitter at
% position (0,0,25). Both, transmitter and receiver are equipped with
% dipole antennas. The last three lines create the large scale parameters
% (LSPs).      

t = track('linear',20,-pi/8);
t.initial_position = [60;0;1.5];
t.interpolate_positions( 128/20 );
t.segment_index       = [1,40,90];
t.scenario            = { 'BERLIN_UMa_LOS' , 'BERLIN_UMa_NLOS' ,...
    'BERLIN_UMa_LOS' };
t.interpolate_positions( s.samples_per_meter );

l = layout( s );
l.tx_array.generate('dipole');
l.rx_array = l.tx_array;
l.tx_position(3) = 25;
l.track = t;

l.visualize;

RandStream.setGlobalStream(RandStream('mt19937ar','seed',5));
p = l.create_parameter_sets;


%% Channel generation and results
% Next, we generate the channel coefficients. Note that here, the initial
% sample density is 2.5. We then interpolate the sample density to 20. It
% would take ten times as long to achieve the same result with setting the
% initial sample density to 20. The interpolation is significantly faster.
% It is done by first setting the speed to 1 m/s (default setting) and then
% creating a distance vector which contains a list of effective sampling
% points along the track.      

c = p.get_channels;
cn = c.merge;

t.set_speed( 1 );
dist = t.interpolate_movement( s.wavelength/(2*20) );
ci = cn.interpolate( dist , 'spline' );


%%
% The next plot shows the power of the first three taps from both, the
% original and the interpolated channel, plotted on top of each other. The
% values are identical except for the fact, that the interpolated values
% (blue line) have 5 times as many sample points.

pwr_orig = 10*log10(squeeze(abs(cn.coeff(1,1,1:3,:))).^2);
nsnap = cn.no_snap;
dist_orig = (0:nsnap-1) * t.get_length/(nsnap-1);
pwr_int = 10*log10(squeeze(abs(ci.coeff(1,1,1:3,:))).^2);

figure
plot( dist_orig,pwr_orig , 'r','Linewidth',2 )
hold on
plot( dist,pwr_int ,'b' )
hold off
axis(  [ min(dist) , max(dist) ,...
    min( pwr_orig( pwr_orig>-Inf ) ) , ...
    max( pwr_orig( pwr_orig>-Inf ) )+10 ] )

xlabel('Distance from start point [m]');
ylabel('Power [dB]');


%%
% The following plot shows the power delay profile (PDP) for the
% interpolated channel. As defined in the track object, it starts with a
% LOS segment, going into a shaded area with significantly more multipath
% fading at around 4 seconds and then back to LOS at around 13 sec.

h = ci.fr( 100e6,512 );
h = squeeze(h);
pdp = 10*log10(abs(ifft(h,[],1).').^2);

figure
imagesc(pdp(:,1:256));
caxis([ max(max(pdp))-50 max(max(pdp))-5 ])
colorbar

cm = colormap('hot');
colormap(cm(end:-1:1,:));

set(gca,'XTick',1:32:255);
set(gca,'XTickLabel',(0:32:256)/100e6*1e6);
xlabel('Delay [\mus]');

set(gca,'YTick',1:ci.no_snap/8:ci.no_snap);
set(gca,'YTickLabel', (0:ci.no_snap/8:ci.no_snap)/ci.no_snap * 20 );
ylabel('Time [s]');


%%
% Now, we create a movement profile. It is defined by a set of value pairs
% in track.movement_profile. The first value represents the time in
% seconds, the second value the position on the track. Here, we start at a
% position of 7 m, i.e. in the second (NLOS) segment. We then go back to
% the beginning of the track. This takes 5 seconds. Then, we wait there for
% 1 second and go to the end of the track, which we reach after additional
% 14 seconds.      
%
% The next step is to interpolate the sample points. This is done by the
% interpolate_movement method. It requires the sample interval (in s) as an
% input argument. Here, we choose an interval of 1 ms which gives us 1000
% samples per second. The plot the illustrates the results.

t.movement_profile = [ 0,7 ; 5,0 ; 6,0 ; 20,20  ]';
dist = t.interpolate_movement( 1e-3 );
ci = cn.interpolate( dist );

nsnap = ci.no_snap;
time = (0:nsnap-1) * t.movement_profile(1,end)/(nsnap-1);

figure
plot( time,dist , 'r' )

xlabel('Time [s]');
ylabel('Position on track [m]');


%%
% The last plot shows the PDP of the interpolated channel with the movement
% profile applied. The channel starts in the second segment with a lot of
% fading, goes back to the first while slowing down at the same time. After
% staying constant for one second, the channel starts running again,
% speeding up towards the end of the track.

h = ci.fr( 100e6,512 );
h = squeeze(h);
pdp = 10*log10(abs(ifft(h,[],1).').^2);

figure
imagesc(pdp(:,1:256));
caxis([ max(max(pdp))-50 max(max(pdp))-5 ])
colorbar

cm = colormap('hot');
colormap(cm(end:-1:1,:));

set(gca,'XTick',1:32:255);
set(gca,'XTickLabel',(0:32:256)/100e6*1e6);
xlabel('Delay [\mus]');

set(gca,'YTick',1:ci.no_snap/8:ci.no_snap);
set(gca,'YTickLabel', (0:ci.no_snap/8:ci.no_snap)/ci.no_snap * 20 );
ylabel('Time [s]');

%%
%close all
disp(['QuaDRiGa Version: ',simulation_parameters.version])

%%
% The following code segment shows a movie of the channel response.

% 
% h = ci.fr( 20e6,128 );
% h = squeeze(h);
% 
% mi = -90; ma = -80;
% while true
%     for n = 1:size(h,2)
%         pdp  = 10*log10(abs(h(:,n)).^2);
%         plot(pdp)
%         ma = max( ma,max([pdp]) );
%         mi = min( mi,min([pdp]) );
%         axis([1,128,mi,ma])
%         title(round(time(n)))
%         drawnow
%     end
% end
