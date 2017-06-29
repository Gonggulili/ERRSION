%% Effects of the Antenna-Orientation
%
% This tutorial shows how to evaluate antenna effects. It creates a simple
% setup with a transmit and a receive antenna facing each other in pure LOS
% conditions. Then, the transmitter is rotated around its x-axis and the
% effect on the received power is studied.
%
% One feature of the model is that it allows to freely orient the antennas
% at the transmitter and receiver. In the following, two cross-polarized
% patch antennas are aligned on the optical axis facing each other. The
% surface normal vectors of the transmit and the receive patch are aligned
% with the LOS. The transmitter is rotated from -90° to 90° around the
% optical axis. The real and imaginary parts of the channel coefficients
% are then calculated for each angle. Each real and imaginary part is
% normalized by its maximum and the results are plotted. The calculation is
% done for both, linearly and crossed polarized elements.

%% Model and Antenna Setup
% Here, we parametrize the simulation. We place the receiver 10 m away from
% the transmitter and chose the scenario "LOSonly". Thus, no NLOS
% components are present. The receiver is set up as a multi-element array
% using both, circular and linear polarization.

clear all
close all

s = simulation_parameters;
s.sample_density = 2;
s.show_progress_bars = 0;

l = layout(s);
l.track.generate('linear',0,0);
l.track.scenario = 'LOSonly';
l.rx_position = [11;0;0];
l.tx_position = [0;0;0];

a = array('lhcp-rhcp-dipole');
a.generate('custom',3,90,90,0);
a.set_grid( (-180:10:180)*pi/180 , (-90:10:90)*pi/180 );
a.rotate_pattern(180,'z',3);
a.copy_element(3,4);
a.rotate_pattern(90,'x',4);

b = a.copy_objects;
a.rotate_pattern(180,'z');

l.tx_array = a;
l.rx_array = b;


%% Iteration over all angles
% Next, we rotate the receive antenna in 10 degree steps around its x-axis
% and calculate the channel  response for each angle.

[~,~,cb] = l.get_channels;      % Get the channel builder object

rot = -120:10:120;
h = zeros(4,4,numel(rot));
for n=1:numel(rot)
    
    cc = a.copy_objects;
    cc.rotate_pattern( rot(n) , 'x');
    
    cb.par.tx_array = cc;
    c = cb.get_channels;
    
    h(:,:,n) = c.coeff(:,:,1,1);
end

%% Linear Polarization results
% Now we plot the results for the linear polarization. There are two
% linearly polarized antennas at the transmitter and two at the receiver.
% Their orientation can either be vertical (denoted as V) or horizontal
% (denoted as H). The channel matrix thus has 8 coefficients, VV, VH, HV
% and HH. Each coefficient is complex-valued. Thus, figure shows 8 curves,
% 4 for the real parts and 4 for the imaginary parts.

set(0,'defaultTextFontSize', 16)
set(0,'defaultAxesFontSize', 16)
set(0,'DefaultAxesFontName','Times')
set(0,'DefaultTextFontName','Times')

figure('Position',[ 100 , 100 , 760 , 400]);
g = h([3,4],[3,4],:);

plot(rot,squeeze(real(g(1,1,:))),'-sk','Linewidth',0.5,'MarkerfaceColor','k','Markersize',12)
hold on
plot(rot,squeeze(real(g(2,2,:))),'-db','Linewidth',0.5,'MarkerfaceColor','b','Markersize',8)
plot(rot,squeeze(real(g(2,1,:))),'-or','Linewidth',0.5,'MarkerfaceColor','r','Markersize',8)
plot(rot,squeeze(real(g(1,2,:))),'-^g','Linewidth',0.5,'MarkerfaceColor','g','Markersize',8)

plot(rot,squeeze(imag(g(1,1,:))),'--sk','Linewidth',0.5,'Markersize',12)
plot(rot,squeeze(imag(g(2,2,:))),'--db','Linewidth',0.5,'Markersize',8)
plot(rot,squeeze(imag(g(2,1,:))),'--or','Linewidth',0.5,'Markersize',12)
plot(rot,squeeze(imag(g(1,2,:))),'--^g','Linewidth',0.5,'Markersize',12)
hold off

xlabel('Rotaion Angle')
ylabel('Normalized Amplitude')
legend('real V-V','real H-H','real H-V','real V-H',...
    'imag V-V','imag H-H','imag H-V','imag V-H','location','EastOutside')


%% Circular Polarization results
% The second plot shows the same for circular polarization. The first
% element is LHCP (denoted as L) and the second is RHCP (denoted as R). As
% expected, all cross-polarization coefficients (RL and LR) are zero.

figure('Position',[ 100 , 100 , 760 , 400]);
g = h([1,2],[1,2],:);

plot(rot,squeeze(real(g(1,1,:))),'-sk','Linewidth',0.5,'MarkerfaceColor','k','Markersize',12)
hold on
plot(rot,squeeze(real(g(2,2,:))),'-db','Linewidth',0.5,'MarkerfaceColor','b','Markersize',8)
plot(rot,squeeze(real(g(2,1,:))),'-or','Linewidth',0.5,'MarkerfaceColor','r','Markersize',8)
plot(rot,squeeze(real(g(1,2,:))),'-^g','Linewidth',0.5,'MarkerfaceColor','g','Markersize',8)

plot(rot,squeeze(imag(g(1,1,:))),'--sk','Linewidth',0.5,'Markersize',12)
plot(rot,squeeze(imag(g(2,2,:))),'--db','Linewidth',0.5,'Markersize',8)
plot(rot,squeeze(imag(g(2,1,:))),'--or','Linewidth',0.5,'Markersize',12)
plot(rot,squeeze(imag(g(1,2,:))),'--^g','Linewidth',0.5,'Markersize',12)
hold off

xlabel('Rotaion Angle')
ylabel('Normalized Amplitude')
legend('real L-L','real R-R','real R-L','real L-R',...
    'imag L-L','imag R-R','imag R-L','imag L-R','location','EastOutside')

%%
%close all
disp(['QuaDRiGa Version: ',simulation_parameters.version])