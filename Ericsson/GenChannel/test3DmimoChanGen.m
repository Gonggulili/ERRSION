% First, we set up the simulation environment.
clc
close all
clear all

Monte=30;
Mtx=4;
Ntx=4;
ant_dist=3;%���߼��;���ٸ�����
ant_pol=1;
U=4;

sysParam.U=U;
set(0,'defaultTextFontSize', 14)    % Set default font size for the plots
set(0,'defaultAxesFontSize', 14)

s = simulation_parameters;          % Basic simulation parameters
s.center_frequency = 2.53e9;        % Center Frequency
s.sample_density = 2.5;               % 4 samples / half wave length

% Setting up the antenna arrays
l=layout(s);                        % create a new layout
l.tx_array.generate('3gpp-3d',0,Mtx,Ntx,s.center_frequency,ant_pol,20,ant_dist)  % create 3gpp-3d 4*4 antenna array for Tx 
% l.rx_array.generate('3gpp-3d',0,2,2)
l.rx_array.generate('dipole')       % create dipole antenna for Rx

% Setting up the Tx and Rx position and track
l.no_tx = 1;                        % number of Tx
l.no_rx = U;                        % number of Rx
l.tx_position(:,1) = [10;20;25];    % position of Tx01
% l.tx_position(:,2) = [80;90;30];    % position of Tx02

par=1;
if(par)
    Threads=8;
    distcomp.feature( 'LocalUseMpiexec', false )
    parpool(Threads)
end

parfor monte=1:Monte
    fprintf('Monte %d\n',monte);
    hmonte(:,:,:,monte)=MIMO3DchanGen(l,s,sysParam);
end
% if(par)
%     parpool close
% end

pathname='/home/ubuntu/gongyi/Ericsson/data/';
filename=[pathname 'MIMO3Dchan' num2str(Mtx*Ntx) 'x' num2str(U) 'Wave' num2str(10*ant_dist) '.mat'];
save(filename,'hmonte'); 





 