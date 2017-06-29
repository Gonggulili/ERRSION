function H=MIMO3DchanGen(l,s,sysParam)
U=sysParam.U;
for u=1:1:U
    l.track(u).generate('linear', 60 ,3*pi/4 ); % track of Rx01
    l.track(u).initial_position = [20*rand(2,1);1*rand(1)]; %�û���ʼλ��          % initial position of Rx01
    % l.track(1).scenario = {'WINNER_UMa_C2_NLOS';'WINNER_Indoor_A1_LOS'};   % Two Scenarios of Rx01 for two Tx
    l.track(u).scenario = {'WINNER_UMa_C2_NLOS'}; 
end
% l.track(2).generate('linear', 40,0 );     % track of Rx02
% l.track(2).initial_position = [ 100;50;0 ];        % initial position of Rx02
% l.track(2).scenario = {'WINNER_UMa_C2_LOS';'WINNER_UMa_C2_NLOS'};   % Two Scenarios of Rx02 for two Tx

l.track.interpolate_positions( s.samples_per_meter );

% l.visualize;                                    % Plot all tracks
% view(-33,45)

% Generate parameter sets and channel coefficients
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));
p = l.create_parameter_sets;         % Create parameter sets
c = p.get_channels;             % generate channel coefficients

% function description: generate discrete channel coefficients, the coefficients is stored in [discrete channel variable].coeff 
%                       in the shape of U*S*N*T matrix.
%                         U: number of receiver elements
%                         S: number of transmitter elements
%                         N: number of paths 
%                         T: number of time samples
% function usage:  get_discrete_channel(channel,fc,sample_max)
%                    channel: channel object generate from QuaDRiGa
%                    fc: center_frequency
%                    [sample_max]: optional, the maximum number of snapshot(time sample)
%                                  used here
%
discrete_c = get_discrete_channel(c,s.center_frequency,1);   
for u=1:1:U
    H(:,u,:)=squeeze(discrete_c(u).coeff);
end
[Nr,U,L]=size(H);
for nr=1:1:Nr
    for u=1:1:U
        h=squeeze(H(nr,u,:));
        H(nr,u,:)=h/norm(h);
    end
end




