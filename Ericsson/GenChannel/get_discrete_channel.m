function [channel] = get_discrete_channel(varargin)

tt_max = 0;          % the max number of snapshot (time samples) used in this function

ni=length(varargin);
if ni>0, if (~isempty(varargin{1})), channel=varargin{1}; end, end
if ni>1, if (~isempty(varargin{2})), fc=varargin{2}; end, end
if ni>2, if (~isempty(varargin{3})), tt_max=varargin{3}; end, end
if ni>3, error('Too many input arguments of function get_discrete_channel()!'), end

K = length(channel);    % number of links


%% H is output from wim.m ,the channel matrix. 
%H is a cell array of size K(number of links).Each element of this cell
%array contains a U*S*N*T matrix.
%% delays is a K*N vector of path delay values,output from wim.m
%   U: number of receiver elements
%   S: number of transmitter elements
%   N: number of paths 
%   T: number of time samples
%% 
N=0;M=50;  % M:length of the sequence y(m), N:length of the sequence x(n)
L=-N:M; %  l=m-n m=0:M  n=0:N
W=30.72e6;  %  bandwidth
hl=cell(1,K);   % baseband channel matrix, include K cells, each cell contains a U*S*L*T



numU = zeros(1,K);
numS = zeros(1,K);
numN = zeros(1,K);
numT = zeros(1,K);
delays = cell(1,K);
for numlink=1:K
    numU(1,numlink)=channel(1,numlink).no_rx;
    numS(1,numlink)=channel(1,numlink).no_tx;
    numN(1,numlink)=channel(1,numlink).no_path;
    if tt_max==0
        numT(1,numlink)=channel(1,numlink).no_snap;
    elseif tt_max>0
        numT(1,numlink)=tt_max;
    else
        error('time samples must be greater than zero');
    end
    delays{1,numlink}=channel(1,numlink).delay;
end


%% change a(t) to baseband equivalent a(t)
H_baseband=cell(1,K);
for numlink=1:K
    for nn=1:numN(:,numlink)
        for tt = 1:numT(:,numlink)
            H_baseband{1,numlink}(:,:,nn,tt)=channel(1,numlink).coeff(:,:,nn,tt)*exp(-1j*2*pi*fc*delays{1,numlink}(nn,tt));  
        end
    end
end

    
   
%% calculate the baseband equivalent hl
 for numlink=1:K
     hl_tmp=zeros(numU(1,numlink),numS(1,numlink),length(L),numT(1,numlink));
     for uu=1:numU(1,numlink)
         for ss=1:numS(1,numlink)
             sum=zeros(length(L),numT(1,numlink));
             for tt=1:numT(1,numlink)
                 for l=-N:M                    
                     for nn=1:numN(1,numlink)
                         a=H_baseband{1,numlink}(uu,ss,nn,tt);
                         sum(l+ N+1,tt)=sum(l+N+1,tt)+a*sinc(l-W*delays{1,numlink}(nn,tt));
                     end
                 end
             end
             hl_tmp(uu,ss,:,:)=sum;
         end
     end
     channel(1,numlink).coeff=hl_tmp;
 end


                
