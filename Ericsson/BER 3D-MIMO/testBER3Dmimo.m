clc
clear all
close all

snrStep=1;
snrVec=0:snrStep:20;  
QAMtype=4;
U=8;
% N=64;
sysParam.QAMtype=QAMtype;
sysParam.U=U;
% sysParam.N=N;
ant_dist=1.5;

pathname='/home/ubuntu/gongyi/Ericsson/data/';

%% QPSK
Monte=1;
par=1;%Matlab2016���ϰ汾�����ô˴���
if(par)
    Threads=8;
    distcomp.feature( 'LocalUseMpiexec', false )
    parpool(Threads)
end

% par=1;
% if(par)
%     Threads=8;
%     distcomp.feature( 'LocalUseMpiexec', false )
%     matlabpool(Threads)
% end

Nr=32;
chan_filename=[pathname 'MIMO3Dchan' num2str(Nr) 'x' num2str(U) 'Wave' num2str(10*ant_dist) '_10x1' '.mat'];
sysParam.chanFile=chan_filename;
tic
for t=1:length(snrVec)
    fprintf('��%d��snr����\n',t);
    snr=snrVec(t);
    mode='SC-GAMP';
    [ser_gamp(t),div_prob(t)]=MIMO3Dser(mode,Monte,snr,Nr,sysParam);%�ڶ������ֵΪ��ɢֵ����Ҫ��η���õ����
    mode='SC-MMSE';
    ser_mmse(t)=MIMO3Dser(mode,Monte,snr,Nr,sysParam);
end
toc
% if(par)
%   delete(gcp('nocreate'))
% end

%����һֱ���?�������ε�
% if(par)
% %     matlabpool close
%     parpool close
% end

pathname='/home/ubuntu/gongyi/Ericsson/data/';
filename='scSER.mat';
matfile=[pathname filename];
save(matfile,'snrVec','ser_gamp','ser_gamp');

figure
semilogy(snrVec,ser_gamp,'r.-')
hold on
semilogy(snrVec,ser_mmse,'k.-')
hold on
axis([snrVec(1) snrVec(end) 10^-5 1])
grid on
legend('64x8 GAMP','64x8 MMSE')
xlabel('SNR dB')
ylabel('SER')


