clc
clear all
close all

snrStep=2;
snrVec=-10:snrStep:15;  
QAMtype=4;
U=8;
N=64;
sysParam.QAMtype=QAMtype;
sysParam.U=U;
sysParam.N=N;

%% QPSK
Monte=100;
par=1;
if(par)
    Threads=2;
    distcomp.feature( 'LocalUseMpiexec', false )
    parpool(Threads)
end

parfor t=1:length(snrVec)
    fprintf('第%d个snr情形\n',t);
    snr=snrVec(t);
    Nr=16;
    mode='SC-GAMP';
    ser_gamp1(t)=scERRIser(mode,Monte,snr,Nr,sysParam);
    mode='SC-MMSE';
    ser_mmse1(t)=scERRIser(mode,Monte,snr,Nr,sysParam);
    
    Nr=64;
    mode='SC-GAMP';
    ser_gamp2(t)=scERRIser(mode,Monte,snr,Nr,sysParam);
    mode='SC-MMSE';
    ser_mmse2(t)=scERRIser(mode,Monte,snr,Nr,sysParam);
end
if(par)
   delete(gcp('nocreate'))
end

pathname='D:\matlab2016(b)\bin\Ericsson\data\';
filename='scSER.mat';
matfile=[pathname filename];
save(matfile,'snrVec','ser_gamp1','ser_gamp2','ser_mmse1','ser_mmse2');

snrVec=snrVec+9;
figure
semilogy(snrVec,ser_gamp1,'r.-')
hold on
semilogy(snrVec,ser_gamp2,'r*-')
hold on

semilogy(snrVec,ser_mmse1,'k.-')
hold on
semilogy(snrVec,ser_mmse2,'k*-')
hold on
axis([snrVec(1) snrVec(end) 10^-5 1])
grid on
legend('16x8 GAMP','64x8 GAMP','16x8 MMSE','64x8 MMSE')
xlabel('SNR dB')
ylabel('SER')


