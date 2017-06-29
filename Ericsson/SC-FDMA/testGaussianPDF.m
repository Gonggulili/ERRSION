clc
clear all
close all

Monte=20000;
N=24;
QAMtype=4;

for monte=1:1:Monte
    q=randQAMsymbols(QAMtype,N);
    t=1/sqrt(N)*fft(q);
    s((monte-1)*N+1:monte*N,1)=t;
end

figure
subplot(121)
hist(real(s),100)
xlabel('Real Parts')
subplot(122)
hist(imag(s),100)
xlabel('Imaginary Parts')
title('N=24')






