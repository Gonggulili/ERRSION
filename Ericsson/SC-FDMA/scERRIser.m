function ser=scERRIser(mode,Monte,snr,Nr,sysParam)
sysParam.Nr=Nr;
N=sysParam.N;
U=sysParam.U;
QAMtype=sysParam.QAMtype;
noiseVar=10^(-snr/10)*U;

gamParam.maxIteration=35;
gamParam.tol=10^-6;
gamParam.dampSig=0.3;
gamParam.sig_mu=0;
gamParam.sig_var=1;
gamParam.noiseVar=noiseVar;
    
for monte=1:1:Monte
    for u=1:1:U
        q(u,:)=randQAMsymbols(QAMtype,N);
        t(u,:)=1/sqrt(N)*fft(q(u,:));
    end
    %ÖðÔØ²¨¼ì²â
    for n=1:1:N 
        H=wgn(Nr,U,0,'complex');
        y=H*t(:,n)+wgn(Nr,1,10*log10(noiseVar),'complex');
        if(strcmp(mode,'SC-GAMP'))
            t_(:,n)=SCGAMPdetect(y,H,gamParam,sysParam);
        elseif(strcmp(mode,'SC-MMSE'))
            t_(:,n)=inv(H'*H+noiseVar*eye(U))*H'*y;
        end
    end
    for u=1:1:U
        q_(u,:)=sqrt(N)*ifft(t_(u,:));
        symb_camp(u,:)=compQAMhardDec(q_(u,:),QAMtype);
        ser_(u,monte)=compSERjudge(symb_camp(u,:),q(u,:));
    end
end
ser=mean(mean(ser_,2));
end

function y=compQAMhardDec(xc_gamp,QAMtype)
Nt=length(xc_gamp);
x_gamp=[real(xc_gamp) imag(xc_gamp)];
amp=wlanQAMamplitude(QAMtype);
symb_gamp=realQAMhardJudge(x_gamp/amp,QAMtype);
y=amp*(symb_gamp(1:Nt)+i*symb_gamp(Nt+1:2*Nt));
end




