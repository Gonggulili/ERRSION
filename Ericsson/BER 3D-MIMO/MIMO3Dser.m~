function [ser,div_prob]=MIMO3Dser(mode,Monte,snr,Nr,sysParam)
sysParam.Nr=Nr;
% N=sysParam.N;

pusch.Modulation = '16QAM';%----------GY
TrBlkSizes=1980; %------------GY
ue.NCellID = 1;  %------------GY
ue.NSubframe = 0; %------------GY
ue.RNTI = 1; %------------GY
pusch.PRBSet = [0:5].';%------------GY
pusch.RV = 0;%------------GY
 

U=sysParam.U;
QAMtype=sysParam.QAMtype;
noiseVar=10^(-snr/10)*U;
gamParam.maxIteration=20;
gamParam.tol=10^-6;
gamParam.dampSig=0.6;
gamParam.sig_mu=0;
gamParam.sig_var=1;
gamParam.noiseVar=max(noiseVar,0.1);
load(sysParam.chanFile)  %�ṩhmonte�����ʱ���ŵ���Ӧ
chanMonte=size(hmonte,4);
div_event=0;
dic=1;

for monte=1:1:Monte 
%     div=zeros(1,N);
%     for u=1:1:U %----------------wsc
%         q(u,:)=randQAMsymbols(QAMtype,N);
%         t(u,:)=1/sqrt(N)*fft(q(u,:));
%     end
%     
    
    for u=1:1:U  %----------------gy
        trBlk = randi([0,1],TrBlkSizes,1);
        cw = lteULSCH(ue,pusch,trBlk );
        puschSym(u,:) = ltePUSCH(ue,pusch,cw);
        N1 = length(puschSym(u,:));
        t(u,:)=1/sqrt(N1)*fft(puschSym(u,:));
    end
    
     div=zeros(1,N1);%----gy
    
    chan_monte=ceil(rand(1)*chanMonte);
%     Ht=vecBit1FFTh(squeeze(hmonte(:,:,:,chan_monte)),N);
    Ht=vecBit1FFTh(squeeze(hmonte(:,:,:,chan_monte)),N1);
    %���ز����
%     for n=1:1:N 
    for n=1:1:N1 
        H=squeeze(Ht(:,:,n));
        y=H*t(:,n)+wgn(Nr,1,10*log10(noiseVar),'complex');
        if(strcmp(mode,'SC-GAMP'))
            [t_(:,n),div(n)]=SCGAMPdetect(y,H,gamParam,sysParam);
        elseif(strcmp(mode,'SC-MMSE'))
            t_(:,n)=inv(H'*H+noiseVar*eye(U))*H'*y;
        end
    end
    
    if(~sum(div))
        for u=1:1:U
            q_(u,:)=sqrt(N)*ifft(t_(u,:));
            symb_camp(u,:)=compQAMhardDec(q_(u,:),QAMtype);
            ser_(u,dic)=compSERjudge(symb_camp(u,:),q(u,:));
        end
        dic=dic+1;
    end
    
    div_event=div_event+sum(div);

end


ser=mean(mean(ser_,2));
% div_prob=div_event/(Monte*N);
div_prob=div_event/(Monte*N1);
end

function y=compQAMhardDec(xc_gamp,QAMtype)
Nt=length(xc_gamp);
x_gamp=[real(xc_gamp) imag(xc_gamp)];
amp=wlanQAMamplitude(QAMtype);
symb_gamp=realQAMhardJudge(x_gamp/amp,QAMtype);
y=amp*(symb_gamp(1:Nt)+i*symb_gamp(Nt+1:2*Nt));
end

function H=vecBit1FFTh(h,N)
[Nt,U,L]=size(h);
for u=1:1:U
    for nt=1:1:Nt
        vec=[squeeze(h(nt,u,:));zeros(N-L,1)];
        H(nt,u,:)=fft(vec);
    end
end
end

function [encbits, info] = lteULSCH(varargin)
    
    [encbits, info] = mwltelibrary(['lteULSCH' varargin]);
    
    if (iscell(info))
        info = [info{:}];
    end
    
end

