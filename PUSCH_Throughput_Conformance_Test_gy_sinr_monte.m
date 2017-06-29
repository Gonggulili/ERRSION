%% PUSCH Throughput Conformance Test

%% Simulation Configuration
times=1;
NFrames = 1;                % Number of frames to simulate at each SNR
SNRIn = 1; % SNR points to simulate
Length = length(SNRIn); %------6.19
Ber =zeros(1,Length);  %------6.19

%% UE Configuration
% User Equipment (UE) settings are specified in a structure form.

ue.TotSubframes = 1; % Total number of subframes to generate a waveform for
ue.NCellID = 10;     % Cell identity
ue.RC = 'A3-7';      % FRC number  A3-7de ����Ϊ10296
ue.NTxAnts = 4;     % Number of transmit antennas

%% Propagation Channel Model Configuration  ԭ����

chcfg.NRxAnts = 4;               % Number of receive antenna
chcfg.DelayProfile = 'EPA';      % Delay profile
chcfg.DopplerFreq = 5.0;         % Doppler frequency    
chcfg.MIMOCorrelation = 'Low';   % MIMO correlation
chcfg.Seed = 100;                % Channel seed    
chcfg.NTerms = 16;               % Oscillators used in fading model
chcfg.ModelType = 'GMEDS';       % Rayleigh fading model type 
chcfg.InitPhase = 'Random';      % Random initial phases
chcfg.NormalizePathGains = 'On'; % Normalize delay profile power
chcfg.NormalizeTxAnts = 'On';    % Normalize for transmit antennas

%% Channel Estimator Configuration

cec.PilotAverage = 'UserDefined'; % Type of pilot averaging 
cec.FreqWindow = 13;              % Frequency averaging windows in REs
cec.TimeWindow = 1;               % Time averaging windows in REs
cec.InterpType = 'cubic';         % Interpolation type
cec.Reference = 'Antennas';       % Reference for channel estimation
% cec.Reference = 'Layers' ;        %---GY

%% Uplink RMC Configuration

% Generate FRC configuration structure for A3-2
frc = lteRMCUL(ue);
rvSeq = frc.PUSCH.RVSeq;
frc.PUSCH.NLayers = 4;  % the layers are shared across the 1 codewords 

% Transport block sizes for each subframe within a frame
trBlkSizes = frc.PUSCH.TrBlkSizes;
codedTrBlkSizes = frc.PUSCH.CodedTrBlkSizes;

%% Setup HARQ Processes 
% Generate HARQ process table
noHarqProcesses = 8;
harqTable = mod(0:noHarqProcesses-1, noHarqProcesses)+1;  

%% Set Propagation Channel Model Sampling Rate
% The sampling rate for the channel model is set using the value returned
% from <matlab:doc('lteSCFDMAInfo') lteSCFDMAInfo>.

info = lteSCFDMAInfo(frc);
chcfg.SamplingRate = info.SamplingRate;     

%% Processing Loop
% Initialize variables used in the simulation and analysis
totalBLKCRC = zeros(numel(SNRIn), NFrames*10);   % Total block CRC vector
bitThroughput = zeros(numel(SNRIn), NFrames*10); % Total throughput vector
resultIndex = 1;        % Initialize frame counter index

%% �����ŵ���������
% pathname='D:\matlab2016(b)\bin\Ericsson\data\';
U  = ue.NTxAnts; %�û���������ƽ̨��Ϊ����������Ŀ
Nr = chcfg.NRxAnts; %����������Ŀ
N1 = info.Nfft;    %�ʹ���й�
Monte=30;  %������
ant_dist=1.5;%���߼�ľ���
% chan_filename=[pathname 'MIMO3Dchan' num2str(Nr) 'x' num2str(U) 'Wave' num2str(10*ant_dist) '.mat'];
% sysParam.chanFile=chan_filename;
% load(sysParam.chanFile)  %�ṩhmonte�����ʱ���ŵ���Ӧ
% chanMonte=size(hmonte,4);%ȡ�õ�����
% chan_monte=ceil(rand(1)*chanMonte);%��ȡ��

%% ����������������������GAMP��ⲿ��6.15��������������������������������������������������

gamParam.maxIteration=20;%��������
gamParam.tol=10^-6;%����������������
gamParam.dampSig=0.6;%�������ӣ���ǿ�ȶ��Ե�
gamParam.sig_mu=0;%�źŵ�����ֲ���ֵ
gamParam.sig_var=0.0000714;%����
% gamParam.noiseVar=max(noiseVar,0.1);%����������
div_event=0;
dic=1;
sysParam.Nr=chcfg.NRxAnts;
sysParam.U=ue.NTxAnts;


%% 
 for jjjj = 1:times
     fprintf('Monte %d\n',jjjj);
    for iiii = 1:Length
        SNRdB = SNRIn(iiii);

        fprintf('\nSimulating at %g dB SNR for a total %d Frame(s)', ...
            SNRdB, NFrames);

        % Calculate required AWGN channel noise
        SNR = 10^(SNRdB/20);
        N = 1/(SNR*sqrt(double(info.Nfft)))/sqrt(2.0);    
        rng('default');

        % Store results for every subframe at SNR point
        bitTp = zeros(1, NFrames*10);  % Intermediate bit throughput vector	
        blkCRC = zeros(1, NFrames*10); % Intermediate block CRC vector         

        % Initialize state of all HARQ processes
        for i = 1:8
            harqProc(i) = hPUSCHNewHARQProcess( ...
                trBlkSizes(i), codedTrBlkSizes(i), rvSeq); %#ok
        end

        offsetused = 0;
        for subframeNo = 0:(NFrames*10-1)

            % Update subframe number
            frc.NSubframe = subframeNo;

            % Get HARQ index for given subframe from HARQ index table
            harqIdx = harqTable(mod(subframeNo, length(harqTable))+1);   

            % Update current HARQ process
            harqProc(harqIdx) = hPUSCHHARQScheduling(harqProc(harqIdx));
            frc.PUSCH.RV = harqProc(harqIdx).rvSeq(harqProc(harqIdx).rvIdx);
            frc.PUSCH.RVSeq = harqProc(harqIdx).rvSeq(harqProc(harqIdx).rvIdx);


            bbb= harqProc(harqIdx).ulschTransportBlk;

            % Create an SC-FDMA modulated waveform
            [txWaveform, txSubframe] = lteRMCULTool( ...
                frc, harqProc(harqIdx).ulschTransportBlk);

            % Transmit an additional 25 samples at the end of the waveform to
            % cover the range of delays expected from the channel modeling
             txWaveform = [txWaveform; zeros(25, ue.NTxAnts)]; %#ok

            % The initialization time for channel modeling is set each subframe
            % to simulate a continuously varying channel
            chcfg.InitTime = subframeNo/1000;

            % Pass data through channel model  -------ԭ����
            rxWaveform = lteFadingChannel(chcfg,txWaveform);  % ----6.6

            % Add noise at the receiver
            v = N*complex(randn(size(rxWaveform)), randn(size(rxWaveform)));
            rxWaveform = rxWaveform+v;

            %--------�����ŵ������Ӧ��ɣ�6.6��------------

    %           Ht=vecBit1FFTh(squeeze(hmonte(:,:,:,chan_monte)),N1);%����ΪɶҪ���Σ��˴���NΪ�����ݴ���С����
    %           
    %           for n=1:1:N1
    %                 H=squeeze(Ht(:,:,n));%��ά��������ı任Ϊ����
    %                 txWaveform_909=txWaveform(:,n);
    %                 rxWaveform=conv(H,txWaveform(:,n))+v;
    %                 
    %                 
                    %yΪ���ն��յ����źţ������ɷ����źž���ŵ��弤��ӦH+��������ģ�WGN�����������
    % 
    %                 if(strcmp(mode,'SC-GAMP'))
    %                     [t_(:,n),div(n)]=SCGAMPdetect(y,H,gamParam,sysParam);%��һ�δ������Ϊ��ͨ��y��H�ָ�����t(:,n)��������㷨��SCGAMPdetect
    %                   elseif(strcmp(mode,'SC-MMSE'))%��һ�δ���Ҳ��Ϊ��ͨ��y��H�ָ�����t(:,n)������������һ���㷨�У�SCMMSEdetect
    %                     t_(:,n)=inv(H'*H+noiseVar*eye(U))*H'*y;
    %                 end
    %           end


            % Calculate synchronization offset
            offset = lteULFrameOffset(frc, frc.PUSCH, rxWaveform);
            if (offset < 25)
                offsetused = offset;
            end

            % SC-FDMA demodulation
            rxSubframe = lteSCFDMADemodulate(frc, ...
                rxWaveform(1+offsetused:end, :));

            % Channel and noise power spectral density estimation
            [estChannelGrid, noiseEst] = lteULChannelEstimate(frc, ... 
                frc.PUSCH, cec, rxSubframe);
            gamParam.noiseVar=max(noiseEst,0.1);%-----by GY

            % Extract REs corresponding to the PUSCH from the given subframe
            % across all receive antennas and channel estimates
            puschIndices = ltePUSCHIndices(frc, frc.PUSCH);
            [puschRx, puschEstCh] = lteExtractResources( ...
                puschIndices, rxSubframe, estChannelGrid);

            % MMSE equalization
             rxSymbols = lteEqualizeMMSE(puschRx, puschEstCh, noiseEst);

            % GAMP equalization---------------6.15����
             [rxSymbols1,div]=SCGAMPdetect(puschRx, puschEstCh,gamParam,sysParam);

             %-------------------------------6.16����
    %           if(~sum(div)) %�߼���Ϊ�棬��div=0ʱ�����������GAMP�㷨����ʱ��mmse�㷨��Զ����
    %             for u=1:1:U
    % %                 q_(u,:)=sqrt(N)*ifft(t_(u,:)); %ʱ��תƵ���渵��Ҷ�仯��ΪʲôҪ����sqrt(N)
    %                 symb_camp(u,:)=compQAMhardDec(q_(u,:),QAMtype);
    %                 ser_(u,dic)=compSERjudge(symb_camp(u,:),q(u,:));
    %             end
    %             dic=dic+1;
    %           end
    %          div_event=div_event+sum(div);
             %-----------------------------------------------------------------

            % Update frc.PUSCH to carry complete information of the UL-SCH
            % coding configuration
            frc.PUSCH = lteULSCHInfo(frc, ...
                frc.PUSCH, harqProc(harqIdx).trBlkSize, 'chsconcat');

            % Decode the PUSCH
            rxEncodedBits = ltePUSCHDecode(frc, frc.PUSCH, rxSymbols1);

            % Decode the UL-SCH channel and store the block CRC error for given
            % HARQ process harqIdx
            trBlkSize = trBlkSizes(mod(subframeNo, 10)+1);
            [rxDecodedBits, harqProc(harqIdx).crc, ...
                harqProc(harqIdx).decState] = lteULSCHDecode(...
                frc, frc.PUSCH, trBlkSize, ...
                rxEncodedBits, harqProc(harqIdx).decState);

            % Store the CRC calculation and total number of bits per subframe
            % successfully decoded
            blkCRC(subframeNo+1) = harqProc(harqIdx).crc;
            bitTp(subframeNo+1) = ...
                harqProc(harqIdx).trBlkSize.*(1-harqProc(harqIdx).crc);

        end     

        % Record the block CRC error and bit throughput for the total number of
        % frames simulated at a particular SNR
        totalBLKCRC(resultIndex, :) = blkCRC;
        bitThroughput(resultIndex, :) = bitTp;
        resultIndex = resultIndex + 1;
        %-------------------------��������BER-----6.16-------------------------
            numErrs = sum(rxDecodedBits ~= bbb) ;                                       
            numBits = length(rxDecodedBits) ;
            Ber(iiii) = numErrs/numBits ;



    end 
 end

 ber=mean(Ber);



%% Display Throughput Results
% Throughput calculation as a percentage
% throughput = 100*(1-mean(totalBLKCRC, 2)).';

% hPUSCHResults(SNRIn, NFrames, trBlkSizes, throughput, bitThroughput);

%% ------------------����BER����6.16-----------------------------
figure;
semilogy(SNRIn,Ber,'r.-');
hold on
axis([SNRIn(1) SNRIn(end) 10^-5 1])
grid on
xlabel('SNR (dB)')
ylabel('ber')
title(' BER');

% isequal(codeword,awcw>0)


%% Appendix
% This example uses the helper functions:
%
% * <matlab:edit('hPUSCHHARQScheduling.m') hPUSCHHARQScheduling.m>
% * <matlab:edit('hPUSCHNewHARQProcess.m') hPUSCHNewHARQProcess.m>
% * <matlab:edit('hPUSCHResults.m') hPUSCHResults.m>

%% Selected Bibliography
% # 3GPP TS 36.104 "Base Station (BS) radio transmission and reception"

displayEndOfDemoMessage(mfilename)



% -----------����6.14----------------------
function H=vecBit1FFTh(h,N)
   [Nt,U,L]=size(h);
  for u=1:1:U
     for nt=1:1:Nt
%        vec=[squeeze(h(nt,u,:));zeros(N-L,1)];          %-----ԭ����
%          H(nt,u,:)=fft(vec);%ʱ��תƵ�򣬽��и���Ҷ�任  %-----ԭ����
        H(nt,u,:)=[squeeze(h(nt,u,:));zeros(N-L,1)];%-----GY
     end
  end
end
%-------------����6.15-----------------------
function [x,div]=SCGAMPdetect(reGrid, Hest,gamParam,sysParam)



%---------����
    if (length(size(Hest))==4 || length(size(reGrid))==3 ||...
        size(Hest,2)>8 || size(reGrid,2)>8)
        fullgrid=1;
      else
        fullgrid=0;
        Hest=reshape(Hest,1,size(Hest,1),size(Hest,2),size(Hest,3));
        reGrid=reshape(reGrid,1,size(reGrid,1),size(reGrid,2));
    end

    noSc = size(reGrid,1);
    noSymb = size(reGrid,2);
    noTxAnt = size(Hest,4);    
    Hest = permute(Hest,[3 4 1 2]);   % Permute to noRxAnts-by-noTxAnts-... 
    reGrid = permute(reGrid,[3 1 2]); % to allow quicker access
    x = zeros(noSc,noSymb,noTxAnt);
%     csi = zeros(noSc,noSymb,noTxAnt);
%     n0eye = n0*eye(noTxAnt);
%--------------------------------------------------------------
for scNo=1:noSc
        for symbNo=1:noSymb
            Nr=sysParam.Nr;
            U=sysParam.U;
            maxIteration=gamParam.maxIteration;
            tol=gamParam.tol;
            noiseVar=gamParam.noiseVar;
            
            dampSig=gamParam.dampSig;
            sig_mu=gamParam.sig_mu;
            sig_var=gamParam.sig_var;

            xhat=wgn(U,1,10*log10(1),'complex');%��ʼֵ��1
            vx=eye(U);%�źŵķ���
            vx=0.0000714*eye(U);%�źŵķ���----by  GY
            shat=zeros(Nr,1);

            iter=0;
            stop=0;
            div=0;
            
            H = Hest(:,:,scNo,symbNo);
%             T = H'*H;
%             C = inv(T+n0eye);
%             csi(scNo,symbNo,:) = 1./real(diag(C));
%             G = C*H'; %#ok<MINV>
            y = reGrid(:,scNo,symbNo);           
%             out(scNo,symbNo,:) = G*r;


            while(iter<maxIteration&&(~stop))
                for nr=1:1:Nr
                     vp(nr,1)=abs(H(nr,:)*vx*H(nr,:)'); %���abs()?�ƺ��Ƶ�û����˳
            %          vp(nr,1)=(H(nr,:)*vx*H(nr,:)');
%                       vp(nr,1)=abs(H(nr,:).*vx.*H(nr,:)'); %�ڲ�����ά�ȱ���һ�£�BY GY
                end
                phat=H*xhat-vp.*shat;
                [shat,vs]=shat_vs_update(y,phat,vp+noiseVar);%���ӽڵ�ľ�ֵ����

                vr=update_rCov(H,vs,U);
                rhat=xhat+vr*(H'*shat); %BUG!H'ǰ�����conj��ȥ������Գɹ�
                [xhat,vx]=fa_fv(rhat,vr,sig_mu,sig_var,U);
                if(iter)
                    xhat=dampSig*xhat_prev+(1-dampSig)*xhat;%ǰһ�εĽ�͵�ǰ��
                    vx=dampSig*vx_prev+(1-dampSig)*vx;
                    if(norm(xhat-xhat_prev)<tol)
                        stop=1;
                    end
                end
                if(isnan(sum(xhat))) %�ж��Ƿ�ɢ
                    div=1;
                    stop=1;
                end
                iter=iter+1;
                xhat_prev=xhat;
                vx_prev=vx;
            end
            x=xhat;
            x1(scNo,symbNo,:)=x;     %---by  GY
            
  
    end
end  
          if (~fullgrid)
              x = reshape(x1(1,:,:),noSymb,noTxAnt);
          end

end

function [shat,vs]=shat_vs_update(y,phat,vp)
shat=(y-phat)./vp;
vs=1./vp;
end


function vr=update_rCov(H,vs,U)
Nr=size(H,1);
for r=1:1:Nr
    Hs(r,:)=vs(r)*H(r,:);
end
for j=1:1:U
    Hj=H(:,j);
    Hsj=Hs(:,j);
    vr(j,j)=inv(Hsj'*Hj); %BUG! miss inv()
%     vr(index,index)=inv(Hj'*Hsj);
end
end

function [xhat,vx]=fa_fv(rhat,vr,sig_mu,sig_var,U)
for j=1:1:U
    [xhat(j,1),vx(j,j)]=sub_fa_fv(rhat(j),vr(j,j),sig_mu,sig_var);
end
end

function [xhat,vx]=sub_fa_fv(rhat,vr,sig_mu,sig_var)
vx=1./(1./vr+1./sig_var);
xhat=(rhat./vr+sig_mu/sig_var).*vx;
end



