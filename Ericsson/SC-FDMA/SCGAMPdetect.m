%% 核心代码。只有相位观测(p)下GAMP信号检测 
function [x,div]=SCGAMPdetect(y,H,gamParam,sysParam)
Nr=sysParam.Nr;
U=sysParam.U;
maxIteration=gamParam.maxIteration;
tol=gamParam.tol;
noiseVar=gamParam.noiseVar;
dampSig=gamParam.dampSig;
sig_mu=gamParam.sig_mu;
sig_var=gamParam.sig_var;

xhat=wgn(U,1,10*log10(1),'complex');
vx=eye(U);
shat=zeros(Nr,1);

iter=0;
stop=0;
div=0;
while(iter<maxIteration&&(~stop))
    for nr=1:1:Nr
        vp(nr,1)=abs(H(nr,:)*vx*H(nr,:)'); %需加abs()?似乎推导没有理顺
%         vp(nr,1)=(H(nr,:)*vx*H(nr,:)');
    end
    phat=H*xhat-vp.*shat;
    [shat,vs]=shat_vs_update(y,phat,vp+noiseVar);

    vr=update_rCov(H,vs,U);
    rhat=xhat+vr*(H'*shat); %BUG!H'前多加了conj；去掉后调试成功
    [xhat,vx]=fa_fv(rhat,vr,sig_mu,sig_var,U);
    if(iter)
        xhat=dampSig*xhat_prev+(1-dampSig)*xhat;
        vx=dampSig*vx_prev+(1-dampSig)*vx;
        if(norm(xhat-xhat_prev)<tol)
            stop=1;
        end
    end
    if(isnan(sum(xhat))) %判断是否发散
        div=1;
        stop=1;
    end
    iter=iter+1;
    xhat_prev=xhat;
    vx_prev=vx;
end 
x=xhat;
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



