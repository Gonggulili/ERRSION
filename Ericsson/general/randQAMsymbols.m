function y=randQAMsymbols(type,symbols)
amp=wlanQAMamplitude(type);
if(type==1)
    symbVec=[-1 1];
elseif(type==2)
    symbVec=[-1 1];
elseif(type==4)
    symbVec=-3:2:3;
elseif(type==6)
    symbVec=-7:2:7;
elseif(type==8)
    symbVec=-15:2:15;
end
L=length(symbVec);
index=ceil(L*rand(1,symbols));
I(:,1)=symbVec(index);
if(strcmp(type,'BPSK'))
    y=amp*I;
else
    index=ceil(L*rand(1,symbols));
    Q(:,1)=symbVec(index);
    y=amp*(I+i*Q);
end