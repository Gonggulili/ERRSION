%% 对解调的QAM符号直接做硬判决
function y=realQAMhardJudge(symb,QAMtype)
N=length(symb);
for n=1:1:N
    y(n,1)=subHardJudge(symb(n),QAMtype);
end
end

function y=subHardJudge(data,QAMtype)
if(QAMtype==1||QAMtype==2)
    if(data>0)
        y=1;
    else
        y=-1;
    end
elseif(QAMtype==4)
    if(data>2)
        y=3;
    elseif(data>0&&data<=2)
        y=1;
    elseif(data<=0&&data>-2)
        y=-1;
    else
        y=-3;
    end
elseif(QAMtype==6)
    if(data>=6)
        y=7;
    elseif(data>=4&&data<6)
        y=5;
    elseif(data>=2&&data<4)
        y=3;
    elseif(data>=0&&data<2)
        y=1;
    elseif(data<0&&data>=-2)
        y=-1;
    elseif(data<-2&&data>=-4)
        y=-3;
    elseif(data>=-6&&data<-4)
        y=-5;
    else
        y=-7;
    end
end
end
        