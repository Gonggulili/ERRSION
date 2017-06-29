function y=wlanQAMamplitude(QAMtype)
if(QAMtype==1)
    y=1;
elseif(QAMtype==2)
    y=1/sqrt(2);
elseif(QAMtype==4)
    y=1/sqrt(10);
elseif(QAMtype==6)
    y=1/sqrt(42);
elseif(QAMtype==8)
    y=1/sqrt(170);
end
    