function y=compSERjudge(QAMhard,QAMsymbol)
err=0;
Users=length(QAMsymbol);
for user=1:1:Users
    symb=QAMhard(user);
    if(real(symb)~=real(QAMsymbol(user))||imag(symb)~=imag(QAMsymbol(user)))
        err=err+1;
    end
end
y=err/Users;