function L=CholDecompose(A)
N=size(A,1);
L=zeros(N,N);
L(1,1)=sqrt(A(1,1));
L(2:N,1)=A(2:N,1)/L(1,1);
for j=2:1:N-1
    L(j,j)=sqrt(A(j,j)-L(j,1:j-1)*L(j,1:j-1)');
    for i=j+1:1:N
        L(i,j)=1/L(j,j)*(A(i,j)-L(i,1:1:j-1)*L(i,1:1:j-1)');
    end
end
L(N,N)=sqrt(A(N,N)-L(N,1:1:N-1)*L(N,1:1:N-1)');
end