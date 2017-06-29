function L=outCholDecompose(A)
N=size(A,1);
for k=1:1:N
    A(k,k)=sqrt(A(k,k));
    A(k+1:N,k)=A(k+1:1:N,k)/A(k,k);
    for j=k+1:1:N
        A(j:1:N,j)=A(j:1:N,j)-(A(j:N,k))*conj(A(j,k));
    end
end
L=zeros(N,N);
for n=1:1:N
    L(n,1:n)=A(n,1:n);
end