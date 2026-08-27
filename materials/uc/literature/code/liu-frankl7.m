Dim=27;%90;%60; %max dimension
A=zeros(2*Dim,Dim,Dim);
B=zeros(Dim,Dim,Dim);
C=zeros(Dim,Dim);
D=zeros(Dim,Dim);
O=zeros(Dim+2,Dim+2); 
for m=1:Dim
    for n=1:Dim
        for k=max(m,n):m+n
            A(k,m,n)=(-1)^(m+n-k)*nchoosek(k,k-n)*nchoosek(n,k-m);
        end
    end
end
for m=1:Dim
    for n=1:Dim
        for k=ceil(max(m,n)/2):min(m,n)
            for d=max(m+n-3*k,0):min(m,n)-k
                B(k,m,n)=B(k,m,n)+(-1)^(m+n)*nchoosek(k,d+3*k-m-n)*nchoosek(-2*k-d+m+n,m-k-d)*nchoosek(n-k,d);
            end
        end
    end
end
for m=1:Dim
    for n=1:Dim
        for k=max(m,n):m+n
            C(m,n)=C(m,n)-A(k,m,n)/k/2^k;
        end
    end
end
for m=1:Dim
    for n=1:Dim
        for k=ceil(max(m,n)/2):min(m,n)
            D(m,n)=D(m,n)-B(k,m,n)/k;
        end
    end
end
O(2:Dim+1,2:Dim+1)=2*C;
O(3:Dim+2,2:Dim+1)=O(3:Dim+2,2:Dim+1)-C;
O(2:Dim+1,3:Dim+2)=O(2:Dim+1,3:Dim+2)-C;
O(3:Dim+2,3:Dim+2)=O(3:Dim+2,3:Dim+2)+C;

O(1:Dim,1:Dim)=O(1:Dim,1:Dim)+D;
O(2:Dim+1,2:Dim+1)=O(2:Dim+1,2:Dim+1)-2*D;
O(3:Dim+2,2:Dim+1)=O(3:Dim+2,2:Dim+1)+D;
O(2:Dim+1,3:Dim+2)=O(2:Dim+1,3:Dim+2)+D;
O(3:Dim+2,3:Dim+2)=O(3:Dim+2,3:Dim+2)-D;

min(eig(O(3:Dim,3:Dim)))



  