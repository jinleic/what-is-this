%phi=(sqrt(5)-1)/2;
%bphi=1-phi;
dt=0.0004;
S=[dt:dt:1-dt];
n=length(S);
T= S;
 a=zeros(size(T));
% for k=1:length(S)
%     t=T(k);
%     if t>1-1/sqrt(2) && t<= 0.5
%         a(k)= sqrt((1-2*(1-t)^2)/2/t/(1-t));
%     end
%     if t>0.5
%         a(k) = 1;
%     end
% end
%a=(1-S)*0;
a=S.*(1-S);%*3.2;
H=zeros(length(S),length(T));
for i=1:length(S)
    for j=1:length(T)
        s=S(i);
        t=T(j);
        bs=1-s;
        bt=1-t;
        %z=bs*bt;
        %z=bs*bt+a(i)*a(j)*(min(bs,bt)-bs*bt);
        z=bs*bt+a(i)*a(j);
        %z=bs*bt/(bs+bt-bs*bt);
        H(i,j)=-z*log2(z)-(1-z)*log2(1-z);
    end
end
%mesh(H)
A=zeros(length(S),3);
for i=1:n
    A(i,1)=1;
    A(i,2)=i/n;
    A(i,3)=a(i);
end
P=eye(n,n)-A*(A'*A)^(-1)*A';
H1=P(4:n,:)*H*P(:,4:n);
%H1=P(3:n,:)*H*P(:,3:n);
H1=(H1+H1')/2;
max(eig(H1))
%max(real(eig(H1)))

%dt=0.0004; a=S.*(1-S);
%1.6311e-14 
%if use H1+H1'/2 2.3685e-14
%if not removing f() eigenspace
%    8.2737
%compare eps= 2.2204e-16, random matrix (2500*2500) 100

%independent coupling:
%1.6099e-14 
%if use H1+H1'/2  2.4206e-14

