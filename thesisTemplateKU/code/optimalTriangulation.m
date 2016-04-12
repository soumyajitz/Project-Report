function  X  = optimalTriangulation( x1,x2,P1,P2,F,e1,e2 )
% [ X,r ] = optimalTriangulation( x1,x2,P1,P2,F,e1,e2 )
%   Optimal solution for MLE, using method proposed by Hartley(P315).
%   Finding optimal solution by finding real solutions of a polynomial equation of degree  6.
%   

if any(size(x1)~=size(x2))
    error('size(x1)~=size(x2)');
end

if size(x1,1)==2  % convert x1,x2 to form of (x,y,1)
    x1=[x1;1];
    x2=[x2;1];
else
    x1(1)=x1(1)./x1(3);x1(2)=x1(2)./x1(3);x1(3)=1;
    x2(1)=x2(1)./x2(3);x2(2)=x2(2)./x2(3);x2(3)=1;
end

T1=[1 0 -x1(1);0 1 -x1(2);0 0 1];T1_inv=[1 0 x1(1);0 1 x1(2);0 0 1];
T2=[1 0 -x2(1);0 1 -x2(2);0 0 1];T2_inv=[1 0 x2(1);0 1 x2(2);0 0 1];
F=T2_inv.'*F*T1_inv;
%x1=[0;0;1];%T1*x1;
%x2=[0;0;1];%T2*x2;
e1=T1*e1;e1=e1/norm(e1(1:2));
e2=T2*e2;e2=e2/norm(e2(1:2));
R1=[e1(1) e1(2) 0;-e1(2) e1(1) 0;0 0 1];
R2=[e2(1) e2(2) 0;-e2(2) e2(1) 0;0 0 1];
%e1=[1;0;e1(3)];%R1*e1;
%e2=[1;0;e2(3)];%R2*e2;
F=R2*F*R1.';
f1=e1(3);f2=e2(3);a=F(2,2);b=F(2,3);c=F(3,2);d=F(3,3);
%F_=[f1*f2*d -f2*c -f2*d;-f1*b a b;-f1*d c d];
%F./F_


f1_2=f1^2;f1_4=f1_2^2;f2_2=f2^2;
k1=a^2+f2_2*c^2;k2=2*(a*b+f2_2*c*d);k3=b^2+f2_2*d^2;
k=a*d-b*c;k4=a*c;k5=b*c+a*d;k6=b*d;
c6=-f1_4*k*k4;
c5=k1^2 - f1_4*k*k5;
c4=2*k1*k2 - 2*f1_2*k*k4 - f1_4*k*k6;
c3=k2^2 + 2*k1*k3 - 2*f1_2*k*k5;
c2=2*k2*k3 - k*k4 - 2*f1_2*k*k6;
c1=k3^2 - k*k5;
c0=- k*k6;
% Threshold to set some very samll cis to 0.
p=[c6 c5 c4 c3 c2 c1 c0];
Mci=max(abs(p));
for i=1:6
    if abs(p(i))<Mci*1e-25  % threshold ratio is chosen to be 1e-20
        p(i)=0;
    end
end
if sum(isnan(p))~=0 || sum(isinf(p))~=0
    X=zeros(4,1);return
end
rs=real(roots(p));
nrs=size(rs,1);% number of roots.
s=zeros(nrs+1,1);
for i=1:nrs
    t=rs(i);
    t2=t^2;
    k1_=(a*t+b)^2;
    k2_=(c*t+d)^2;
    s(i)=t2/(1+f1_2*t2)+k2_/(k1_+f2_2*k2_);
end
s(nrs+1)=1/f1_2+c^2/(a^2+f2_2*c^2);
[m I]=min(s);
if I==nrs+1
    l1=[f1;0;-1];
    l2=[-c*f2;a;c];
else
    t=rs(I);
    l1=[t*f1;1;-t];
    l2=[-f2*(c*t+d);a*t+b;c*t+d];
end
x1_=[-l1(1)*l1(3);-l1(2)*l1(3);l1(1)^2+l1(2)^2];
x2_=[-l2(1)*l2(3);-l2(2)*l2(3);l2(1)^2+l2(2)^2];
%r=sqrt((x1_(1)^2+x1_(2)^2)/x1_(3)^2+(x2_(1)^2+x2_(2)^2)/x2_(3)^2);% residual
%r=sqrt(m);
x1_=T1_inv*R1.'*x1_;
x2_=T2_inv*R2.'*x2_;

Ax=zeros(4,4);
Ax(1,:)=x1_(1)*P1(3,:)-x1_(3)*P1(1,:);
Ax(2,:)=x1_(2)*P1(3,:)-x1_(3)*P1(2,:);
Ax(3,:)=x2_(1)*P2(3,:)-x2_(3)*P2(1,:);
Ax(4,:)=x2_(2)*P2(3,:)-x2_(3)*P2(2,:);
[u s v]=svd(Ax);
X=v(:,end);

end

