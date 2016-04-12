function [F,inliers] =RobustFundaMatrix(x1,x2)
numpts=size(x1,2);
[Fi,e1i,e2i]=fundmatrix(x1,x2);
Pai=eye(3,4);Pbi = vgg_P_from_F(Fi);
Xe=zeros(4,numpts);
err=zeros(1,numpts);
et=err;
for i=1:numpts
    Xe(:,i)=optimalTriangulation(x1(:,i),x2(:,i),Pai,Pbi,Fi,e1i,e2i);
    err(i)=sqrt(norm(pflat(Pai*Xe(:,i))-x1(:,i))^2+norm(pflat(Pbi*Xe(:,i))-x2(:,i))^2);
    et(i)=sum(pflat(Pai*Xe(:,i))-x1(:,i))+sum(pflat(Pbi*Xe(:,i))-x2(:,i));
end
m=median(abs(err));
sigma=1.4826*(1+5/(numpts-8))*m;
index=find(abs(err)<3*sigma);
[F,e1,e2]=fundmatrix(x1(:,index),x2(:,index));
Pa=eye(3,4);Pb = vgg_P_from_F(F);
for i=1:numpts
    Xe(:,i)=optimalTriangulation(x1(:,i),x2(:,i),Pa,Pb,F,e1,e2);        
    err(i)=sqrt(norm(pflat(Pa*Xe(:,i))-x1(:,i))^2+norm(pflat(Pb*Xe(:,i))-x2(:,i))^2);
    et(i)=sum(pflat(Pai*Xe(:,i))-x1(:,i))+sum(pflat(Pbi*Xe(:,i))-x2(:,i));
end
m=median(abs(err));
sigma=1.4826*(1+5/(numpts-8))*m;
index=find(abs(err)<3*sigma);
F=fundmatrix(x1(:,index),x2(:,index));
inliers=index;