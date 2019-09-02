clear all 
clc

kx=1;
ky=1;

load coor.txt;
load noc.txt;
coor=coor(:,2:end);
noc=noc(:,2:end);
f=@(x) -x.^2+10*x;

bc=[1:11,12,22,23,33,34,44,45,55,56:66];
Tbc=[zeros(1,19),f(coor(56:66,1)')];


nn=size(coor,1);
ne=size(noc,1);
Ka=zeros(nn);

for k=1:ne
      xel(1:3)=coor(noc(k,:),1);
	  yel(1:3)=coor(noc(k,:),2) ;
	  detj= (xel(1)-xel(3)) * (yel(2)-yel(3)) - (xel(2)-xel(3)) * (yel(1)-yel(3));
      
      dN_dx(1,1)=(yel(2)-yel(3))/detj;
	  dN_dx(2,1)=(yel(3)-yel(1))/detj;
	  dN_dx(3,1)=(yel(1)-yel(2))/detj;

      dN_dy(1,1)=(xel(3)-xel(2))/detj;
      dN_dy(2,1)=(xel(1)-xel(3))/detj;
      dN_dy(3,1)=(xel(2)-xel(1))/detj;
      
      ke=(kx*dN_dx*dN_dx'+ky*dN_dy*dN_dy')*detj*.5; 
      q=noc(k,:);
      Ka(q,q)=Ka(q,q)+ke; 
end

Ka(bc,:)=0;
Fa=-Ka(:,bc)*Tbc';
Ka(:,bc)=0;
Fa(bc)=Tbc;


for i=1:size(bc,2)
  Ka(bc(i),bc(i))=1.;
end
 
T=Ka\Fa;
x=coor(:,1);
y=coor(:,2);
disp('        x          y         T')
disp([x,y,T])




