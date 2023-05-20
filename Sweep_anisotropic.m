
%FAST SWEEPING ALGORITHMS FOR A CLASS OF HAMILTON¡VJACOBI EQUATIONS
%SIAM J. NUMER. ANAL. c 2003 Society for Industrial and Applied Mathematics
%Vol. 41, No. 2, pp. 673¡V694
% [D,L1error,X,Y,Errorscat]=Sweep_anisotropic(101,2,2);
%D: time map 
%Bx, By: boundary 
% metric tensor = [a,c;c,b]
function [D,L1error,X,Y,Errorscat2]=Sweep_anisotropic(gridnum,Bx,By)
%clear 
n=gridnum;
%  l=2;w=2;
%  Bx=2;By=2;
%  n=51;
x = linspace(-Bx,Bx,n);
y = linspace(By,-By,n);
[X,Y] = meshgrid(x,y);

%==================================================
    
neigh = [[-1;0] [1;0] [0;1] [0;-1] [-1 ;-1] [1; -1] [-1; 1] [1 ;1]];
boundary = @(x)mod(x-1,n)+1;
ind2sub1 = @(k)[rem(k-1, n)+1; (k - rem(k-1, n) - 1)/n + 1];
sub2ind1 = @(u)(u(2)-1)*n+u(1);
Neigh = @(k,i) sub2ind1( boundary(ind2sub1(k)+neigh(:,i)) );

I=find(X==0 & Y== 0);
%D:distance map
D = zeros(n)+Inf;
D(I) = 0;

%reference point:four corner
% refere = [(n+1)/2,(n+1)/2];
refere = [1,1;size(X,1),1;1 size(X,2);size(X,1) size(X,2)];
for Seild = 1:size(refere,1)
l2metric = sqrt((X-X(refere(Seild,1),refere(Seild,2))).^2+(Y-Y(refere(Seild,1),refere(Seild,2))).^2);
l2metric = reshape(l2metric,1,size(l2metric,1)*size(l2metric,2));
[l2Order,l2index]= sort(l2metric);
[iniPosrr,iniPoscc]=find(l2index==I);
l2index(iniPoscc)=[];
SweepDir{Seild}=l2index;%Sweeping direction
end

L1error = [Inf Inf];
 
 while (isnan(L1error(end)) | isinf(L1error(end)) | abs(L1error(end)-L1error(length(L1error)-1)) >10^-10)

Dn=D;
tic
for Sweep =1:4   
%for ii=2:length(l2index)%29:-1:1%29:nodex
tt = SweepDir{Sweep};
for ii=1:size(SweepDir{Sweep},2)  
i = tt(ii);
%ind2sub1(i);
 posxm = ind2sub1(i);

Tabc = [];
 xm=[8 2;3 8; 3 7;7 1; 5 1; 4 5;4 6; 6 2];
if posxm(1)>1 && posxm(1)<n && posxm(2)>1 &&posxm(2)<n
    SS = 1:8;
elseif posxm(1)==1 && posxm(2)==1
    SS = [1 2];
elseif posxm(1)==1 && posxm(2)==n
    SS = [7 8];
elseif posxm(1)==n && posxm(2)==1
    SS = [3 4];
elseif posxm(1)==n && posxm(2)==n
    SS = [5 6];

elseif posxm(1)==1
    SS = [1 2 7 8];
elseif posxm(1)==n
    SS = [3 4 5 6];
elseif posxm(2)==1
    SS = [1 2 3 4];
elseif posxm(2)==n
    SS = [5 6 7 8];
end

for j= SS
%i = sub2ind1([51,50]);    
% for j = 1
uABC = [];
xa = X(Neigh(i,xm(j,1)));
xb = X(Neigh(i, xm(j,2)));
ya = Y(Neigh(i,xm(j,1)));
yb = Y(Neigh(i, xm(j,2)));
xc = X(i); yc = Y(i);
%xc=4.5;yc=2

a =1;
c= 0.9;
b= 1;

M = [a -c;-c b];%anisotropic metric

Tc = D(i);Ta = D(Neigh(i,xm(j,1))); Tb = D(Neigh(i, xm(j,2)));
la = sqrt( (xc-xb)^2+(yc-yb)^2 );
lb = sqrt( (xc-xa)^2+(yc-ya)^2 );
lc = sqrt( (xb-xa)^2+(yb-ya)^2 );
n11=(xc-xa)/lb;n12 = (yc-ya)/lb;
n21 = (xc-xb)/la;n22= (yc-yb)/la;
r = acos((la^2+lb^2-lc^2)/2/la/lb);

Pinv = 1/sin(r)^2*[n11-n21*cos(r),n21-n11*cos(r);n12-n22*cos(r),n22-n12*cos(r)];
g1 = Pinv(1)/lb+Pinv(3)/la;
g2 = -(Pinv(1)/lb*Ta+ Pinv(3)/la*Tb);
g3 = (Pinv(2)/lb+Pinv(4)/la);
g4 = -(Pinv(2)/lb*Ta+ Pinv(4)/la*Tb);
w1 = a*g1^2+b*g3^2-2*c*g1*g3;
w2 = 2*a*g1*g2+2*b*g3*g4-2*c*(g1*g4+g2*g3);
w3 = a*g2^2+b*g4^2-2*c*g2*g4;
TC = [ (-w2+sqrt(w2^2-4*w1*(w3-1)) )/2/w1,(-w2-sqrt(w2^2-4*w1*(w3-1)) )/2/w1 ];
%check for two root
for kk = 1:2
    % kk=2
    if  w2^2-4*w1*(w3-1)>0 && TC(kk)>0
        %isreal(TC(kk)) ==1 && TC(kk)>0 && isnan(TC(kk))~=1 && isinf(TC(kk))~=1 
        pq = [g1*TC(kk)+g2; g3*TC(kk)+g4];%Pinv*[(Tc(kk)-TA)/lb; (Tc(kk)-TB)/la ];
        d = M*pq;
      %if d(2)/d(1)< max([(yc-yb)/(xc-xb),(yc-ya)/(xc-xa)]) && d(2)/d(1)>min((yc-yb)/(xc-xb),(yc-ya)/(xc-xa)) && d(2)< max([yc-yb,yc-ya]) && d(2)>min([yc-yb,yc-ya])
      if j==1 && d(2)>0 && d(1)<0
         Dir  = 1;
      elseif  j==2 && d(2)>0 && d(1)<0
         Dir=1;
      elseif j==3 && d(2)<0 && d(1)<0
          Dir=1;
      elseif j==4 && d(2)<0 && d(1)<0
          Dir=1;
      elseif j==5 && d(2)<0 && d(1)>0
          Dir=1;
      elseif j==6 && d(2)<0 && d(1)>0
          Dir=1;
      elseif j==7 && d(2)>0 && d(1)>0
          Dir=1;
      elseif j==8 && d(2)>0 && d(1)>0
          Dir=1;
      else
          Dir=0;
      end
        
   if Dir==1       
      if xb-xa==0 
          Xnew = xb;
          Ynew =  yc+d(2)/d(1)*(Xnew-xc);
       elseif yb-ya==0
          Xnew = (ya-yc-(yb-ya)/(xb-xa)*xa+d(2)/d(1)*xc)/(d(2)/d(1)-(yb-ya)/(xb-xa));
          Ynew = yb;
       else
       Xnew = (ya-yc-(yb-ya)/(xb-xa)*xa+d(2)/d(1)*xc)/(d(2)/d(1)-(yb-ya)/(xb-xa));
       Ynew = yc+d(2)/d(1)*(Xnew-xc);
       end
     if  isnan(Xnew)~=1 && Xnew<=max(xa,xb) && Xnew>=min(xa,xb)  && Ynew<=max(ya,yb) && Ynew>=min(ya,yb)
     uABC = [uABC,TC(kk)];
     end
   end
    end
end
uABC = min(uABC);
%one value for single triangle 
if isempty(uABC )==1         
     uABC = min([D(i),Ta+ sqrt([xc-xa;yc-ya].'*inv(M)*[xc-xa;yc-ya]), Tb+ sqrt( [xc-xb;yc-yb].'*inv(M)*[xc-xb;yc-yb])] );
else
     uABC = min([D(i),uABC]);
end
%8 value from 8 triangle
Tabc =[Tabc, uABC];
end
D(i) = min(Tabc);
end
end
toc
Order=[];
for hh=1:size(D,2)-2
    Order= [Order,(2:size(D,1)-1)+hh*n];
end
  for ii= 1:length(Order)
   for ieSS = 1:8
       ie = Order(ii);
   Errors_eachtriangle =( abs(D(ie)-Dn(ie))+abs((D(Neigh(ie,xm(ieSS,1))))-Dn(Neigh(ie,xm(ieSS,1))))+abs((D(Neigh(ie,xm(ieSS,2)))-Dn(Neigh(ie,xm(ieSS,2))))) )/3;
   Errorscat(ii,ieSS) =Errors_eachtriangle/8;
   Errorscat2(ii) = sum(Errorscat(ii,ieSS));
   end
  end

 L1error = [L1error,sum(sum(Errorscat))/size(Order,2)];

 %if length(L1error)>200
 %    break;
 %end
 end

figure;
[C,h]=contour(X,Y,D)
set (h, 'LineWidth', 3);
set(gca,'xtick',linspace(-2,2,5))
set(gca,'ytick',linspace(-2,2,5))
set(gca,'fontsize',18)
xlabel('x');
ylabel('y')
