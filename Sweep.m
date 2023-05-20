
%The fast sweeping method in two dimensions. 
%S. Luo, and J. Qian, "Fast Sweeping Methods for Factored Anisotropic 
%Eikonal Equations: Multiplicative and Additive Factors," J Sci Comput 52, 360-382 (2012)
%(a) shows the arrival times after sweeping up-right. 
%(b) is an image of the arrival times after a further sweep up-left. 
%(c) displays the arrival times, T, after a further sweep down right. 
%(d) shows the final state of the arrival times, which is obtained after
%a further sweep down-left. A fifth sweep is required which does not
%update any values for the method to reach its termination condition.
clear 
% 网格数
n = 50;
W = ones(n);
x0=[n/2;n/2];
% interested domain
x = linspace(-1,1,n);
y = linspace(-1,1,n);
[X,Y] = meshgrid(x,y);
neigh = [[-1;0] [1;0] [0;1] [0;-1] [-1 ;-1] [1; -1] [-1; 1] [1 ;1]];
boundary = @(x)mod(x-1,n)+1;
ind2sub1 = @(k)[rem(k-1, n)+1; (k - rem(k-1, n) - 1)/n + 1];
sub2ind1 = @(u)(u(2)-1)*n+u(1);
Neigh = @(k,i) sub2ind1( boundary(ind2sub1(k)+neigh(:,i)) );


I = sub2ind1(x0);%19*40+20
%u:distance map
u = zeros(n)+Inf;
u(I) = 0;

[nodex,nodey]=size(u);
h=1;

[nodex,nodey]=size(u);
Order_j = [1:nodex;1:nodex;nodex:-1:1;nodex:-1:1];
Order_i = [1:nodey;nodey:-1:1;1:nodey;nodey:-1:1];
% sweeping 算法
for O = 2:4 % 不同扫描顺序 %%%调整following number 1~4改变扫描图
for i=Order_i(O,1:50) % y
for j=Order_j(O,:)    % column
% The Gauss Seidel loop
% Only carry it out if the node hasn't been set initially
if ~(isnan(u(i,j)))
if ((i>1) && (i<nodey))
uymin = min(u(i-1,j),u(i+1,j));
elseif (i==1)
uymin = u(2,j);
else
uymin = u(nodey-1,j);
end
if ((j>1)&&(j<nodex))
uxmin = min(u(i,j-1),u(i,j+1));
elseif (j==1)
uxmin = u(i,2);
else
uxmin = u(i,nodex-1);
end
if (abs(uxmin-uymin)>=h)
ubar = min(uxmin,uymin)+h;
else
ubar = (uxmin+uymin+sqrt(2*h^2-(uxmin-uymin)^2))/2;
end
u(i,j)=min(u(i,j),ubar);
end %if
end
end
end
figure;
[C,h]=contour(x,y,u)
set (h, 'LineWidth', 3);
set(gca,'xtick',-1:0.2:1)
set(gca,'ytick',[-1:0.2:1])
set(gca,'fontsize',10)
% xlabel('x');
% ylabel('y')
