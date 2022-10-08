clc
clear all
close all

x=[ 0.9958    0.9978    0.9040    0.7209    0.4571    0.1240   -0.2638 ...
   -0.6891   -1.1324   -1.5728   -1.9888   -2.3590   -2.6634   -2.8838 ...
   -3.0053   -3.0162   -2.9094   -2.6823   -2.3373   -1.8820   -1.3286 ...
   -0.6943   -0.0000    0.7299    1.4686    2.1878    2.8588    3.4536 ...
    3.9461    4.3133    4.5359    4.5998    4.4962    4.2227    3.7830 ...
    3.1874    2.4527    1.6013    0.6608   -0.3368   -1.3564   -2.3608 ...
   -3.3124   -4.1743   -4.9116   -5.4933   -5.8931   -6.0907   -6.0726 ...
   -5.8331   -5.3746   -4.7073   -3.8498   -2.8281   -1.6748   -0.4287    0.8674    2.1675];

y=[ 0.3160    0.4685    0.6221    0.7660    0.8901    0.9855    1.0442    1.0605 ...
    1.0299    0.9506    0.8226    0.6484    0.4327    0.1822   -0.0945   -0.3872 ...
   -0.6845   -0.9744   -1.2445   -1.4828   -1.6779   -1.8199   -1.9003   -1.9132 ...
   -1.8547   -1.7237   -1.5221   -1.2546   -0.9284   -0.5537   -0.1427    0.2905 ...
    0.7305    1.1607    1.5648    1.9265    2.2308    2.4642    2.6155    2.6764 ...
    2.6413    2.5085    2.2796    1.9600    1.5585    1.0875    0.5621    0.0000 ...
   -0.5792   -1.1548   -1.7054   -2.2102   -2.6494   -3.0050   -3.2616   -3.4070   -3.4331   -3.3354];

m=0.025;  %mass of the particle
t0=0.15;  %actually we don't need these in code but since it's given I wrote anyway
tf=3;   %final data time
delta_t=0.05;  %sampling frequency
N=length(x);   %I could find same result by (tf-t0)/delta_t+1

fx=zeros(1,N);  %for each point we want force values so I create 1*N matrix
fy=zeros(1,N);  %since motion is in 2D we have fx and fy
% fx(1)=0;  
% fx(end)=0;
% fy(1)=0;
% fy(end)=0;  % I commented these lines because these values are already
% zeros, don't need to assign them 0 again

for i=2:N-1   %first and last points's forces cannot be calcultaed by these method because it needs 2 adjacent point too
    fx(i)=m*((x(i+1)-2*x(i)+x(i-1))/delta_t^2);  %I used the O(h^2) method in the book
    fy(i)=m*((y(i+1)-2*y(i)+y(i-1))/delta_t^2);  %we want derivative with respect to time so we need x and y
                                            %values in time domain, but luckly sampling is periodic
                                            %and x(i) actually means x(t0+(i-1)*delta_t), so we are good in terms of domains                                                                                                                                 
end

plot(x,y)
hold on
quiver(x,y,fx,fy,0) %this built-in function draws [fx fy] vector at corresponding [x y] point



