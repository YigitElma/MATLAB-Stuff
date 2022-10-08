clc
clear all
close all

%for zate values given in the question
zeta = 0:0.2:1;
%common time matrix for step responses
time = 0:0.1:50;
N = length(time);
%for 3D plotting step responses
zeta1 = linspace(0,1,N);

for i=1:N   %this loop generates responses for 3D plot
    y(i) = tf(1, [ 1 2*zeta1(i) 1 ]);
    z = step(y(i),time); 
    k(:,i) = z;
end
for i=1:6   %this loop generates plot for given zeta values
    h = tf(1, [ 1 2*zeta(i) 1 ]);
    hold on
    grid on
    step(h,time)
end
figure
s = surf(zeta1,time,k);
s.EdgeColor = 'flat';
xlabel('zeta')
ylabel('time')
zlabel('response')