%ME303 LAB
%Numerical integration using Simpson's rule
%This method gives usually better results than trapezoidal rule

clc
clear all
close all

f=@(x) 2+sin(2*sqrt(x));
a=1;
b=6;
M=10;
h=(b-a)/(2*M);
x=a:h:b;
s=zeros(1,M);

for i=1:M
    x1=x(2*i-1);
    x2=x(2*i);
    x3=x(2*i+1);
    f1=f(x1);
    f2=f(x2);
    f3=f(x3);
    s(i)=(h/3)*(f1+4*f2+f3);
end

S=sum(s)

dx=0.001;
xp=a:dx:b;
figure
hold on
grid on
plot(xp,f(xp),'b','linewidth',3)
plot(x,f(x),'*r','linewidth',2)
plot([a b],[0 0],'k','linewidth',2)

for i=1:M    % for each sub-area of integration we define 3 points
    x1=x(2*i-1);
    x2=x(2*i);
    x3=x(2*i+1);
    f1=f(x1);
    f2=f(x2);
    f3=f(x3);
    xpi=x1:dx:x3;
    if mod(i,2)==0                % This bit colors the areas of numerical integration in alternating colors
        area(xpi,f(xpi),'facecolor',[0.1 0.5 0.9]); % "area" function paints the area under a curve
    else                                            % You can set the color of the areas using the argument "facecolor" and specifying the rgb values
        area(xpi,f(xpi),'facecolor',[0.1 0.1 0.6]);
    end
    plot([x1 x1],[0 f1],'--k');
    plot([x3 x3],[0 f3],'--k');
    text(x2,(f1+f3)/4,num2str(round(s(i),3)),'FontSize',12,'color',[1, 0.5, 0],'VerticalAlignment', 'top','HorizontalAlignment','center');
    % Using "text" function to show the value of the numerical integration for each area
end
text((a+b)/2,max(f(xp)),['Integral=' num2str(S)],'FontSize',20,'color',[1, 0.5, 0],'VerticalAlignment', 'top','HorizontalAlignment','center');



