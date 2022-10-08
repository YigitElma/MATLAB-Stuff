clc
clear all
close all

omega=1;
ksi=[0 0.1 0.7 1 5];

for i=1:5
    H(i)=tf([2], [1 2*omega*ksi(i) 1]);  
    fprintf('poles of H%d are\n', i)  %i print out the poles by finding roots of denominator
    roots([1 2*omega*ksi(i) 1])
    %poles of H1 are pure complex so step response is undamped
    %poles of H2 are complex with negative reel part so underdamped
    %poles of H3 are complex with negative reel part so underdamped
    %poles of H4 are negative reel and same so critically damped
    %poles of H5 are negative reel and different so overdamped
end
step(H(1),H(2),H(3),H(4),H(5),20)
figure()
ksi2=0.8;
TS=4/ksi2/omega  %settling time
TR=(1-0.4167*ksi2+2.917*ksi2^2)/omega  %rise time
SSE=1/3  %Steady state error=1/(1+G(s)) as s->0 ==> sse=(s^2+1.6*s+1)/(s^2+1.6*s+3)
DCgain=2/omega^2  
OS=exp(-pi*ksi2/sqrt(1-ksi2^2))*100 %percent overshoot
tf2=tf([2],[1 2*omega*ksi2 1]);
step(tf2,20)
stepinfo(tf2)
