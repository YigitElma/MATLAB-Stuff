clc
clear
close all

cetin_1_max = [40
50
60
72
];
cetin_1_per = [5
11
19
27
];

cetin_2_max = [40
58
70
77
];
cetin_2_per = [6
21
31
38
];

deney_us_max = [30
40
50
60
70
75.5
];
deney_us_per = [0
8.5
16
24
33
38.5
];

ansys_max = [37.5
41
46
51
55.5
60
63.5
67
74
78
];
ansys_per = [6
8.5
12
16
19
23
26
29
35
38
];

ansys = @(x) 0.798*x - 24.52 ;
cet_2 = @(x) 0.8567*x - 28.47 ;
us = @(x) 0.8345*x - 25.27 ;
cet_1 = @(x) 0.698*x - 23.24 ;

x = 30:0.1:90 ;

hold on
% plot(ansys_max,ansys_per,'r--')
% plot(cetin_2_max,cetin_2_per,'b--')
% plot(deney_us_max,deney_us_per,'k--')

% plot(x,ansys(x),'r')
plot(x,us(x),'k')
plot(x,cet_2(x),'b')

% plot(ansys_max,ansys_per,'r*','MarkerSize',10)
plot(cetin_2_max,cetin_2_per,'b^','MarkerSize',10)
plot(deney_us_max,deney_us_per,'k.','MarkerSize',25)

plot(x,cet_1(x),'g')
plot(cetin_1_max,cetin_1_per,'g.','MarkerSize',25)

grid on
xlabel('Maximum Bending Angle (degree)','FontSize',13,'FontWeight','bold')
ylabel('Permanent Bending Angle (degree)','FontSize',13,'FontWeight','bold')
legend('Experiment 1 [\theta_{per}=0.8567\theta_{max} - 28.47]','Experiment 2 [\theta_{per}=0.8345\theta_{max} - 25.27]',...
    'Exp 1 Data','Exp 2 Data','Experiment 3 [\theta_{per}=0.698\theta_{max} - 23.24]','Exp 3 Data','FontSize',13,'FontWeight','bold')





