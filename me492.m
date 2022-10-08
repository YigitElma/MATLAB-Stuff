clc
clear 
close all


moment = [380
390
400
410
420
422.5
425
427.5
429
430
432.5
435
437.5
440
];

max_def = [1.1009
1.1771
1.3074
1.4686
1.6456
1.6962
1.7602
1.8438
1.9058
1.9501
2.0865
2.2605
2.4589
2.6958
];

plastic = [0.18876
0.24198
0.35037
0.49041
0.64713
0.69303
0.7526
0.8326
0.8925
0.9356
1.0697
1.2421
1.4397
1.6765
];
for i=1:length(plastic)
    spring_back_ratio(i) = (max_def(i)-plastic(i))/max_def(i);
end

plot(moment,max_def)
hold on
plot(moment,plastic)
hold on
% plot(moment,spring_back_ratio)
grid on
xlabel("MOMENT (Nmm)")
ylabel("DEFLECTION (mm)")
legend("Maximum deflection","Plastic Deflection")
title("SPRINGBACK EFFECT")



