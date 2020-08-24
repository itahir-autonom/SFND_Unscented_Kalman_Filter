data_lidar = load('NIS_lidar.txt');
data_radar = load('NIS_radar.txt');

x_radar=7.815.*ones(length(data_lidar),1);
x_lidar=5.991.*ones(length(data_lidar),1);

figure(1);
plot(data_lidar);
title("NIS Lidar");
hold on;
plot(x_lidar,'color','r','linestyle','--','linewidth',5);


figure(2)
plot(data_radar);
title("NIS Radar");
ylim([1 12])
set(gca,'YTick',0:4:12)
hold on;
plot(x_radar,'color','r','linestyle','--','linewidth',5);