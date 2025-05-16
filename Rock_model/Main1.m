clc; clear; close all;

%%
% define the uncracked rock's properties
Vp1   = 4.0;
Vs1   = 2.2;
Dens1 = 2.4;

% define parameters for inclusions
alpha1 = 0.01;
alpha2 = 0.02;
alpha3 = 0.03;
alpha4 = 0.04;
alpha5 = 0.05;
alpha6 = 0.06;
alpha7 = 0.07;
alpha8 = 0.08;
alpha9 = 0.09;
alpha10 = 0.1;

Dens2 = 1.0;
Vp2   = 1.5;
Vs2   = 0;

% [Vp_KT,Vs_KT]=crack_KTB(Dens1,Vp1,Vs1,PHI,alpha,Dens2,Vp2,Vs2,4);

%%
P = 0:0.002:0.5;
Vs_1 = zeros(1,length(P));
Vs_2 = zeros(1,length(P));
Vs_3 = zeros(1,length(P));
Vs_4 = zeros(1,length(P));
Vs_5 = zeros(1,length(P));
Vs_6 = zeros(1,length(P));
Vs_7 = zeros(1,length(P));
Vs_8 = zeros(1,length(P));
Vs_9 = zeros(1,length(P));
Vs_10 = zeros(1,length(P));


for i=1:length(P)
    PHI = P(i);
    [~,Vs_KT1]=crack_KTB_1(Dens1,Vp1,Vs1,PHI,alpha1,Dens2,Vp2,Vs2,4);
    [~,Vs_KT2]=crack_KTB_1(Dens1,Vp1,Vs1,PHI,alpha2,Dens2,Vp2,Vs2,4);
    [~,Vs_KT3]=crack_KTB_1(Dens1,Vp1,Vs1,PHI,alpha3,Dens2,Vp2,Vs2,4);
    [~,Vs_KT4]=crack_KTB_1(Dens1,Vp1,Vs1,PHI,alpha4,Dens2,Vp2,Vs2,4);
    [~,Vs_KT5]=crack_KTB_1(Dens1,Vp1,Vs1,PHI,alpha5,Dens2,Vp2,Vs2,4);
    [~,Vs_KT6]=crack_KTB_1(Dens1,Vp1,Vs1,PHI,alpha6,Dens2,Vp2,Vs2,4);
    [~,Vs_KT7]=crack_KTB_1(Dens1,Vp1,Vs1,PHI,alpha7,Dens2,Vp2,Vs2,4);
    [~,Vs_KT8]=crack_KTB_1(Dens1,Vp1,Vs1,PHI,alpha8,Dens2,Vp2,Vs2,4);
    [~,Vs_KT9]=crack_KTB_1(Dens1,Vp1,Vs1,PHI,alpha9,Dens2,Vp2,Vs2,4);
    [~,Vs_KT10]=crack_KTB_1(Dens1,Vp1,Vs1,PHI,alpha10,Dens2,Vp2,Vs2,4);


    Vs_1(i) = Vs_KT1;
    Vs_2(i) = Vs_KT2;
    Vs_3(i) = Vs_KT3;
    Vs_4(i) = Vs_KT4;
    Vs_5(i) = Vs_KT5;
    Vs_6(i) = Vs_KT6;
    Vs_7(i) = Vs_KT7;
    Vs_8(i) = Vs_KT8;
    Vs_9(i) = Vs_KT9;
    Vs_10(i) = Vs_KT10;
end

%%
x_shade1 = [0.181 0.181 0.229 0.229];
y_shade1 = [0.0 4.5 4.5 0.0];

x_shade2 = [0.00 0.00 0.5 0.5];
y_shade2 = [1.4 2.0 2.0 1.4];

%%
figure
hold on;
% fill(x_shade1, y_shade1, [0.3, 0.5, 0.7], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
fill(x_shade2, y_shade2, [0.8, 0.4, 0.4], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

plot(P,Vs_1,'b:o', 'LineWidth', 2, 'MarkerSize', 5,'color',[0.3, 0.1, 0.7], 'MarkerFaceColor', 'k');
plot(P,Vs_2,'b:o', 'LineWidth', 2, 'MarkerSize', 5,'color',[0.3, 0.2, 0.7], 'MarkerFaceColor', 'k');
plot(P,Vs_3,'b:o', 'LineWidth', 2, 'MarkerSize', 5,'color',[0.3, 0.3, 0.7], 'MarkerFaceColor', 'k');
plot(P,Vs_4,'b:o', 'LineWidth', 2, 'MarkerSize', 5,'color',[0.3, 0.4, 0.7], 'MarkerFaceColor', 'k');
plot(P,Vs_5,'b:o', 'LineWidth', 2, 'MarkerSize', 5,'color',[0.3, 0.5, 0.7], 'MarkerFaceColor', 'k');
plot(P,Vs_6,'b:o', 'LineWidth', 2, 'MarkerSize', 5,'color',[0.3, 0.6, 0.7], 'MarkerFaceColor', 'k');
plot(P,Vs_7,'b:o', 'LineWidth', 2, 'MarkerSize', 5,'color',[0.3, 0.7, 0.7], 'MarkerFaceColor', 'k');
plot(P,Vs_8,'b:o', 'LineWidth', 2, 'MarkerSize', 5,'color',[0.3, 0.8, 0.7], 'MarkerFaceColor', 'k');
plot(P,Vs_9,'b:o', 'LineWidth', 2, 'MarkerSize', 5,'color',[0.3, 0.9, 0.7], 'MarkerFaceColor', 'k');
plot(P,Vs_10,'b:o', 'LineWidth', 2, 'MarkerSize', 5,'color',[0.3, 1.0, 0.7], 'MarkerFaceColor', 'k');

legend('','Alpha=0.01','Alpha=0.02','Alpha=0.03','Alpha=0.04', ...
    'Alpha=0.05','Alpha=0.06','Alpha=0.07','Alpha=0.08', ...
    'Alpha=0.09','Alpha=0.10');
xticks(0:0.05:0.5);
yticks(0:0.5:4.0);
yticklabels((0:0.5:4.0)*1000);
xlim([0 0.5]);
ylim([0.5 2.5]);
ylabel('S-wave velocity (m/s)'); xlabel('Fracture density');
grid on; box on;
set(gca,'linewidth',2,'fontsize',18.5,'fontname','Arial','TickLength',[0.007 0.01]);
set(gcf,'unit','centimeters','position',[10,10,24,15]);

set(gcf,'PaperSize',[15 12]);
print(gcf,'Fig4a','-dpdf','-r0');
drawnow;




