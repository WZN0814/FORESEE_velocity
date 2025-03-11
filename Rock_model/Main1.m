clc; clear; close all;

%%
% define the uncracked rock's properties
Vp1   = 4.0;
Vs1   = 2.2;
Dens1 = 2.4;

% define parameters for inclusions
alpha = 0.01;
Dens2 = 1.0;
Vp2   = 1.5;
Vs2   = 0;

%%
P = 0:0.002:0.05;
Vs = zeros(1,length(P));
for i=1:length(P)
    PHI = P(i);
    [Vp_KT,Vs_KT]=crack_KTB_1(Dens1,Vp1,Vs1,PHI,alpha,Dens2,Vp2,Vs2,4);
    Vs(i) = Vs_KT;
end

%%
x_shade1 = [0.006 0.006 0.025 0.025];
y_shade1 = [0.5 2.5 2.5 0.5];

x_shade2 = [0.00 0.00 0.05 0.05];
y_shade2 = [1.4 2.0 2.0 1.4];

%%
figure
hold on;
fill(x_shade1, y_shade1, [0.3, 0.5, 0.7], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
fill(x_shade2, y_shade2, [0.8, 0.4, 0.4], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot(P,Vs,'b:o', 'LineWidth', 2, 'MarkerSize', 5,'color','b', 'MarkerFaceColor', 'b');
xticks(0:0.01:0.05);
yticks(0.5:0.5:2.5);
yticklabels((0.5:0.5:2.5)*1000);
xlim([0 0.05]); ylim([0.5 2.5]);
ylabel('S-wave velocity (m/s)'); xlabel('Fracture density');
grid on; box on;
set(gca,'linewidth',2,'fontsize',20,'fontname','Arial','TickLength',[0.007 0.01]);
set(gcf,'unit','centimeters','position',[10,10,15,11]);