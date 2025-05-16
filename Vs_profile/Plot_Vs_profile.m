 clc;clear;close all;

%%
Location = load('foresee_ch_loc.txt');
V = load('Vs_profile.mat');
V_smooth = V.V_smooth;

%%
sx = zeros(1,2137);
sy = zeros(1,2137);
for i = 1:2137
    for j = 1:150-1
        if isnan(V_smooth(j,i)) && ~isnan(V_smooth(j+1,i))
            sx(i) = i; sy(i) = j;
        end
    end
end

figure
h = imagesc(V_smooth);
hold on;
scatter(sx,sy,10,'k','filled');
set(h,'alphadata',~isnan(V_smooth));
colormap(flipud(color_map('spectral')));
yticks(-13:160/8:147);
yticklabels(((-13:160/8:147)+303-350)*-1+350);
ylim([27 147]);
clim([1100 2500]);
h = colorbar; set(h, 'YDir', 'reverse');
h.Label.String = 'S-wave velocity (m/s)';
h.YTick = (1100:400:2500);
xlabel('Channel'); ylabel('Elevation (m)');
grid on; box on;
set(gca,'linewidth',2,'fontsize',20,'fontname','Arial','TickLength',[0.007 0.01]);
set(gcf,'unit','centimeters','position',[0,15,45,10]);

%%
X = Location(:,2)';
Y = Location(:,1)';
X = repmat(X, 150, 1);
Y = repmat(Y, 150, 1);
Z = (1:150)';
Z = repmat(Z, 1, 2137);

xx = X(:);
yy = Y(:);
zz = Z(:);
vv = flipud(V_smooth);
vv = vv(:);

%%
figure
scatter3(xx,yy,zz,50,vv,'filled');
colormap(flipud(color_map('spectral')));
clim([1100 2500]);
view([-80 40]);
box on;
zticks(3:160/8:163);
zticklabels(((147:-160/8:-13)+303-350)*-1+350);
zlim([3 123]);

xlim([-77.880 -77.850]);
xticks(-77.880:0.01:-77.850);
ylim([40.788 40.815]);
yticks(40.788:0.004:40.815);

zlabel('Elevation (m)');
h = colorbar; set(h, 'YDir', 'reverse');
h.Label.String = 'S-wave velocity (m/s)';
h.YTick = (1100:400:2500);
set(gca,'linewidth',2,'fontsize',20,'fontname','Arial','TickLength',[0.007 0.01]);
set(gcf,'unit','centimeters','position',[0,2,50,15]);
