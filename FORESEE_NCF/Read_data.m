clc, clear, close all;

%% load data
files = dir('*.mat');

numArray = zeros(numel(files), 1);
for i = 1:numel(files)
    numberPart = regexp(files(i).name, '(\d+)', 'match');
    if ~isempty(numberPart)
        numArray(i) = str2double(numberPart{1});
    else
        numArray(i) = 0;
    end
end

[~, idx] = sort(numArray);
files = files(idx);

%%
dx = 2;
dt = 0.008;

for k = 1:length(files)
    fileName = files(k).name;  
    fprintf('Loading file: %s\n', fileName);
    
    data0 = load(fileName);
    seis = data0.data';
    [nz,nx] = size(seis);

    figure
    imagesc(seis);
    colormap(gray);
    xlabel('Channel'); 
    ylabel('Time (s)');
    
    title_name = strcat('No.',fileName);
    title(title_name(1:length(fileName)-1));
    
    box on;
    clim([-0.01 0.01]);
    yticks(1:50:nz);
    yticklabels(((1:50:nz)-1)*dt);
    
    set(gca,'Linewidth',2,'fontsize',20,'Fontname','Arial');
    set(gcf,'unit','centimeters','position',[10,10,15,15]);
    drawnow;
end

