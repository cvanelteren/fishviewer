% Load dat set
clear all;

% set paths; assuming im running from dd-CRP master
addpath('../ZebrafishRecordings/cell_view_package'); % scripts for i/o datasets
datasetPath = '../ZebrafishRecordings/photo_dataset1'                       % change for indexing into datasets
addpath(datasetPath)

load(strcat(datasetPath, '/cell_resp_dim_lowcut.mat'));
load(strcat(datasetPath,'/cell_info.mat'));

coordinates = zeros(cell_resp_dim(1), 3); % x,y,z
for i = 1: length(cell_info)
    xy    = cell_info(i).center;
    z     = cell_info(i).slice;      % assuming this is the z-coordinate
    coordinates(i,:) = [xy, z];
end

% visualize the coordinates
% figure()
% scatter3(coordinates(:,1), coordinates(:,2), coordinates(:,3));
% z-coordinates need cleaning, not clear as to what belongs where
%%
dims = max(coordinates);
% reduce the data by super-voxels
voxelDims           = [16, 16, 5];
newCoordinates      = zeros(size(coordinates));           %empty storage
grid = {};
for n = 1: numel(voxelDims)
    grid{n}              = 0 : voxelDims(n) : dims(n) + voxelDims(n);        %create bins
    %     grid{n}                = linspace(0, dims(n), dims(n) / voxelDims(n));
    newCoordinates(:, n) = discretize(coordinates(:, n), grid{n},'IncludedEdge', 'right');
end
%%
data = read_LSstack_fast_float(strcat(datasetPath, '/cell_resp_lowcut.stackf'), cell_resp_dim);
%%
tmp = {};
i = 0;
uniques = unique(newCoordinates, 'rows');
for xyz = uniques'
    i      = i + 1;
    member = find(ismember(newCoordinates, xyz', 'rows'));
    tmp{i} = mean(data(member, :),1);
end
data = cell2mat(tmp');      % empty rawData
% data = reshape(data, i, cell_resp_dim(2));
data = zscore(data, 0, 2);
tmp  = {};                  % empty tmp
%%
% The subset of the data is selected according to the top 5 percent
% of the correlation coefficients with the behavior vector
% The data was smoothed with a gaussian kernel (sigma = 10)


% steps :
%   First we need to construct A;
%     In the demo script provided this is represented as a sparse code in a
%     struct file structure.
%   Next I will visualize the data

nCell       = size(data, 1);
% coordinates = [xx;yy;zz]; % x y need to be shifted to make it in same plane
% minC        = min(coordinates');

% reset the coordinates to zero minimum
% for ci = 1 : size(coordinates,1)
%     coordinates(ci,:) = coordinates(ci,:) - minC(ci);
% end
% only use the ones from the extracted cells
% redCoords = coordinates(:, indices)
% redCoords   = coordinates' + 1;
newCoordinates = uniques;
% Create weight matrix to remove any cells that are far way
% weights                          = 1./dist(redCoords', redCoords, 'euclidean');
weights                          = squareform(1 ./ pdist(newCoordinates,'euclidean'));
weights(eye(size(weights)) == 1) = nan;       % remove diagonal
% diag(weights)
weights                          = num2cell(weights, 2); %cellarlize
%%
A  = {1:nCell};
A  = A(ones(nCell,1)) ;
% create the Adjacency matrix
for i = 1:nCell
    % create space code of connection to connection
    indx = weights{i} < 1;
    A{i}(indx) = [];
    weights{i}(indx) = [];  %  remove if not adjacent
end

% weights = {};
size(weights)
%%
% from the provided script (not altered)
opts.steps      = 20; % default 20 ; reduced for time
opts.hyp.a0     = 2;
opts.hyp.b0     = 1;
opts.hyp.mu0    = 0;
opts.hyp.kappa0 = 1;

% run the algorithm
[MAP, samples]  = PMC_ddCRP_NG(data',A,opts); % data clustered along the 2nd dim



%% visualize the results
% clear all;
% load cleanedDataClustered
% load clusterRedData
Z = double(bsxfun(@eq, MAP.Pi, 1:max(MAP.Pi)));
Q = Z*Z';
nClus        = max(MAP.Pi);
clusterSizes = zeros(nClus,1);
clusters     = zeros(size(Z,1),1);
threshold    = 150;
% thresholded  = zeros(size(Z,1),1);
j = 1;
for i = 1 : nClus
    findThem                = find(MAP.Pi == i);
    clusterSizes(i)         = numel(findThem);
    clusters(findThem)      = i;
end

figure()
histogram(clusterSizes,10);
xlabel('Cluster size')
ylabel('Count');

%% plot only large clusters
close all
threshold = 150;
biggest   = find(clusterSizes > threshold);
n = numel(biggest);
c         = colorcube(n+5); % design a high contrasting colormap
figure()
hold on;
% j = 1;

x = double(newCoordinates(:,1));
y = double(newCoordinates(:,2));
z = double(newCoordinates(:,3));
k = boundary(x,y,z);
for i = 1: n
    idx = find(Z(:, biggest(i)));
    colorToPlot = c(i, :);
    
    scatter3(x(idx),...
        y(idx),...
        z(idx), 10, colorToPlot,'filled');
    %     trisuf(k, x,y,z)
    %     scatter3(x,y,z,.1,[.9,.9,.9])
    %     scatter(coordinates(idx,1), coordinates(idx, 2), 'filled')
end
colormap(c)
xlabel('x'); ylabel('y'); zlabel('z');
xlim([1, max(x)]); ylim([1, max(y)]); zlim([0, max(z)])
val  = max(coordinates,[],1);
colorbar(gca)
title('Empirical found cluster assignment')


%% make a time lapse of cluster activity;
% per cluster compute the min max and mean of activity over time
% modulate the alpha (transparancy) of the color of the cluster per time
% step
% met dif zou je de grens kunnen opzoeken door eerst een map te maken per
% cluster en die te vullen met 1 en dan dif te doen dan krijg je de grens
%%
load(strcat(datasetPath, '/frame_turn'));
% close all
behavior = frame_turn(1 : size(data,2), 17);
behavior = smooth(behavior);
c  = corr(MAP.ClusterTCs, behavior);
[sortC, idxC] = sort(abs(c), 'descend' ); 
plotTop       = 20;
selection     = idxC(1:plotTop);

% selection = find(clusterSizes > 20 & clusterSizes < 150 );
% selection = find(clusterSizes > 0);


clusterData = MAP.ClusterTCs(:, selection);
relData     = [min(clusterData,[], 1); median(clusterData,1); max(clusterData,[],1)];
nT          = size(clusterData,1);
nC          = size(clusterData,2);
dims        = max(newCoordinates,[],1);

tmp = nan(dims(1) + 1, dims(2) + 1);
nZ  = numel(unique(z));

c         = colormap(colorcube(numel(selection)));  % outline colors of relevant clusters

% Visualize different layers with behavior 
f = figure(1);      
clim = [-2,2];                % keep consistent color values
h    = {};                      % storage for image of layers

hold('on');

% setup the figure; for every layer make axes update is below
for layer = 1 : nZ
    subplot(4,2,layer);
    hold('on');
    xlim([1, dims(1)+1]); ylim([1, dims(2)]);
    h{layer} = imagesc(tmp',clim);
    colormap(parula);
    for ci =  1 : numel(selection)
        member = find(Z(:, selection(ci)));
        if ~isempty(member)
            disp(numel(member));
            xx = double( x(member) );
            yy = double( y(member) );
            k  = boundary(xx,yy);                           % boundary of clusters
            
            plot(xx(k), yy(k), '-.',...
                'color',c(ci,:), 'linewidth', 1,...
                'markersize', 2);
        end
    end
end

subplot(4,2,7:8)                                            % behavior plot (animated)
hold('on');                                                 % matlab shenaningans
plot(behavior,'-b');
% hold('on');
hh = plot([0,0],[min(behavior), max(behavior)], '--r');

ylim([min(behavior), max(behavior)]);
xlim([0, length(behavior)])
% hold('off');

disp(size(h));
% hold('on');
% hold('off')
% time animation
% for every time point
for t = 1 : nT
    % update layers
    for layer =  1 : nZ
        idx = find(z == layer);
        dd = h{layer}.CData;
        % update voxels
        for i = idx'
            dd(y(i), x(i)) = data(i, t);
        end
        h{layer}.CData = dd;            % set data
    end
    %     subplot(4,2,7:8)
    hh.XData = [t,t];
    suptitle(sprintf('t = %d', t));
    drawnow()                           % update figure
end
% hold('off');
% disp('done'


%% relate to behavior vector
% load behavior_vector
lr = frame_turn(1:size(MAP.ClusterTCs,1),17);  % I think the first contains left right

c  = corr(MAP.ClusterTCs, behavior);
fig = figure();
hist(c)
xlabel('R'); ylabel('Counts')
title('Cluster time courses behavior correlation')
saveas(fig, 'Corr_clus_beh_hist.jpg')

% plot the largest R with behavior time courses
[sortC, idxC] = sort(abs(c), 'descend')
plotTop = 10;

fig = figure();
hold('on');
plot(lr,'k--');
labels = {'Behavior'}
for idx = 1 : plotTop
    plot(MAP.ClusterTCs(:,idxC(idx)))
    labels{idx + 1} =  sprintf('R = %.3f', c(idxC(idx)));
end
legend(labels);
xlabel('Time [step]')
hold('off')
saveas(fig, 'Corr_clus_beh.jpg')

% visuzalize the clusters positions in the brain (x,y) projection
cmap = parula(plotTop);

fig = figure();
set(gca,'Color',[0.8 0.8 0.8]);
hold('on');

% uz = unique(z);
% for i = uz'
%     idx = find(z == i);
%     tx  = x(idx);
%     ty  = y(idx);
%     tz  = z(idx);
%     k = boundary(tx,ty);
%     plot3(tx(k), ty(k), tz(k),'k')  % plot the map;
% end
scatter3(x,y,z, 3, 'k');
for idx = 1 : plotTop
    members = find(Z(:, idxC(idx)));
    scatter3(x(members), y(members),z(members), 35, cmap(idx,:), 'filled');
end
labels{1} = 'Outline';
legend(labels)
xlabel('x'); ylabel('y');
hold('off')
saveas(fig, 'Corr_clus_beh_map.jpg')

%%
% threshold = 100;
% 
% 
% layer = 2;
% 
% selection = find(clusterSizes > 20 & clusterSizes < 150 );
% selection = find(clusterSizes > 50);
% 
% 
% clusterData = MAP.ClusterTCs(:, selection);
% relData     = [min(clusterData,[], 1); median(clusterData,1); max(clusterData,[],1)];
% nT          = size(clusterData,1);
% nC          = size(clusterData,2);
% 
% dims        = max(newCoordinates,[],1);
% 
% 
% % plotMatrix  = zeros(dims(1), dims(2));
% c           = colorcube(nC + 2);
% figure(2)
% tmp = nan(dims(1), dims(2));
% hold('on');
% 
% h = imagesc(tmp, [-2,2]);
% 
% colorbar;                                                       % colorbar
% xlabel('x'); ylabel('y')                                        % axis labels
% axis square
% 
% cLayerCounter = 0;
% % draw cluster boundaries
% for ci = 1: nC
%     member = find(Z( :, selection(ci)));
%     %         member = find(Z(z == layer, selection(ci)));
%     if ~isempty(member)
%         disp(numel(member))
%         xx = double( x(member) );
%         yy = double( y(member) );
%         k  = boundary(xx,yy);                           % boundary of clusters
%         
%         %         plot(xx(k)-1, yy(k)-1, '.',...
%         %             'color',c(ci,:), 'linewidth', 3,...
%         %             'markersize', 10);
%         plot(xx, yy, '.',...
%             'color',c(ci,:), 'linewidth', 3,...
%             'markersize', 10);
%         cLayerCounter = cLayerCounter + 1;
%         
%     end
%     %     pause(1)
%     
% end
% fprintf('Number of clusters %d\n',cLayerCounter);
% xlim([1, dims(1)+1]); ylim([1, dims(2)]);
% 
% % loop over time points and visualize the data with the boundaries
% % visualized
% % idx = find(z)'
% 
% idx = find(z ==2);
% for t = 1 : nT
%     for tt = idx'
%         tmp(x(tt),y(tt)) = data(tt,t);
%     end
%     h.CData = tmp';
%     
%     title(sprintf('t=%d',t))
%     drawnow()
%     %     M(t)    = getframe(gcf);
% end
% hold('off')

