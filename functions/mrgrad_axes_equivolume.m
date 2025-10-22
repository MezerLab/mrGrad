function vout = mrgrad_axes_equivolume(mask,varargin)
%--------------------------------------------------------------------------
% function receives ROI mask, number of PC (1, 2, or 3), and number of
% segments and returns the segmented coordinates
%--------------------------------------------------------------------------
% INPUTS:
%
%     mask:   ROI mask or masked qmri map
%
%--------------------------------------------------------------------------
% OPTIONAL INPUTS:
%
%     'PC':       followed by a number of principal component (e.g. 1)
%                 default: 1
%
%     'Nsegs':    followed by a number of wanted segments
%                 default: 7
%
%     'figures':  followed by 1 or 0 in order to produce figures or not.
%                 default: 0
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2021
%--------------------------------------------------------------------------

% specify PC
[found, pc, varargin] = argParse(varargin, 'PC');
if ~found; pc = 1; end
% specify Nsegs
[found, Nsegs, varargin] = argParse(varargin, 'Nsegs');
if ~found; Nsegs = 7; end
% figures?
[found, isfigs, varargin] = argParse(varargin, 'figures');
if ~found; isfigs = 0; end
% pc change axis
[found, maxchange, varargin] = argParse(varargin, 'maxchange');
if ~found; maxchange = []; end
% alternative mask for axes computation
[~, alternative_mask, varargin] = argParse(varargin, 'apply_alternative_axes');

%% Optional: computed axes on a different ROI, then apply on ROI
% If the user chooses to compute axes of a ROI A, then apply it on ROI B:
if ~isempty(alternative_mask)
    % get 3D coordinates of ROI
    [x1,y1,z1] = ind2sub(size(alternative_mask),find(alternative_mask));
    X = [x1,y1,z1];

    % perform PCA
    PCA.coeff = pca(X, 'Algorithm','svd');
    PC = PCA.coeff(:,pc);
end

%% Get 3D coordinates of ROI
[x1,y1,z1] = ind2sub(size(mask),find(mask));
X = [x1,y1,z1];
cent = mean(X);

%% Perform PCA on ROIs coordinates
% If we did not compute axes on a different ROI, we use the default option
% of computing the axes of the actual ROI of analysis
if isempty(alternative_mask)
    PCA.coeff = pca(X, 'Algorithm','svd');
    PC = PCA.coeff(:,pc);
end
% % perform SVD instead of full PCA
% X_centered = X-cent;
% [~,~,V] = svd(X_centered);
% PCA.coeff = V;
% PC = PCA.coeff(:,pc);

%% find data edges
% calculate 2 points on the PC line that are very far from centroid
% (outside data), then find the closest data point for each.

% computing two points on the PC line, outside from data
% line is described by:
% (x(t),y(t),z(t)) = point_on_line + t * direction
% equiv: (x(t),y(t),z(t)) = cent + t * PC'
% equiv:
% x = cent(1) + PC(1) * t
% y = cent(2) + PC(2) * t
% z = cent(3) + PC(3) * t
large_num = 1000;
% px = linspace(cent(1)-large_num,cent(1)+large_num,2); 
% py = cent(2) + (PC(2)/PC(1))*(px-cent(1));
% pz = cent(3) + (PC(3)/PC(1))*(px-cent(1));

pz = linspace(cent(1)-large_num,cent(1)+large_num,2); 
py = cent(2) + (PC(2)/PC(3))*(pz-cent(3));
px = cent(1) + (PC(1)/PC(3))*(pz-cent(3));
% finding a data point closest to each of the line points
data_edges = zeros(2,3);
for pp=1:2
    p = [px(pp),py(pp),pz(pp)];
    p = repmat(p,length(X),1);
    d = rssq(p-X,2);
    [~,mindIDX] = min(d);
    data_edges(pp,:) = X(mindIDX(1),:);
end
%% what is the coordinate that changes the most with the PC axis?
% In order to make sure that the axis segmenting in all subjects is done in the
% "same" direction of their PC lines, we check for each PC with what x,y,z
% axis it changes the most, and then make sure that the axis segmentaytion
% happens from the "lower" part of the PC line to the "higher" part of the
% PC line (w.r.t the coordinate with fastest change).
% since the coordinate with fastest change can differ across subjects, we
% look at all subjects first and then choose the coordinate with fastest
% change in the majority of subjects. so the axis segmentation itself
% happens without correction of direction, then after collecting all
% subjects we flip those that need flipping. this happens in the main
% function 'mrGrad.m', where the function 'RG_flip.m' is called.
%--------------------------------------------------------------------------

% % save axis with fastest change
% [~,max_change_axis] = max(abs([px(1)-px(end),py(1)-py(end),pz(1)-pz(end)]));
%% FIG - Data with PC line and edges
if isfigs
figure;
h2=scatter3(X(:,1),X(:,2),X(:,3),'o');
    h2.MarkerFaceColor = [0, 0.4470, 0.7410];
    h2.MarkerEdgeColor = 'none';
    h2.MarkerFaceAlpha = 0.4;
    h2.MarkerEdgeAlpha = 0.4;
    h2.SizeData = 100;

hold on;
l=line(px,py,pz);
l.LineWidth = 2;
hold on;
for pp=1:2
h1(pp)=plot3(data_edges(pp,1),data_edges(pp,2),data_edges(pp,3),'*r');
hold on;
end
axis equal;
xlim([min(X(:,1))-10,max(X(:,1))+10]);
ylim([min(X(:,2))-10,max(X(:,2))+10]);
zlim([min(X(:,3))-10,max(X(:,3))+10]);
h1(1).MarkerSize = 12;
h1(2).MarkerSize = 12;
h1(1).LineWidth = 2;
h1(2).LineWidth = 2;
xlabel('x','FontSize',22);
ylabel('y','FontSize',22);
zlabel('z','FontSize',22);
view([260,50]);

l=legend('data points','PC1 line','data edges w.r.t PC1');
l.FontSize = 22;
end
%% intersection of outer planes with PC line
l_start = intersec_plane_line(PC, data_edges(1,:), cent);
l_end = intersec_plane_line(PC, data_edges(2,:), cent);

% if axis change direction provided, use provided
if ~isempty(maxchange) &&  l_start(maxchange(pc)) > l_end(maxchange(pc))
    tmp = l_start;
    l_start = l_end;
    l_end = tmp;
end
%% Define new points of 1st PC line: all intersections with planes
% % line is described by:
% % (x(t),y(t),z(t)) = point_on_line + t * direction
% % equiv: (x(t),y(t),z(t)) = cent + t * PC'
% % equiv:
% % x = cent(1) + PC(1) * t
% % y = cent(2) + PC(2) * t
% % z = cent(3) + PC(3) * t

x1 = l_start(1);
x2 = l_end(1);
% x1 = min(l_start(1), l_end(1));
% x2 = max(l_start(1), l_end(1));
lx = [x1,x2];%linspace(x1,x2,Nsegs+1); % Nsegs+1 points on PC line
ly = cent(2) + (PC(2)/PC(1))*(lx-cent(1));
lz = cent(3) + (PC(3)/PC(1))*(lx-cent(1));
%% FIG - scatter with segmented PC line
if isfigs
figure;
p1=scatter3(X(:,1),X(:,2),X(:,3),'o');
hold on;
l1=line(lx,ly,lz);
hold on;
p2=plot3(lx,ly,lz,'*g');
hold off
figproper(l1);
    p1.MarkerFaceColor = [43,140,190]/255;
    p1.MarkerEdgeColor = 'none';
    p1.MarkerFaceAlpha = 0.4;
    p1.MarkerEdgeAlpha = 0.4;

% p1.Color = [43,140,190]/255;
p2.Color = [254,178,76]/255;
p2.MarkerSize = 12;
p2.LineWidth = 2;
end
%% compute all planes
% Plane equation: PC(1)*(x-p(1))+PC(2)*(y-p(2))+PC(3)*(z-p(3)) = 0
gmax = 1000;
gmin = -gmax;
int = (gmax-gmin)/10;
[gx, gy] = meshgrid(gmin:int:gmax);
gz = cell(Nsegs+1,1);
for seg = 1:2 %:Nsegs+1
    p = [lx(seg) ly(seg) lz(seg)];
    gz{seg} = -(PC(1)*(gx-p(1))+PC(2)*(gy-p(2)))/PC(3) + p(3);
end
%% FIG
if isfigs
figure;
p1=scatter3(X(:,1),X(:,2),X(:,3),'o');
hold on;
l1=line(lx,ly,lz);
hold on;
p2=plot3(lx,ly,lz,'*g');
for seg=1:Nsegs+1
    hold on;
    surf(gx, gy, gz{seg},'FaceAlpha',0.5, 'EdgeColor', 'none');
end
hold off;
figproper(l1);
    p1.MarkerFaceColor = [43,140,190]/255;
    p1.MarkerEdgeColor = 'none';
    p1.MarkerFaceAlpha = 0.4;
    p1.MarkerEdgeAlpha = 0.4;

% p1.Color = [43,140,190]/255;
p2.Color = [254,178,76]/255;
p2.MarkerSize = 12;
p2.LineWidth = 2;
xlim([min(X(:,1))-10,max(X(:,1))+10]);
ylim([min(X(:,2))-10,max(X(:,2))+10]);
zlim([min(X(:,3))-10,max(X(:,3))+10]);
end
%% calculate distance of data points from planes, classify to segments
% PCs are unit normal vectors, so their length is 1 (rssq(PC) = 1);

% define number of voxels in each segment (should be equal, except the last one
% that can be a little different)
numVox = round(length(X)/Nsegs);
numVox_last = length(X) - (Nsegs-1)*numVox;
numVox = [repmat(numVox,Nsegs-1,1); numVox_last];

% find numVox data points that are closest to 1st plane.
% the distance between data point and plane is given by:
% distance_eq = (on_plane_point - data_point) * PC;
tmp_data = X;
segment_data = cell(Nsegs,1);
p_onplane = repmat(l_start,length(tmp_data),1); % point on plane
norm_vec = repmat(PC',length(tmp_data),1);
distances = abs(sum((p_onplane-tmp_data).*norm_vec, 2));
[a,b] = sort(distances);
for j=1:Nsegs
    pntsIDX = b(1:numVox(j));
    segment_data{j} = tmp_data(pntsIDX,:);
    b(1:numVox(j))=[];
end
%% linearizing coordinates
linearInd = cell(Nsegs,1);
for j=1:Nsegs
    b = segment_data{j};
    linearInd{j} = sub2ind(size(mask), b(:,1),b(:,2),b(:,3));
end
%% FIG
if isfigs
figure;
cmap = [178,24,43; 214,96,77; 244,165,130; 253,219,199; 209,229,240; 146,197,222; 67,147,195; 33,102,172]/255;
h = cell(Nsegs,1);
for j=1:Nsegs
    inpoints = segment_data{j};
    h{j}=scatter3(inpoints(:,1),inpoints(:,2),inpoints(:,3),'o');
        h{j}.MarkerFaceColor = cmap(j,:);
        h{j}.MarkerEdgeColor = 'none';
        h{j}.MarkerFaceAlpha = 0.4;
        h{j}.MarkerEdgeAlpha = 0.4;
        hold on;
end
l1=line(lx,ly,lz);
for seg=1:Nsegs+1
    hold on;
    surf(gx, gy, gz{seg},'FaceAlpha',0.5, 'EdgeColor', 'none');
end
hold on;
plot3(lx(j),ly(j),lz(j),'*y')
hold on;
plot3(PC(1),PC(2),PC(3),'*g')
hold off;
figproper(l1);
xlim([min(X(:,1))-10,max(X(:,1))+10]);
ylim([min(X(:,2))-10,max(X(:,2))+10]);
zlim([min(X(:,3))-10,max(X(:,3))+10]);
h1=title(['Segmenting ROI along PC',num2str(pc),' axis, single subject']);
h1.FontSize=22;
end
%% OUTPUT
vout.all_Inds = X;
vout.segment_inds = segment_data;
vout.segment_inds_linear = linearInd;
vout.N_segments = Nsegs;
vout.segment_size = [];
vout.segment_Nvoxels = numVox;
vout.planes.gx = gx;
vout.planes.gy = gy;
vout.planes.gz = gz;
vout.PC_coeff = PCA.coeff;
vout.analysisPC = pc;
% vout.max_change_axis = max_change_axis;
vout.analysisPC_line = [lx; ly; lz];
vout.data_centroid = cent;
end
%% FUNCTIONS
function interpoint = intersec_plane_line(PC, p, cent)
% function recive a point p in which a plane orthogonal to PC goes
% through, and returns the intersection point of the plane and the line of PC

% line equation
syms t
x = cent(1) + PC(1) * t;
y = cent(2) + PC(2) * t;
z = cent(3) + PC(3) * t;

% plane equation
eq = PC(1)*(x-p(1))+PC(2)*(y-p(2))+PC(3)*(z-p(3)) == 0;
t = double(solve(eq, t));
interpoint = cent + t*PC';
end

function figproper(l1)
xlabel('x','FontSize',22);
ylabel('y','FontSize',22);
zlabel('z','FontSize',22);
axis equal;
if exist('l1','var')
    l1.LineWidth = 2;
    l1.Color = [254,178,76]/255;
end
end