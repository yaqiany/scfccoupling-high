%% Decoupling
% Structure-function divergence across large-scale functional network gradients based
% on low frequency eigenmodes.

%% prediction model
% A multilinear model was used to predict the functional connection profile
% of every node based on the 2nd-14th eigenmodes of structural Laplacian
% network.
rsq = zeros(N, 1); % node-wise  vector
% define fc response (y) and sc predictors (x)
for jj = 1:N    
    y = fc(:, jj);   
    x1 = zscore(BEC_sort(:,2:14));  % standardize predictors of 2nd-7th eigenmodes    
    % fit multiple regression (OLS, main effects only) exclude self-connections for all variables. N.B., 
    lm = fitlm(x1, y, 'Exclude', jj);
    rsq(jj) = lm.Rsquared.Ordinary; % record R-square for parcellation ii, node jj
end
rsq1=sqrt(rsq);

%% Fig2: Boxplot of region 
% plot the boxplot of R across regions where the region is sorted by mean R
load('rsn_mapping.mat');
ii=5;N=1000;
R_name=rsn_names;                         % total region name 
R_n=length(R_name);                       % total region number
Rn_Ind=rsn_mapping{ii};                   % R_name index for each node
Rn_name=cell(N,1);                        % R_name for each node
R_R=zeros(R_n,1);                         % region median R
R_A=zeros(R_n,1);                         % region average R
for i=1:R_n
    R_R(i)=median(rsq1(Rn_Ind==i));
    R_A(i)=mean(rsq1(Rn_Ind==i));
    Rn_name(Rn_Ind==i)=R_name(i);
end

[temp,I]=sort(R_R,'descend');             % sort region by R mean
a=rsq1;b=Rn_name;
c=R_name(I);                              % sorted region name                    
Color=[219 2 10;231 95 27;238 146 43;246 191 65;246 236 84;202 222 169;147 205 137;76 177 99]/255;
Color=flipud(Color);
%Color=Color(I);d=Color(Other_O,:);
figure;boxplot(a,b,'Orientation','horizontal','GroupOrder',c,'Colors',Color,'BoxStyle','filled','Whisker',0,'OutlierSize',2.5,'Symbol','.','MedianStyle','target');


%% Fig2-1: R_region compared to null
% calculate 'region2num' in session 'Fig2'
ci=Rn_Ind;
% number of communities / classes
nci = max(ci);
% label-permuting null model, set number of permutations
nperm = 10000;
% initialize region-wise mean R-square for each permutation
ci_perm = zeros(nci, nperm);
for i = 1:nperm
    % permute community assignments (without replacement)
    p = randperm(N);
    % dummy-code permuted community assignments
    dumdum = dummyvar(ci(p));

    % get mean R-square for permuted community assignments
    ci_perm(:, i) = rsq1' * dumdum ./ sum(dumdum);
end

R_R1=R_A(I);
ci_perm1=ci_perm(I,:); %ci_perm1(7,10000)
figure;
edges= 0.4:0.01:0.6;
for i=1:7
subplot(2,4,i)
h=histogram(ci_perm1(i,:),12,'Normalization','probability');
h.FaceColor = '#0072BD';
h.EdgeColor = '#0072BD';
set(gca, 'FontSize',12 ); 
set(gca, 'Fontname','Times New Roman' ); 
hold on;
plot([R_R1(i) R_R1(i)],[0 0.28],'LineWIdth',2);
end


%% Fig2-3: plot pie of high accuracy nodes proportion across regions
a=ci(rsq1>mean(rsq1));        % find nodes with prediction accuracy higher than average level
ord_R=I;                      % region order(low-frequency)
cm=[0 65 133;73 98 153;0 122 184;69 147 199;102 175 212;153 202 227;181 216 239;178 226 241;237 246 251]/255;
figure;dumdum = dummyvar(a);s_dum=sum(dumdum);
labels=R_name(I);             % region order
subplot 121
pie(s_dum(ord_R));
subplot 122
pie(s_dum(ord_R),labels);colormap(cm);
cm=[0 144 72;64 172 118;128 200 163;157 186 96;186 213 133;216 217 109;229 230 157;242 242 206;249 249 231]/255;
figure;dumdum = dummyvar(ci);s_dum=sum(dumdum);
subplot 121
pie(s_dum(ord_R));
subplot 122
pie(s_dum(ord_R),labels);colormap(cm);

%% Fig2-4: region map(z-score compared to node label shuffle)
dumdum = dummyvar(ci);
tmp=(repmat(rsq1,1,R_n)-mean(ci_perm'))./std(ci_perm');
z_rsq=tmp.*dumdum;
z_rsq1=sum(z_rsq,2);

a=[z_rsq1,coor,dumdum];
a=sortrows(a,1);
% colormap
mycolorpoint=[[48 122 184];[71 146 206];[189 216 238];[255 255 255];[248 180 182];[239 100 103];[226 46 46]]/255; % blue to red
mycolorposition=[1 400 466 533 592 652 N];
mycolormap_r=interp1(mycolorposition,mycolorpoint(:,1),1:N,'linear','extrap');
mycolormap_g=interp1(mycolorposition,mycolorpoint(:,2),1:N,'linear','extrap');
mycolormap_b=interp1(mycolorposition,mycolorpoint(:,3),1:N,'linear','extrap');
C=[mycolormap_r',mycolormap_g',mycolormap_b'];%color
% figure
for i=1:R_n
figure
region_i=i;
x=a(:,2);y=a(:,3);z=a(:,4);
S=ones(N,1)*100;                    % size
C1=ones(N,3)*175/255;               % all nodes are grey
C2=C(a(:,region_i+4)>0,:);
x2=a(a(:,region_i+4)>0,2);y2=a(a(:,region_i+4)>0,3);z2=a(a(:,region_i+4)>0,4);
S2=S(a(:,region_i+4)>0);
scatter3(x,y,z,S,C1);hold on;
scatter3(x2,y2,z2,S2,C2,'filled');
colormap(C);title(R_name(i));
%axis off 
end

%% Figure 2-4
% The spatial distribution of R. Nodes are colored in proportion to
% R; nodes with weak and high structure-function correspondence are colored by
% blue and red respectivily. High correspondence is observed in primary
% sensory and motor cortices, while lower correspondence is observed in
% transmodal cortex

figure
x=coor(:,1);y=coor(:,2);z=coor(:,3);
N1=ceil(max(rsq1)*10000)-floor(min(rsq1)*10000);
S=ones(N,1)*200;
mycolorpoint=[[255 135 158];[255 191 204];[250 235 214];[239 255 255];[64 104 224];[20 48 135]]/255;%color point red-blue
mycolorpoint=flipud(mycolorpoint); 
mycolorposition=[1 1000 round(N1/2)-100 round(N1/2)+100 N1-1000 N1];%rsql
mycolormap_r=interp1(mycolorposition,mycolorpoint(:,1),1:N1,'linear','extrap');
mycolormap_g=interp1(mycolorposition,mycolorpoint(:,2),1:N1,'linear','extrap');
mycolormap_b=interp1(mycolorposition,mycolorpoint(:,3),1:N1,'linear','extrap');
C=[mycolormap_r',mycolormap_g',mycolormap_b'];colormap(C);
scatter3(x,y,z,S,C(round((rsq1-min(rsq1))*10000)+1,:),'filled')
axis off;colorbar('Ticks',[])
%% Figure 3-1
% The figure shows the node-wise R values are anticorrelated with
% position along functional gradient, suggesting that structure and
% function closly correspond in unimodal cotex but diverge in transmodal
% cortex
figure;
y = rsq1;
[rho, pval] = corr(xx, y); 
lm = fitlm(xx, y);
xhat = linspace(min(xx), max(xx), 100);
yhat = lm.Coefficients.Estimate(1) + (lm.Coefficients.Estimate(2) * xhat);
plot(xx, y, '.','Markersize',10); hold on
plot(xhat, yhat,'LineWidth',2)
title(['rho = ' num2str(rho) ', p = ' num2str(pval)]);
axis square
xlabel('gradient');ylabel('R');

%% Reference
% [1] Situating the default-mode network along a principal gradient of macroscale cortical organization
%
% [2] Gradients of structure-function tethering across neocortex
%
% [3] Hierarchical Connectome Modes and Critical State Jointly Maximize Human Brain Functional Diversity
% 
% [5] Mapping functional brain networks from the structural connectome Relating the series expansion and eigenmode approaches

