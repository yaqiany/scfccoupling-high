%% load data
% Lad `Mats` cell array containing group-consensus structural (sc) and
% functional (fc) networks and 3D coordinates. For more information on methods used to generate group-consensus
% networks see Betzel et al., (2019), Net Neurosci
% clear all;clc;
load('Consensus_Connectomes.mat');
Mats=LauConsensus.Matrices;
nscale = 5;                  % 5 parcellation scales for the Lausanne atlas
ii=5;                        % select the highest resolution(1000*1000)
W= Mats{ii, 1};              % structure weighted matrix
Sc= Mats{ii, 1};Sc(Sc>0)=1;  % structure binary matrix
coor = Mats{ii, 4};          % x,y,z node coordinates | nx3 matrix
N=size(W,1);                 % node number
fc = Mats{ii, 3};            % group-consensus resting-state functional network | nxn node matrix

A=W;A=-A;
for i=1:N
A(i,i)=-sum(A(i,:));
end
B=A/max(eig(A));   % structural Laplacian matrix
[BEC,BE]=eig(B);   % BEC-eigenvectors, BE-eigenvalues
%将结构矩阵的特征值�?特征向量按从小到大排�?
[a,b]=sort(diag(BE));
BE_sort=diag(a);
BEC_sort=BEC(:,b);

mappedX = diffusion_maps(fc, 10, 0.5, 1);   % diffusion map algorithm on functional network
x = mappedX(:, 1);                          % node weights for first eigenvector/component
xx = x * -1;                                 % reverse weights such that positive = top of hierarchy
