clear;clc;
interaction = cell2mat(struct2cell(load(['D:\SZU Learning\project\projectsRBP-AS\RBP_AS_Assciation_code\data\interaction.mat'])));
FS=cell2mat(struct2cell(load(['D:\SZU Learning\project\projectsRBP-AS\RBP_AS_Assciation_code\data\FS.mat'])));
FSP=cell2mat(struct2cell(load(['D:\SZU Learning\project\projectsRBP-AS\RBP_AS_Assciation_code\data\FSP.mat'])));
SS=cell2mat(struct2cell(load(['D:\SZU Learning\project\projectsRBP-AS\RBP_AS_Assciation_code\data\SS.mat'])));
SSP=cell2mat(struct2cell(load(['D:\SZU Learning\project\projectsRBP-AS\RBP_AS_Assciation_code\data\SSP.mat'])));
%interaction:the known RBP-AS event association
% FS:the functional similarity between RBP(i) and RBP(j)
% FSP:Functional similarity weighting matrix
% SS:the module similarity between AS(i) and AS(j).
% SSP:module similarity weighting matrix

SS=SS(4000:4283,4000:4283);
SSP=SSP(4000:4283,4000:4283);
CD44=interaction(:,4000:4283);
y=CD44;
RBP_AS_Y=CD44;
AS_Function_S=SS;
AS_Sem_S=SSP;
RBP_Function_S=FS;
RBP_Sequences_Needle_S=FSP;
% infer the regulating RBPs for CD44
K1 = [];
K1(:,:,1)=RBP_Function_S;
K1(:,:,2)=RBP_Sequences_Needle_S;
K2 = [];
K2(:,:,1)=AS_Function_S;
K2(:,:,2)=AS_Sem_S;
y = RBP_AS_Y;
K1(:,:,3)=kernel_gip(y,1, 1);
K2(:,:,3)=kernel_gip(y,2, 1);
[weight_v1] = FKL_weights(K1,y,1,200);
K_COM1 = combine_kernels(weight_v1, K1);		

[weight_v2] = FKL_weights(K2,y,2,200);
K_COM2 = combine_kernels(weight_v2, K2);
	
score = LapRLS_shan(K_COM1,K_COM2,y, 2^(-5),20,1);
	
	
