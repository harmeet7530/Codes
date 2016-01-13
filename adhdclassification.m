%% ----------------------------------------------------
% Author : Harmeet Singh
% Title : Machine Based Classification of ADHD patients using Graph
% Theoretical Approach
% Date : December 10, 2015

clear all
close all


%% Feature based - Clustering Coefficient

clear all
close all

%Define directories
mainDir = 'C:\Users\Harmeet\Documents\ADHD200_CC200\';
temp = dir(mainDir);
temp = temp(3:2082,1);

i=1:4:numel(temp);
subjects = temp(i);
fvector = [];

% Feature Extraction using Brain Connectivity Toolkit
for i=1:numel(subjects)
    
    cmatrix = load(subjects(i).name);
    
    cmatrix = abs(cmatrix);

    C = clustering_coef_bd(cmatrix);

    feature = [C'];
    
    fvector(i,:) = feature;
end

save('fvector_C.mat','fvector');


%% Normalization of Feature Vector
fvector_norm = (fvector-min(fvector(:)))/...
    (max(fvector(:))-min(fvector(:)));


%% Load Ground Truth

GT = load('classes.txt');

%% Binary Classification of Ground Truth
idx = find(GT ~= 1);
GT(idx) = 2;

TrainSubjects = fvector_norm(1:100,:);
TestSubjects = fvector_norm(101:520,:);

TrainGT = GT(1:100,1);
TestGT = GT(101:520,1);

numLabels = max(GT);

figure; hist(TrainGT)
title('Training Ground Truth')
figure; hist(TestGT)
title('Testing Ground Truth')




%% Random Forest

rf_results = rf_naive(TestSubjects,...
    TrainSubjects, TrainGT);

figure;
hist(rf_results)
title('RF Results')

%% Classificaton Tree

tree = ClassificationTree.fit(TrainSubjects,TrainGT);
maxprune = max(tree.PruneList);
treePrune = prune(tree,'level',maxprune-3);
view(treePrune,'mode','graph');
pred = predict(tree,TestSubjects);
subplot(2,1,1);
hold off; plot(TestGT,'g','linewidth',2); hold on; plot(pred,'b','linewidth',2);
subplot(2,1,2);
plot(TestGT-pred,'k','linewidth',3);
mse = (1/length(TestGT))*sum((TestGT-pred).^2);  %0.14
fprintf('MSE for regular classification tree: %f \n ',mse);