function [rf_results] = rf_naive(Xtest, Xtrain, ytrain)


    % Input number of trees in the forest
    nTrees = 50;

    % Train TreeBagger (Decision Forest)
    tic;
    rf_model = TreeBagger(nTrees, Xtrain, ytrain,...
        'Method', 'classification',...
        'OOBPred','On');
    toc;


    % Preliminary Stuff
    oobErrorBaggedEnsemble = oobError(rf_model);
    figure;
    plot(oobErrorBaggedEnsemble);
    xlabel 'Number of grown trees';
    ylabel 'Out-of-Bag Classification Error';

    % Apply model to test data
    tic;
    rf_results = rf_model.predict(Xtest);
    toc;
    
    rf_results = cell2mat(rf_results);
    rf_results = str2num(rf_results);

end