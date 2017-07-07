dataset = loadData();

fminBest = cell(1, length(dataset));
AlgoNames = {'FSOLVE'; 'LSQNONLIN'; 'LSQNONLIN_LM'; 'FMINUNC'; 'FMINUNC_QN'; 'FMINCON'; 'FMINCON_INTERIOR'; 'FMINCON_SQP'; 'FMINCON_ACTIVESET'};
FunNames = {'_SOLVE', '_RUNGE', };
RowNames = cell(1, length(AlgoNames)*length(FunNames));
k = 1;
for j = 1:length(FunNames)
    for i=1:length(AlgoNames)
            RowNames{k} = [AlgoNames{i} FunNames{j}];
            k = k + 1;
    end
end
for testNum = 1:length(dataset)
    kFirstC = squeeze(testRes{testNum}.kFirst(1,:,:))';
    kFirstGOFC = squeeze(testRes{1}.kFirstGOF(1,:,:))';

    kFirstR = squeeze(testRes{testNum}.kFirst(2,:,:))';
    kFirstGOFR = squeeze(testRes{testNum}.kFirstGOF(2,:,:))';

    dataC = [kFirstC kFirstGOFC];
    dataR = [kFirstR kFirstGOFR];
    data = [dataC; dataR];
    data(:, end+1) = sqrt(abs(data(:,2).*data(:,3)));
    tt = array2table(data, 'RowNames', RowNames, 'VariableNames',{'K_max','K_s','K_I', 'r2', 'rmse', 'c_zv'});
    display('------------------------------------------------------------------------------------------');
    display(['                                   ' dataset{testNum}.title]);
    tt=sortrows(tt,-4);
    
    fminBest{testNum} = tt(1:10,:);
    display(tt(1:10,:));
    
    if testNum == 1
       winMap = containers.Map(tt.Properties.RowNames, 1:length(RowNames));
    else
        for i = 1:length(RowNames)
            winMap(tt.Properties.RowNames{i}) = winMap(tt.Properties.RowNames{i}) + i;
        end
    end
    
  
    
end

winT=array2table(cell2mat(winMap.values)' ./ length(RowNames), 'RowNames', winMap.keys', 'VariableNames', {'Scores'});
winT = sortrows(winT,1);
disp('Best Algo');
display(winT(1:10, :));
fminChart = winT(1:10, :);
    %FUNSET = [TYPE_C, TYPE_RUNGE];
    %     kFirst = zeros(length(FUNSET), 3, length(ALGOSET));
    %     kFirstGOF  = zeros(length(FUNSET), 2, length(ALGOSET));
    %     for fun=FUNSET(2)
    %         for algo=ALGOSET([2:3 6 8:length(ALGOSET)])
    %             k = findMinimum(C, T, X, [], lb, ub, A, b, fun, algo);
    %             Cfit = dCdT(k, T, C0, ssv, 1e-5);
    %
    %             if length(C) == length(Cfit)
    %                 [r2,rsme]=rsquare(C, Cfit);
    %                 kFirstGOF(fun, :, algo) = [r2,rsme];
    %             end
    %             kFirst(fun, :, algo) = k;
    %         end
    %     end
    %     testRes{testNum}.kFirst = kFirst;
    %     testRes{testNum}.kFirstGOF = kFirstGOF;
    
%     for i = 1:length(lamdaSet)
%         lamdaReg = lamdaSet(i);
%         index=(testNum-1)*length(lamdaSet)+i;
%         dataout(index, 1) = testNum;
%         dataout(index, 2) = lamdaReg;
%         k = findMinimum(C, T, X, [], lb, ub, A, b, FUNCTION_TYPE, ALGORITHM);
%         dataout(index, 3:5) = k;
%         Cfit = dCdT(k, T, C0, ssv, 1e-5);
%         if length(C) == length(Cfit)
%             [r2,rsme]=rsquare(C, Cfit);
%             dataout(index, 6:7) = [r2,rsme];
%         end
%     end
%tt = array2table(dataout, 'VariableNames',{'dataset', 'lamda','K_max','K_s','K_I', 'r2', 'rmse'});
%tt = sortrows(tt,{'dataset','r2'},{'ascend','descend'});
%display(tt);