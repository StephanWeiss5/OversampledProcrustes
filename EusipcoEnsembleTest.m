% EusipcoEnsembleTest.m
%
% Matlab script file to generate the ensemble simulation for [1].
%
% [1] S. Weiss, S.J. Schlecht, M. Moonen: "Best Least Squares Paraunitary
%     Approximation of Matrices of Analytic Functions," submitted to IEEE
%     Trans. Signal Process., Mar. 2025.

filename = 'EusipcoEnsembleResults.txt';

MVals = [32];
%MVals = [64];
SeedVals = (101:105);
for m = 1:length(MVals),
  disp(sprintf('Matrix dimension: %d',MVals(m)));
  for i = 1:length(SeedVals),
    disp(sprintf('Seedvalue: %d',SeedVals(i)));
    Results = zeros(2,8);
    Results(:,1) = [1; 2];
    Results(1,2:8) = ProcrustesLargeMatrixSingleTest(SeedVals(i),MVals(m));
    Results(2,2:8) = ProcrustesNewLargeMatrixSingleTest(SeedVals(i),MVals(m));
    
    %--------------------------------------------------------
    %  write results to a file
    %--------------------------------------------------------
    if exist(filename,'file') ~= 2,
      disp('new results file created');
      dlmwrite(filename,Results);
    else  
      disp('results appended');
      dlmwrite(filename,Results,'-append');
    end;
  end;  
end;
