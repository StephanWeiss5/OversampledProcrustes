function Results = ProcrustesLargeMatrixSingleTest(SeedVal,M);
% ProcrustesLargeMatrixSingleTest.m
% 
% Perform one simulation run within an ensemble of tests using randomised
% matrices with ground truth Procrustes solution. The randomisation is 
% initialised by the seed value SeedVal. The function returns a number of
% metrics for the simulations in [1] in the variable Results:
%    Results(1)       seed value
%    Results(2)       M
%    Results(3)       DFT length
%    Results(4)       length of paraunitary matrix
%    Results(5)       paraunitarity error
%    Results(6)       LS mismatch
%    Results(7)       computation time 
%
% Input parameter:
%       SeedVal            seed value for random number generator
%       M                  spatial dimension of matrix
%
% Output parameter:
%       Results              vector containing various metrics
%
% [1] S. Weiss, S.J. Schlecht, M. Moonen: "Best Least Squares Paraunitary 
%     Approximation of Matrices of Analytic Functions", submitted to IEEE
%     Trans. Signal Process., March 2025.

%-------------------------------------------------------------------
%   parameters
%-------------------------------------------------------------------
N = 64;                        %   order of allpass design
Nfft = 1024;
L = 1;
DispMode = 'off';

Results = zeros(1,7);
Results(1)= SeedVal;
Results(2)= M;

%-------------------------------------------------------------------------------
%  matrix A(z)
%-------------------------------------------------------------------------------
rand('seed',SeedVal);
Sigma = zeros(M,M,3);
for m=1:M/2,                            % 1/2 of the SVs with 2 zero crossings
  dummy = randn(1,1,3)+1i*rand(1,1,3);
  dummy = dummy + conj(dummy(1,1,3:-1:1));
  Sigma(m,m,:) = dummy/norm(squeeze(dummy));
end;
for m = M/2+1:M,                        % 1/2 of the SVs with 1 zero crossing
  theta = rand(1,1);
  if (abs(theta)-.5)<0.05, theta = theta+.25; end;
  Sigma(m,m,1:2) = [exp(-1i*theta*pi) exp(1i*theta*pi)]/sqrt(2);
end;   
U = PUPolyMatRand(M,1,SeedVal,'complex');                
V = PUPolyMatRand(M,1,SeedVal-1,'complex');                
A = PolyMatConv(U,PolyMatConv(Sigma,ParaHerm(V)));   
%A = Sigma;

%------------------------------------------------------------------------------
%  Procrustes solution --- using the TSP method
%------------------------------------------------------------------------------
A2 = zeros(M,M,2*size(A,3)-1);
A2(:,:,1:2:end) = A;
B = zeros(M,M,1); 
B(:,:,1) = eye(M);
tstart = tic;
[Q,~,~,~,Results(3),~] = PUProcrustesNew(A,B,Nfft,0,16);
Results(7) = toc(tstart);
%Q = Q(:,:,1:2:end);
QQ = PolyMatConvPH2Fast(Q,Q);
Lq = size(QQ,3);
QQ(:,:,(Lq+1)/2) = QQ(:,:,(Lq+1)/2) - eye(size(Q,1));
Results(5) = PolyMatNorm(QQ);

%------------------------------------------------------------------------------
%  truncate outer zeros of paraunitary matrix
%------------------------------------------------------------------------------
MaxElement = zeros(size(Q,3),1);
% element size
for i = 1:size(Q,3),
   MaxElement(i) = max(max(abs(Q(:,:,i))));
end;   
% find any leading small components
StartIndex = 1; EndIndex = size(Q,3);
while MaxElement(StartIndex)<1e-8, StartIndex = StartIndex+1; end;
while MaxElement(EndIndex)<1e-8, EndIndex = EndIndex-1; end;
Q = Q(:,:,StartIndex:EndIndex);  
Results(4) = size(Q,3);

[~,Results(6),~,~] = PolyMatAlign(A,Q);

%------------------------------------------------------------------------------
%  assign outputs
%------------------------------------------------------------------------------
%Results(1) = SeedVal;
%Results(2) = M;
%Results(3) = Nfft;
%Results(4) = size(Q,3);
%Results(5) = paraunitarity;
%Results(6) = error;
%Results(7) = CompTime;

