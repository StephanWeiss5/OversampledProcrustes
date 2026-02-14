function [Q,D,Sfcorr2,OSflag,Nfft,e] = PUProcrustes(B,A,Nfftmax,Delta,NOrd);
%[Q,D,S,OSFlag] = PUProcrustes_v5(B,A,Nfft,Delta,N);
% 
%  Q= PUProcrustes(B) finds the closest paraunitary matrix Q to the polynomial
%  matrix B such that (Q-B) is minimised in the least squares sense. B must be 
%  a square, full rank matrix.
%
%  Q= PUProcrustes(B,A) finds the closest paraunitary matrix such that the norm
%  (AQ-B) is minimised in the least squares sense. With A optional, its default 
%  value is an identity matrix.
%
%  Q= PUProcrustes(B,A,Nfft) uses Nfft frequency bins to solve the Procrustes 
%  problem. If not supplied, Nfft is the next-larger power of two to four times
%  the support of A'*B.
%
%  Q= PUProcrustes(B,A,Nfft,Delta) assumes that the Procrustes problem us non-
%  causal, and that e.g. a paraunitary approximation of z^Delta A^P(z)B(z) is
%  sought. If Delta is small compared to Nfft, the algorithm works acceptably
%  even with Delta=0.
%
%  Input parameters:
%      B      NxMxL1 matrix
%      A      NxMxL2 matrix  (optional, default A=eye(M) )
%      Nfft   number of frequency bins
%      Delta  advance in case A^P B is non-causal
%      N      order of allpass filters
%
%  Output parameters:
%      Q      MxMxL3 paraunitary matrix
%      D      sign changes applied to bins
%      S      extracted singular values
%      OSflag (1 if OK, 0 if oversampling by 2 is required)

%  S. Weiss and S. Schlecht, 2025

% try to find handle on sign change in case singular values go negative - 5/9/23
% sign change / permutation detection according to lastest draft 'c'
%      and refined search for zero-crossings -- 3/4/24
% trying to include iteration for DFT size --- 4/4/24
% version 5 --- trying to sort some issue with bin 0 --- 9/4/24
%           --- also changed allpass design and iteration with OSflag, which can be 
%                       incorrect for lower DFT lengths --- 20/4/24
% version 6 --- order of allpass filters as input

QuietMode=true;

%------------------------------------------------------------------------------
% parameters and initialisations 
%------------------------------------------------------------------------------
[N,M,~] = size(B);
if nargin==1,
   if N~=M, 
      error('incorrect input dimensions');
   end;
   R = B;
else
   [M2,N2,~]=size(A);
   if M2~=M, error('mismatching dimensions of input matrices'); end;
   R = PolyMatConv(ParaHerm(A),B);
end;      
if nargin<3,
   Nfftmax = 2048; 
end;
if nargin<4,
   Delta = 0;
end;   
EpsSV = 1e-10;                             % threshold for singular values to be considered zero or close to zero
AllpassOrd = NOrd;                         % order for allpass filter design
Nfft = 2^(ceil(log2(size(R,3))+1));        % initial DFT size
Nfft = 512;
      
%------------------------------------------------------------------------------
% initial bin-wise SVD 
%------------------------------------------------------------------------------
Rf = zeros(M,M,Nfft);                        % ultimately DFT domain of A^P(z)B(z)
Rf(:,:,1:size(R,3)) = R;                     % define causal transfer function matrix
Rf = fft(circshift(Rf,-Delta,3),Nfft,3);     % overwrite with its Fourier coefficients
Uf = zeros(M,M,Nfft);
Vf = zeros(M,M,Nfft);
Sf = zeros(M,Nfft);
for k = 1:Nfft,
   [Uf(:,:,k),s,Vf(:,:,k)] = svd(Rf(:,:,k));
   Sf(:,k) = diag(s);
end;
crit = 0; IterNum = 0; Q = zeros(M,M,1); Q(:,:,1) = eye(M);

%------------------------------------------------------------------------------
%
%  iterate for DFT size 
%
%------------------------------------------------------------------------------
while crit == 0,
   Nfft = 2*Nfft;                            % double DFT size
   IterNum = IterNum+1;
   if QuietMode==false, disp(sprintf(' iteration %d with DFT length %d',[IterNum,Nfft])); end;

   %---------------------------------------------------------------------------
   % bin-wise SVD in added bins 
   %---------------------------------------------------------------------------
   Rf = zeros(M,M,Nfft);                     % ultimately DFT domain of A^P(z)B(z)
   Rf(:,:,1:size(R,3)) = R;                  % define causal transfer function matrix
   Rf = fft(circshift(Rf,-Delta,3),Nfft,3);  % overwrite with its Fourier coefficients
   Uff = Uf; Vff = Vf; Sff = Sf;
   Uf = zeros(M,M,Nfft); Uf(:,:,1:2:Nfft) = Uff;
   Vf = zeros(M,M,Nfft); Vf(:,:,1:2:Nfft) = Vff;
   Sf = zeros(M,Nfft); Sf(:,1:2:Nfft) = Sff;
   for k = 2:2:Nfft,
      [Uf(:,:,k),s,Vf(:,:,k)] = svd(Rf(:,:,k));
      Sf(:,k) = diag(s);
   end;
   
   
   %---------------------------------------------------------------------------
   %  determine zero singular values and their corresponding subspaces 
   %---------------------------------------------------------------------------
   Lk = sum(Sf>EpsSV,1);        % # of non-zero singular values per bin
   L = max(Lk);                 % # estimated number of non-vanishing singular values
   % we need to address zeros in non-vanishing singular values, either by removing 
   %    that bin (see Sebastian's approach) or by moving the bin (Marc Moonen/Faizan)
   ZeroBinIndex = find((L-Lk)>0.1);     % indices of bins with at least one spectral zero;
   % we assume that there is a sufficient number of bins without zero crossings (>Nftt/2)
   if length(ZeroBinIndex)>=Nfft/2,
      %  we would now increase the DFT length; for the moment, let's terminate
      error('too many spectral zeros in non-vanishing singular values');
   end;   
   NonZeroBinIndex = find((L-Lk)<=0.1);   % complement to vanishing singular values
   kbins = 1:Nfft;                        % bin locations
   NonZeroBins = kbins(NonZeroBinIndex);
   % some diagnostics:
   if QuietMode==false,
      disp(sprintf('   number of vanishing singular values: %d',M-L));
      disp(sprintf('   bins with spectral zeros in non-vanishing singular values: %d',length(ZeroBinIndex))); 
   end;   
   Ufcorr = Uf(:,1:L,NonZeroBins);
   Vfcorr = Vf(:,1:L,NonZeroBins);  
   Sfcorr = Sf(1:L,NonZeroBins);
         
   %---------------------------------------------------------------------------
   % checking for algebraic multiplicities in non-vanishing singular values
   %---------------------------------------------------------------------------
   % eliminate bins with non-trivial algebraic multiplicities
   %    otherwise these create unwanted ambiguities when looking for sign changes
   DiffEps = 1e-10;
   NoAlgMultIndices=find(sum(abs(diff(Sfcorr,1,1))<DiffEps,1)<0.1);
   NonZeroBins2 = NonZeroBins(NoAlgMultIndices);
   if QuietMode==false,
      disp(sprintf('   bins with non-trivial algebraic multiplicities of singular values: %d',...
           length(NonZeroBins)-length(NonZeroBins2))); 
   end;        
   Ufcorr2 = Ufcorr(:,:,NoAlgMultIndices);
   Vfcorr2 = Vfcorr(:,:,NoAlgMultIndices);
   Sfcorr2 = Sfcorr(:,NoAlgMultIndices);

   %---------------------------------------------------------------------------
   % find permutations, and check for sign changes in singular values
   %---------------------------------------------------------------------------
   D = ones(L,length(NonZeroBins2)); 
   for k = 2:length(NonZeroBins2),
      % check for best alignment to previous bin
      GU = Ufcorr2(:,:,k-1)'*Ufcorr2(:,:,k);   
      GV = Vfcorr2(:,:,k-1)'*Vfcorr2(:,:,k); 
      % determine permutation of singular values and singular vectors from bin (k-1) to bin k
      P = BestPermutation(GU,GV);  
      % in order to track sign changes across bins, we need to get the association of singular values
      %    and vectors right
%      if k ==2, save; end;
      D(:,k) = sign(real(diag(diag(GU*GV'))))*D(:,k-1);
      Ufcorr2(:,:,k) = Ufcorr2(:,:,k)*P';          % contains samples of analytic left- and
      Vfcorr2(:,:,k) = Vfcorr2(:,:,k)*P';          %     right-singular vectors
      Sfcorr2(:,k) = D(:,k).*(P*Sfcorr2(:,k));     % contains now samples of the analytic singular values
      % the sign change D is not applied to any of the singular vectors ---- instead, it is performed 
      % below via an allpass switching function
   end;
   if QuietMode==false,
      figure(100); clf;
      plot((NonZeroBins2-1)/Nfft,real(Sfcorr2)','-*'); hold on; 
      xlabel('norm. freq.'); ylabel('approx. singular values');
   end;    
  
   %---------------------------------------------------------------------------
   %  check for 4-pi periodicity
   %---------------------------------------------------------------------------
   GU = Ufcorr2(:,:,end)'*Ufcorr2(:,:,1);      % checking for a sign change between last and first bin
   GV = Vfcorr2(:,:,end)'*Vfcorr2(:,:,1); 
   P = BestPermutation(GU,GV);  
   D0 = sign(real(diag(diag(GU*GV'))))*D(:,end);
   % there should be no sign mismatch between the interval margins
   if norm(D0-D(:,1))>0.1,
      disp('the matrix possesses singular values with a periodicity greater than 2 pi');
      OSflag=0;
   else
      OSflag=1;   
   end;   

   %---------------------------------------------------------------------------
   %  sign switch using allpass functions
   %---------------------------------------------------------------------------
   if OSflag==1,
      % singular values should all be positive in the first bin
      % determine zero crossings
      Delays = zeros(size(Sfcorr2,1),1);
      At = zeros(size(Sfcorr2,1),Nfft);
      Af = zeros(size(Sfcorr2,1),Nfft);
      % determine time domain coefficients of singular values in order to
      %   refine the search for any zero crossings
      if length(NonZeroBins2)<Nfft,
         T = dftmtx(Nfft);
         Tinv = pinv(T(NonZeroBins2,[3*Nfft/4+1:Nfft, 1:Nfft/4+1]));
         s = (Tinv*Sfcorr2.').';
      else
         dummy = ifft(Sfcorr2,Nfft,2);
         s = dummy(:,[3*Nfft/4+1:Nfft, 1:Nfft/4+1]);
      end;      
      for m = 1:size(Sfcorr2,1),                % iterate for singular values
         % find sign changes
         ZerosIndex = find(abs(diff(D(m,:))));
         if mod(length(ZerosIndex),2)~=0, 
            disp('number of zero crossings should be even; oversampling flag failed');
         end;  
         if length(ZerosIndex)>1,               % singular value has zero crossings
            if QuietMode==false,
                disp(sprintf('   singular value has %d zero crossings ... designing allpass',length(ZerosIndex)));
            end;    
            Omegas = zeros(length(ZerosIndex),1); 
            for i = 1:length(ZerosIndex),
               % time domain singular values for zero-crossing interpolation
               if ZerosIndex(i)==1,                       % catch the case where the zero
                                                          %    crossing is detected at the 
                                                          %    first bin            
                  Omegas(i) = ZeroCrossing(s(m,:),2*pi/Nfft*(NonZeroBins2(ZerosIndex(i))-1),...
                      2*pi/Nfft*(NonZeroBins2(ZerosIndex(i))) );
               else
                  Omegas(i) = ZeroCrossing(s(m,:),2*pi/Nfft*(NonZeroBins2(ZerosIndex(i)-1)),...
                      2*pi/Nfft*(NonZeroBins2(ZerosIndex(i))) );
               end;       
            end;
%            if QuietMode==false, plot(Omegas/2/pi,zeros(size(Omegas)),'ko'); drawnow; end;
            % previous solution based in modulated sections of real allpasses
            % [b,a] = SingularValueAllpassSwitch(Omegas,AllpassOrd);  
            % Delays(m) = (AllpassOrd-1)*length(ZerosIndex)/2;            
            [b,a] = SingularValueAllpassSwitchCompact(Omegas,AllpassOrd);  
            Delays(m) = AllpassOrd-length(Omegas)/2;            
            At(m,:) = impz(b,a,Nfft).';         % non-causal allpass filter coefficients
         else                                   % singular values has NO zero crossings
            if QuietMode==false, disp('   no zero crossings found for singular value'); end;
            At(m,1) = 1; 
            Delays(m) = 0;
         end; 
         Af(m,:) = fft(circshift(At(m,:),-Delays(m),2),Nfft,2);  
      end;            
      MaxDelay = max(Delays);   
      D = D.*Af(:,NonZeroBins2);

      %---------------------------------------------------------------------------
      %  Procrustes solution in reduced number of bins
      %---------------------------------------------------------------------------
      Qfcorr2 = zeros(M,M,length(NonZeroBins2));
      for k = 1:length(NonZeroBins2),
         Qfcorr2(:,:,k) = Ufcorr2(:,:,k)*diag(D(:,k))*Vfcorr2(:,:,k)';
      end;
   
      %---------------------------------------------------------------------------
      %  time domain Procrustes reconstruction ... Sebastian's trick
      %---------------------------------------------------------------------------
      if length(NonZeroBins)<Nfft,
         Qtcorr2 = PolyMatPartialIDFT(Qfcorr2,Nfft,NonZeroBins2);
      else
         Qtcorr2 = ifft(Qfcorr2,Nfft,3); 
      end;   
      LQ = size(Qtcorr2,3);  
      Qnew = PUPolyMatTrim(Qtcorr2(:,:,[LQ/2+1:LQ,1:LQ/2]),0);        % wrap around to have time zero centred
      %---------------------------------------------------------------------------
      %  check difference to previous iteration
      %---------------------------------------------------------------------------
      [~,e,~,~] = PolyMatAlign(Qnew,Q);
      Q = Qnew;
      if QuietMode==false,
         disp(sprintf('   iteration %d complete; mismatch to previous iteration: %2.12f',[IterNum e]));
      end;   
      if (e < 1e-8) || (Nfft>=Nfftmax)
         crit = 1;
      end;
   else      % OSflag == 0
      if Nfft>=Nfftmax,
         crit = 1;
      end;   
   end;             
end;

   
%------------------------------------------------------------------------------
% function BestPermutation()
%------------------------------------------------------------------------------
function P = BestPermutation(X,Y);
%  X is close to a permutation matrix; BestPermutation returns the permutation
%  matrix closest to X and Y in terms of searching the maximum absolute element
%  per row.
%  The search is performed on X. Y is used for control --- a mismatch returns an
%  error.
M = size(X,1);
P = zeros(M,M); 
for m = 1:M,
   [~,i] = max(abs(X(m,:)));
   [~,j] = max(abs(Y(m,:)));
   P(m,i) = 1; 
   if i~=j,
      disp('mismatch in permutation matrices');
      return;
   end;   
end;   
   
%------------------------------------------------------------------------------
%  function ZeroCrossing() 
%------------------------------------------------------------------------------
function Omega = ZeroCrossing(s,Omega1,Omega2);
% Refine the localisation of a zero crossing between the frequencies Omega1
% and Omega2, for a function with a real-values spectrum given by its time-
% domain coefficients s.

Ls = length(s); 
nn = (-(Ls-1)/2:(Ls-1)/2)';
S1 = real(s*exp(-1i*Omega1*nn));
S2 = real(s*exp(-1i*Omega2*nn));
% bisection method for refinement of zero location
if S2<S1,             % function is falling; convert to rising
   s = -s; 
end;   
for i = 1:10,
   Omega = (Omega1+Omega2)/2;
   S12 = real(s*exp(-1i*Omega*nn));
   if S12<0;
      Omega1 = Omega;
   else 
      Omega2 = Omega;   
   end;
end     

