function [b,a] = SingularValueAllpassSwitchCompact(Theta,N,sigma);
% [b,a] = SingularValueAllpassSwitchCompact(Theta,N,sigma);
% 
% Compact design of an allpass filter that can be used to construct a switching
% function at frequencies Omega. The overall order of the allpass filter is N,
% and the design can be weighted by the time domain samples of the singular 
% value in sigma.
%
% The design is a complex-valued modification of Mathias Lang's approach 
% published in "Matt's DSP Blog" in 2016.
% 
% Input parameters:
%    Theta        vector with even number of switching frequencies 0< Omega <2pi
%    N            order of allpass filter
%    sigma        time domain coefficients of singular value for weighting (row
%                       vector, optional)
%
% Output parameters:
%    b            numerator coefficients of complex allpass
%    a            denominator coefficients of complex allpass

% Sebastian J. Schlecht and Stephan Weiss, UoS, 20/4/2024


%-----------------------------------------------------------------------------
%  parameters
%-----------------------------------------------------------------------------
Nfft = 4*1024;
Omega = (0:Nfft)'/Nfft*2*pi;
if nargin<3,
   w = ones(Nfft+1,1);
else
   w = sqrt(abs(fft(sigma,Nfft)));
   w = w(:);
   w = [w; w(1)];
end;
   
%-----------------------------------------------------------------------------
%  target filter
%-----------------------------------------------------------------------------
Lth = length(Theta);
Phi = -(N-Lth/2)*Omega;
for lth = 1:Lth,
   Phi = Phi - pi*(Omega>Theta(lth));
end;
   
%-----------------------------------------------------------------------------
%  allpass implemention analogous to Mathias (`Matt') Lang, DSP Blog 2022
%-----------------------------------------------------------------------------
nn = (1:N);
theta = sum(Theta-pi)/2; 

Ar = exp(1i*theta) * exp(-1i*Omega*(N-nn)) - exp(1i*(Phi(:,ones(1,N)) - Omega*nn));
Ai = -1i*( exp(1i*theta) * exp(-1i*Omega*(N-nn)) + exp(1i*(Phi(:,ones(1,N)) - Omega*nn)) );
%  weighting analogoues to Matt Lang
b = diag(w)*(exp(1i*Phi) - exp(1i*(theta-Omega*N)));
A = diag(w)*[Ar Ai];
a0 = pinv([ real(A); imag(A) ])* [ real(b); imag(b) ];
% weighted least squares should work differently --- applied as a regularisation

%-----------------------------------------------------------------------------
%  output parameters
%-----------------------------------------------------------------------------
a = [1; a0(1:N)+1i*a0(N+1:2*N)];
b = exp(1i*theta)*flipud(conj(a));
