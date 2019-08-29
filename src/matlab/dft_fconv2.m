function [y] = dft_fconv2(x0, x1)
% Computation of fast convolution with N/2-point DFT (type 2)
% Input:
%	x0 - N/2-point complex-valued sequence 
%		x0 = DFT[ s(2k-1)+i*s(2k) ], where s is N-point real-valued sequence
%	x1 - N/2-point complex-valued sequence of impulse response after DFT
%		x1 = DFT[ m(2k-1)+i*m(2k) ], where m is N-point real-valued sequence
% Output:
%	y - N/2-point complex-valued sequence, y = x0*x1
N=numel(x0)*2;

w0=zeros(1,N/2);
w1=zeros(1,N/2);

h0=zeros(1,N/2);
h1=zeros(1,N/2);

y=zeros(1,N/2);

for k=1:N/2
	w0(k)=0.25*(3-exp(-4i*pi*(k-1)/N));
	w1(k)=0.25*(1+exp(-4i*pi*(k-1)/N));
end

x1(N/2+1)=x1(1);
for k=1:N/2
	h0(k)=w0(k)*x1(k) + w1(k)*conj(x1(N/2-k+2));
	h1(k)=w1(k)*(x1(k) - conj(x1(N/2-k+2)));
end

x0(N/2+1)=x0(1);
for k=1:N/2
	y(k)= x0(k)*h0(k) + conj(x0(N/2-k+2))*h1(k);
end

end

