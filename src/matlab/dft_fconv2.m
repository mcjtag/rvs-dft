function [y] = dft_fconv2(x, h)
% Computation of fast convolution with N/2-point DFT (type 2)
% Input:
%	x - N/2-point complex-valued sequence 
%		x = DFT[ s(2k-1)+i*s(2k) ], where s is N-point real-valued sequence
%	h - N/2-point complex-valued sequence of impulse response after DFT
%		h = DFT[ m(2k-1)+i*m(2k) ], where m is N-point real-valued impulse
%		response
% Output:
%	y - N/2-point complex-valued sequence, y = x*h
N=numel(x)*2;

w0=zeros(1,N/2);
w1=zeros(1,N/2);

h0=zeros(1,N/2);
h1=zeros(1,N/2);

y=zeros(1,N/2);

for k=1:N/2
	w0(k)=0.25*(3-exp(-4i*pi*(k-1)/N));
	w1(k)=0.25*(1+exp(-4i*pi*(k-1)/N));
end

h(N/2+1)=h(1);
for k=1:N/2
	h0(k)=w0(k)*h(k) + w1(k)*conj(h(N/2-k+2));
	h1(k)=w1(k)*(h(k) - conj(h(N/2-k+2)));
end

x(N/2+1)=x(1);
for k=1:N/2
	y(k)= x(k)*h0(k) + conj(x(N/2-k+2))*h1(k);
end

end

