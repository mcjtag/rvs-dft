function [y] = dft_fconv1(x0, x1)
% Computation of fast convolution with N/2-point DFT (type 1)
% Input:
%	x0 - N/2-point complex-valued sequence 
%		x = DFT[ s(2k-1)+i*s(2k) ], where s is N-point real-valued sequence
%	x1 - N-point real-valued sequence
% Output:
%	y - N/2-point complex-valued sequence, y = x0*x1
N=numel(x0)*2;

h0=zeros(1,N/2);
h1=zeros(1,N/2);
w0=zeros(1,N/2);
w1=zeros(1,N/2);
w2=zeros(1,N/2);

y=zeros(1,N/2);

for k=1:N/2
	w0(k)=0.5*(1-sin(2*pi*(k-1)/N));
	w1(k)=0.5*(1+sin(2*pi*(k-1)/N));
	w2(k)=0.5i*cos(2*pi*(k-1)/N);
end

for k=1:N/2
	h0(k)=w0(k)*x1(k)+w1(k)*conj(x1(N/2-k+2));
	h1(k)=w2(k)*(x1(k)-conj(x1(N/2-k+2)));
end

x0(N/2+1)=x0(1);
for k=1:N/2
	y(k)=x0(k)*h0(k)+conj(x0(N/2-k+2))*h1(k);
end

end