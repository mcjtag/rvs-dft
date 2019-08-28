function [y] = dft_half(x)
% Computation of the N/2-point DFT of a N-point real sequence (Recovery)
% Input:
% 	x - complex N/2-point sequence: x = DFT[ s(2k-1) + i*s(2k) ], where s is a N-point real-valued sequence
% Output:
%	y - complex recovered N-point sequence: y = DFT[s]
%
N=numel(x)*2;
w0=zeros(1, N/2);
w1=zeros(1, N/2);
y=zeros(1, N);

for k=1:N/2
	w0(k)=(1-1i*exp(-1i*2*pi*(k-1)/N))/2;
	w1(k)=(1+1i*exp(-1i*2*pi*(k-1)/N))/2;
end

x(N/2+1)=x(1);

for k=1:N/2
	y(k)= x(k)*w0(k)+conj(x(N/2-k+2))*w1(k);
end

y(N/2+1)=real(x(1))-imag(x(1));
for i=2:N/2
	y(N-i+2)=conj(y(i));
end

end

