function [xr, xi] = dft_split(s)
% Computation of the DFT of two real sequences (Splitting)
% Input:
%	s = DTF(x0 + i*x1)
% Output: 
%	xr = DFT(x0)
%	xi = DFT(x1)

N = numel(s);
xr = zeros(1, N);
xi = zeros(1, N);

xr(1) = real(s(1)) + 1i*0;
xi(1) = imag(s(1)) + 1i*0;
xr(N/2+1) = real(s(N/2+1));
xi(N/2+1) = imag(s(N/2+1));

for k=2:N/2
	xr(k) = (s(k)+conj(s(N-k+2)))/2;
	xi(k) = (s(k)-conj(s(N-k+2)))/2i;
end

for k=2:N/2
	xr(N/2+k) = conj(xr(N/2-k+2));
	xi(N/2+k) = conj(xi(N/2-k+2));
end

end

