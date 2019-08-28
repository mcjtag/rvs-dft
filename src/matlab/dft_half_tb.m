% Computation of the N/2-point DFT of a N-point real sequence

clear all;

N=1024;
fs=1000;
dt=1/fs;
t=0:dt:(N-1)*dt;
f0=50;
f1=75;
f2=100;
f3=200;
f4=400;

% Generate one real-valued N-point sequence
x=sin(2*pi*t*f0)+cos(2*pi*t*f1)+cos(2*pi*t*f2)+cos(2*pi*t*f3)+sin(2*pi*t*f4);
X=fft(x); % N-point FFT
y=zeros(1, N/2);
% Generate complex-valued N/2-point sequence
for k=1:N/2
	y(k)=x(2*k-1)+1i*x(2*k);
end
Y = fft(y); % N/2-point FFT

% Recovery N-point sequence
y = dft_half(Y);

figure; hold on; plot(abs(X)); grid;
figure; hold on; plot(abs(Y)); grid;
figure; hold on; plot(abs(y)); grid;

% Standard deviation (error) of real and imaginary components 
errR = std(real(y)-real(X));
errI = std(imag(y)-imag(X));