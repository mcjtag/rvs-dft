% Computation of fast convolution with N/2-point DFT (type 2)
clear all;

N=1024;
fs=1000;
dt=1/fs;
t=0:dt:(N-1)*dt;
f0=30;
f1=50;
f2=200;

% Generate real-valued sequence
x=sin(2*pi*t*f0)+sin(2*pi*t*f1)+sin(2*pi*t*f2);

% Generate impulse response (LP filter)
fd = designfilt('lowpassiir','FilterOrder',8,'PassbandFrequency',100,'PassbandRipple',0.2,'SampleRate',1000);
h=impz(fd);
h(numel(h):N)=0;
h=h';

% Calculate reference sequence with the N-point DFT
YE=fft(x).*fft(h);
ye=ifft(YE);

% Generate complex-valued sequence
xc=zeros(1,N/2);
for k=1:N/2
	xc(k)=x(2*k-1)+1i*x(2*k);
end
X=fft(xc); % N/2-point DFT

% Generate complex-valued sequence of impulse response
hc=zeros(1,N/2);
for k=1:N/2
	hc(k)=h(2*k-1)+1i*h(2*k);
end
H=fft(hc);  % N/2-point DFT

% Fast convolution function
Y=dft_fconv2(X, H);
yc=ifft(Y); % N/2-point DFT

y=zeros(1,N);
% Recover result sequence from complex-valued to real-valued
for k=1:N/2
	y(2*k-1)=real(yc(k));
 	y(2*k)=imag(yc(k));
end

% Plot original and filtered signals 
figure; hold on; plot(abs(fft(x))); plot(abs(YE)); grid;
figure; hold on; plot(abs(fft(x))); plot(abs(fft(y))); grid;

% Standard deviation (error) of N/2 and N calculations 
err=std(y-ye);

