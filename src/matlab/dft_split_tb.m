% Computation of the DFT of two real sequences
clear all;

N=1024;
fs=1000;
dt=1/fs;
t=0:dt:(N-1)*dt;
f0=20;
f1=40;

% Generate two real signals
s0=cos(2*pi*t*f0);
s1=sin(2*pi*t*f1);
sc=s0+1i*s1;

% Compute the DFT
s0f=fft(s0);
s1f=fft(s1);
scf=fft(sc);

% Split a complex-valued sequence
[xr,xi]=dft_split(scf);

% Plotting the complex-valued and two splitted real-valued sequences
figure; hold on; plot(abs(scf)); plot(abs(xr)); plot(abs(xi)); grid;

% Standard deviation (error) of real and imaginary components
errRr = std(real(xr)-real(s0f));
errRi = std(imag(xr)-imag(s0f));
errIr = std(real(xi)-real(s1f));
errIi = std(imag(xi)-imag(s1f));
