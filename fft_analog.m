%Emmanuel
clc;
clear all;
close all;


T=[1,2];
Ts1=0.01*T(1);
Ts2=0.01*T(2);
t1=-4*T(1):Ts1:4*T(1);
t2=-4*T(2):Ts2:4*T(2);

x1=sig(t1,T(1));
x2=sig(t2,T(2));

figure();
[f1,fft1]=FFT_Analog(x1,1/Ts1,'single-sided');
plot(f1,fft1);
xlabel('frequency');
ylabel('Normalised FFT');
grid on


figure();
[f2,fft2]=FFT_Analog(x2,1/Ts2,'double-sided');
plot(f2,fft2);
xlabel('frequency');
ylabel('Normalised FFT');
grid on


function [freq,fvals] = FFT_Analog(x,Fs,type)
    N=2^nextpow2(length(x));
    X=abs(fft(x,N));
    X=X./sum(X);
    if type == "single-sided"
        fvals=X(1:N/2);
        freq=Fs*(1:N/2)/(N/2);
    elseif type == "double-sided"
        fvals=X;
        freq=Fs*(-N/2:N/2-1)/(N/2);
    end
end


function xvals = sig(t,T)
    xvals=double(sin(t*pi/T))./(t*pi/T);
    xvals(t==0)=1;
end
