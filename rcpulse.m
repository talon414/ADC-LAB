%Emmanuel
clc;
clear all;
close all;

T=1;
L=100;
alpha=[0.4,0.6,0.8];
t=-4*T:1/L:4*T;
N=length(t)/(1/L);

for i=1:length(alpha)
    figure();
    [time,rv]=RCfunc(t,T,alpha(i),L,N);
    plot(t,rv);
    grid on
end


function [tvals,rvals] = RCfunc(t,T,alp,L,N)
    rvals=(double(sin(t*pi/T)).*cos(pi*alp*t/T))./((t*pi/T).*(1-(2*alp*t/T)));
    tvals=-N/2:1/L:N/2;
    rvals(t==0)=1;
    rvals(t==abs(T/(2*alp)))=alp*sin(pi/(2*alp))/2;
end

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
