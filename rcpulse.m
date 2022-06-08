%Emmanuel
clc;
clear all;
close all;

T=1;
L=100;
alpha=[0,0.4,0.6,0.8,1];
t=-4*T:1/L:4*T;
N=2^nextpow2(length(t));
RV=[];

for i=1:length(alpha)
    figure();
    [time,k1]=RCfunc(t,T,alpha(i),L,N);
    plot(t,k1);
    title('roll of',alpha(i));
    grid on
    RV(i,:)=k1;
    figure();
    subplot(211);
    [f2,fftr2]=FFT_Analog(RCspec(time,T,alpha(i)),1/L,'double-sided');
    plot(f2,fftr2)
    grid on
    subplot(212)
    [f1,fftr]=FFT_Analog(RV(i,:),1/L,'double-sided');
    plot(f1,fftr)
    grid on
end

function specvals = RCspec(f,T,alp)
    specvals=[];
    specvals(f<=(1-alp)/(2*T))=T;
    specvals((f>=(1-alp)/(2*T)) & (f<=(1-alp)/(2*T)))=T/2*(1+cos(pi*T/alp*(f((f>=(1-alp)/(2*T)) & (f<=(1-alp)/(2*T)))-(1-alp)/(2*T))));
    specvals(f>(1+alp)/(2*T))=0;
end

function [tvals,rvals] = RCfunc(t,T,alp,L,N)
    rvals=(double(sin(t*pi/T)).*cos(pi*alp*t/T))./((t*pi/T).*(1-(2*alp*t/T)));
    tvals=-N/2:N/2-1;
    rvals(t==0)=1;
    rvals(t==abs(T/(2*alp)))=alp*sin(pi/(2*alp))/2;
end

function [freq,fvals] = FFT_Analog(x,Fs,type)
    N=2^nextpow2(length(x));
    X=abs(fftshift(fft(x,N)));
    X=X./sum(X);
    if type == "single-sided"
        fvals=X(1:N/2);
        freq=Fs*(0:N/2-1)/(N);
    elseif type == "double-sided"
        fvals=X;
        freq=Fs*(-N/2:(N/2)-1)/(N/2);
    end
end
