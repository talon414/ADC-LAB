clc;
clear all;
close all;

N=1e3;
M=[2,4];
B=[0.3,0.5,0.8];
ii=1;
lvls=-(M(ii)-1):2:(M(ii)-1);
d=lvls(randi(length(lvls),1,N));
L=8;
v=[];
for i=1:N
    v=[v,d(i),zeros(1,L-1)];
end

Ts=1;
Nsym=8;
t=-Nsym/2:1/L:Nsym/2;
pt=P_sig(t,B(1),Ts);

s=conv(v,pt);

EbN0dB=15;
EsN0dB=10*log10(log2(M))+EbN0dB;

D=L*Nsym/2;

P_s=Pwr(L,s);

noise_psd=0.5*P_s/power(10,EsN0dB(2)/10);

r=s+power(noise_psd,1/2).*randn(1,length(s));

scap=conv(r,pt);
%scap=scap((2*D+1):L:(length(scap)-2*D));
%Ep=sum(pt.*pt);
%scap=scap./Ep;

plotEyeDiagram(scap,L,2*L,2*D,100);

function idag = plotEyeDiagram(scap,L,m,offset,t)
    N1=m*t;
    idag=reshape(scap(offset+1:offset+N1),m,t);
    tim=(0:1:m-1)/L;

    plot(tim,idag);

end


function pwr = Pwr(L,s)
    pwr=L*sum(s.*s)/length(s);
end

function pvals = P_sig(t,b,Ts)
    %pvals(t) = (sin(pi.*t*(1-beta)/Ts)+4*beta.*t*cos(pi*t*(1+beta)./Ts))./(Ts*(pi.*t*(1-power(4*beta*t/Ts,2))/Ts));
    for i=1:length(t)
        pvals(i)=(sin(pi*t(i)*(1-b)/Ts)+4*b*t(i)*cos(pi*t(i)*(1+b)/Ts))/(Ts*(pi*t(i)*(1-power(4*b*t(i)/Ts,2))/Ts));
    end
    pvals(t==0)=(1+b*(4/pi-1))/Ts;
    pvals(t==abs(Ts/(4*b)))=b*((1+2/pi)*sin(pi/(4*b))+(1-2/pi)*cos(pi/(4*b)))/(power(2,1/2)*Ts);
end
