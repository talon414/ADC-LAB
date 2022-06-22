clc;
clear all;
close all;

N=1e5;
M=4;
L=8;
lvls=[-3,-1,1,3];
d=lvls(randi(length(lvls),1,N));
v=[];
for i=1:N
    v=[v d(i) zeros(1,L-1)];
end


b=0.3;
Ts=1;
t=-4:1/8:4;
p_t=P_sig(t,b,Ts);
Nsym=8;
figure();
plot(t,p_t);
grid on

s=conv(v,p_t);

EbN0dB=[-2,12];
EsN0dB=10*log10(log2(M))+EbN0dB;

D=L*Nsym/2;

P_s=Pwr(L,s);

noise_psd=0.5*P_s/power(10,EsN0dB(2)/10);

r=s+power(noise_psd,1/2).*randn(1,length(s));

scap=conv(r,p_t);
vcap=scap((2*D+1):L:(length(scap)-2*D));


dcap=Demod(vcap,lvls);
figure();
subplot(211)
plot(s(1:500))
ylabel('s(t)')
grid on
subplot(212)
plot(scap(1:500))
ylabel('scap(t)')
grid on
figure();
subplot(211)
plot(r(1:100))
ylabel('r(t)')
grid on
subplot(212)
stem(vcap(1:20))
ylabel('vcap(t)')
grid on
figure();
subplot(211)
stem(d(1:20))
ylabel('d(t)')
grid on
subplot(212)
stem(dcap(1:20))
ylabel('dcap(t)')
grid on

function dem = Demod(sig,lvl)
    for i=1:length(sig)
        den=[];
        for j=1:length(lvl)
            den(j)=abs(1000*(lvl(j)-sig(i)));
        end
        dem(i)=lvl(den==min(den));
    end
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
