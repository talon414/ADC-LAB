%Emmanuel
clc;
clear all;
close all;

f=2;
fs=16;
t=0:1/fs:1;
L=[256,128,64,32,16];
R=log2(L);
Snr=[];
X=[];
x=255*(1+sin(2*pi*f*t))/2;

for lvl=1:length(L)
    figure();
    X=Quant(x,L(lvl));
    subplot(3,1,1)
    stem(t,de2bi(X))
    xlabel('time');
    subplot(3,1,2)
    plot(t,x,'Color','green','LineWidth',1);
    xlabel('time');
    legend('original signal');
    grid on
    subplot(3,1,3)
    plot(t,X,'Color','red','LineWidth',0.4);
    xlabel('time');
    grid on
    title(L(lvl));
    legend('reconstructed signal');
    Snr=[Snr,SNRdB(X,x,fs)];
end
figure();
plot(flip(R),flip(Snr));
xlabel('bits');
ylabel('SNR(dB)');
grid on

function X = Quant(x,lvl)
    if lvl==256
        X=round(x);
    elseif lvl==128
        X=NEQ(x,2);
    elseif lvl==64
        X=NEQ(x,4);
    elseif lvl==32
        X=NEQ(x,8);
    elseif lvl==16
        X=NEQ(x,16);
    end
end


function neq = NEQ(x,mul)
    for i=1:length(x)
        if rem(x(i),mul)~=0
            x(i)=round(x(i)/mul)*mul;
        end
    end
    neq=x;
end

function snr = SNRdB(X,x,fs)
    differ=(X-x).^2;
    Px=power(255/2,2)/2;
    Pn=sum(differ)/fs;
    snr=10*log10(Px/Pn);
end
