%% TP Detection
%Detection d'un signal constant dans un bruit

clc
clear all
close all


A =.5;
N = 1000;
M = 1000;


%% HYPOTHESE H0
for i=1 : M
    x  = randn(1,N);
    T(i) = mean(x);
end

mean(T);
var(T);
% figure();
% hist(T,100);
% title('HYPOTHESE H0')


%% HYPOTHESE H1
for i=1 : M
    x  = randn(1,N)+A;
    U(i) = mean(x);
end
% figure();
% hist(U,100);
% title('HYPOTHESE H1')

figure();
subplot(2,1,1)
hist(U,100);
title('HYPOTHESE H1')
subplot(2,1,2)
hist(T,100);
title('HYPOTHESE H0')
%% Calcul de PFA et PD

gama=[0.0:0.01:1.0]
PFA =zeros(1,length(gama));
for n = 1:length(gama)
    for i= 1:M
        if(T(i)>gama(n))PFA(n) = PFA(n)+1;
        end
    end
    PFA(n) = PFA(n)/M; 
   
end

PD =zeros(1,length(gama) );
for n = 1:length(gama) 
    for i= 1:M
        if(U(i)>gama(n))
            PD(n) = PD(n)+1;
        end
    end
    PD(n) = PD(n)/M;
end

figure()
axis([-.5,1,0,1.5])
hold on
plot(PFA,PD)
grid on
xlabel('PFA');
ylabel('PD');
title('ROC')


%% Test signal connu

load('signalconnu.mat','signalconnu');
gama2 = 0.05
for j= 1:50
    Tsig(j) = mean(signalconnu(1+(j-1)*100:j*100));
    if(Tsig(j)>gama2)
        H(j) = 1;
    else
        H(j) =0;
    end
end    
    
figure()
subplot(3,1,1);
plot(H,'-ob')
subplot(3,1,2);
plot(signalconnu)
title('signal connu, pfa=0.1')
subplot(3,1,3);
plot(Tsig)
title('stat signal connu, pfa=0.1')


h = get(0,'children');
for i=length(h):-1:1
  saveas(h(i), ['detectionExo1_' num2str(length(h)+1-i)], 'png');
end



















