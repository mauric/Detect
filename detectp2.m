%% TP Detection - Exo 2
%Detection d'un signal constant dans un bruit

clc
clear all
close all


A =.5;
sigmacarre = 2 % selon relation marqué en point 3
N = 100;
M = 1000;


%% HYPOTHESE H0
for i=1 : M
    x  = sqrt(sigmacarre)*randn(1,N);
    A = mean(x);
    sigma1 = 1/N*x*x';
    sigma2 = 1/N*(x-A)*(x-A)';
    Lg0(i) = N*log(sigma1/sigma2);
end

mean(Lg0);
var(Lg0);
% figure();
% hist(T,100);
% title('HYPOTHESE H0')


%% HYPOTHESE H1
for i=1 : M
    x  = sqrt(sigmacarre)*randn(1,N)+A;
    A = mean(x);
    sigma1 = 1/N*x*x';
    sigma2 = 1/N*(x-A)*(x-A)';
    Lg1(i) = N*log(sigma1/sigma2);
end


figure();
subplot(2,1,1)
hist(Lg1,100);
grid on
title('HYPOTHESE H1')
subplot(2,1,2)
hist(Lg0,100);
grid on
title('HYPOTHESE H0')

%% Calcul de PFA et PD

gama=[0.0:0.1:10];
PFA =zeros(1,length(gama));
for n = 1:length(gama)
    for i= 1:M
        if(Lg0(i)>gama(n))PFA(n) = PFA(n)+1;
        end
    end
    PFA(n) = PFA(n)/M; 
   
end

PD =zeros(1,length(gama) );
for n = 1:length(gama) 
    for i= 1:M
        if(Lg1(i)>gama(n))
            PD(n) = PD(n)+1;
        end
    end
    PD(n) = PD(n)/M;
end

figure()
hold on
plot(PFA,PD)
grid on
xlabel('PFA');
ylabel('PD');
title('ROC')

%% Test signal inconnu
load('signalinconnu.mat','signalinconnu');
gama2 = 4;
for j= 1:50
     sig = signalinconnu(1+(j-1)*100:j*100);
     A = mean(sig)
    sigma1 = 1/N*sig*sig';
    sigma2 = 1/N*(sig-A)*(sig-A)';
    Lsig(j) = N*log(sigma1/sigma2);
 
    if(Lsig(j)>gama2)
        H(j) = 1;
    else
        H(j) =0;
    end
end    
    
figure()
subplot(3,1,1);
plot(H,'-ob')
grid on
subplot(3,1,2);
plot(signalinconnu)
grid on
title('signal inconnu, pfa=0.1')
subplot(3,1,3);
hold on
plot(Lsig)
plot(gama2*ones(1,size(Lsig,2)))
title('stat signal inconnu, pfa=0.1')





















