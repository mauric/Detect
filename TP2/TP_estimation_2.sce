// TP Estimation

close 
clear


A = 1;
phi = 1;
fo = 0.25;
sigma = [0.25 0.5 1];
N = [10 50 100 500 800 10000];
for i = 1 : 3
    for j = 1 : 6
        t = (0:N(j)-1)';
        x = A*cos(2*%pi*fo*t+phi)+ (grand(N(j),'mn',0,sigma(i)^2))';
        
        A_ch= sum(x.*cos(2*%pi*fo*t+phi))/sum(cos(2*%pi*fo*t+phi).^2);
        phi_ch = atan(sum(-x.*sin(2*%pi*fo*t)),sum(x.*cos(2*%pi*fo*t)));
        
        A_ch_v(j) = A_ch;
        phi_ch_v(j) = phi_ch;
        
        errA(i,j) = abs(A_ch - A);
        errPhi(i,j) = abs(phi_ch - phi);
    end
end

mean(A_ch_v)
mean(phi_ch_v)


figure(1);
subplot(2,3,1);
title("estimation error of A for sigma = 0.25");
plot(N,errA(1,:),'-*');
subplot(2,3,4);
title("estimation error of phi for sigma = 0.25");
plot(N,errPhi(1,:),'-*');
subplot(2,3,2);
title("estimation error of A for sigma = 0.5");
plot(N,errA(2,:),'-*');
subplot(2,3,5);
title("estimation error of phi for sigma = 0.5");
plot(N,errPhi(2,:),'-*');
subplot(2,3,3);
title("estimation error of A for sigma = 1");
plot(N,errA(3,:),'-*');
subplot(2,3,6);
title("estimation error of phi for sigma = 1");
plot(N,errPhi(3,:),'-*');


figure(2);
subplot(2,3,1);
title("estimation error of A for sigma = 0.25");
plot(N/(2*sigma(1)^2),errA(1,:),'-*');
subplot(2,3,4);
title("estimation error of phi for sigma = 0.25");
plot(N/(2*sigma(1)^2),errPhi(1,:),'-*');
subplot(2,3,2);
title("estimation error of A for sigma = 0.5");
plot(N/(2*sigma(2)^2),errA(2,:),'-*');
subplot(2,3,5);
title("estimation error of phi for sigma = 0.5");
plot(N/(2*sigma(2)^2),errPhi(2,:),'-*');
subplot(2,3,3);
title("estimation error of A for sigma = 1");
plot(N/(2*sigma(3)^2),errA(3,:),'-*');
subplot(2,3,6);
title("estimation error of phi for sigma = 1");
plot(N/(2*sigma(3)^2),errPhi(3,:),'-*');
