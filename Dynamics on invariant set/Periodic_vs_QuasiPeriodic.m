%% numerical simulation of a simple example system for (7.8)
% We cosider the system
%   \partial_t A = \eps^2*(A - (1+5i) A|A|^2 + (1+5i) A|B_o|^2)
%                           + \eps^4 B_o |A|^4 e^{i t}
%   \partial_t B_o = \eps^2*(B_o - B_o|B_o|^2 - (2-i) B_o|A^2|) 
%                           + \eps^4*A|A|^5 e^{-i t}

%% periodic

eps = 0.1;
Tend = 100000;
omegaA = 5;

vecField = @(t,X) [eps^2*X(1).*(1 - (1 + 5*1i)*abs(X(1)).^2 + (1+5*1i)*abs(X(2)).^2) + eps^4*X(2).*abs(X(1)).^4.*exp(1i * t);...
                    eps^2*X(2).*(1 - abs(X(2)).^2 - (2-1i)*abs(X(1)).^2) + eps^4*X(1).*abs(X(1)).^4.*exp(-1i * t)];

Xinit = [1,0];
[t,x] = ode45(vecField,[0,Tend],Xinit);

figure;
plot(angle(exp(1i*t)),angle(x(:,1).*exp(-1i*omegaA*eps^2*t)),'.');
xlabel("$\theta$","FontSize",28,"Interpreter","latex");
ylabel("$\phi_A$","FontSize",28,"Interpreter","latex");
title("$\varepsilon = 0.1$","FontSize",28,"Interpreter","latex");
ylim([-pi,pi]);
xlim([-pi,pi]);
set(gca,'FontSize',20)
saveas(gcf,'periodic.png')

figure;
subplot(2,1,1);
plot(t,abs(x(:,1)));
xlabel("$t$","Interpreter","latex");
ylabel("$|A|$","Interpreter","latex");
subplot(2,1,2);
plot(t,abs(x(:,2)));
xlabel("$t$","Interpreter","latex");
ylabel("$|B|$","Interpreter","latex");
set(gca,'FontSize',20)
saveas(gcf,'amplitudes-eps01.png');


%% quasi-periodic

eps = 0.09;
Tend = 100000;
omegaA = 5;

vecField = @(t,X) [eps^2*X(1).*(1 - (1 + 5*1i)*abs(X(1)).^2 + (1+5*1i)*abs(X(2)).^2) + eps^4*X(2).*abs(X(1)).^4.*exp(1i * t);...
                    eps^2*X(2).*(1 - abs(X(2)).^2 - (2-1i)*abs(X(1)).^2) + eps^4*X(1).*abs(X(1)).^4.*exp(-1i * t)];


Xinit = [1,0];
[t,x] = ode45(vecField,[0,Tend],Xinit);

figure;
plot(angle(exp(1i*t)),angle(x(:,1).*exp(-1i*omegaA*eps^2*t)),'.');
xlabel("$\theta$","FontSize",28,"Interpreter","latex");
ylabel("$\phi_A$","FontSize",28,"Interpreter","latex");
title("$\varepsilon = 0.09$","FontSize",28,"Interpreter","latex");
ylim([-pi,pi]);
xlim([-pi,pi]);
set(gca,'FontSize',20)
saveas(gcf,"quasiperiodic.png");

figure;
subplot(2,1,1);
plot(t,abs(x(:,1)));
xlabel("$t$","Interpreter","latex");
ylabel("$|A|$","Interpreter","latex");
subplot(2,1,2);
plot(t,abs(x(:,2)));
xlabel("$t$","Interpreter","latex");
ylabel("$|B|$","Interpreter","latex");
set(gca,'FontSize',20)
saveas(gcf,'amplitudes-eps009.png');