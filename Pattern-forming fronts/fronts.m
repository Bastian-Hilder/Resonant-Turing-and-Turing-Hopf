%% script to plot the pattern interfaces found close to a 1:1 resonant Turing and Turing–Hopf instability

close all

makevideo = 0;

%% coefficients
tildec = 0.5;
gammaA = -2;
gammaB = -2;

gamma1i = 1;
gamma2i = 1;
gamma7i = 1;
gamma8i = 1;

Tend = 15;

%% vecfield

omegaA = gamma1i;
omegaB = gamma8i;

% X = (rA, rB, phiA, phiB)
vecfield = @(~,X) [-X(1) + X(1).^3 - gammaA*X(1).*X(2).^2;...
                    (1/tildec)*(X(2) - X(2).^3 + gammaB*X(2).*X(1).^2);...
                    -omegaA + gamma1i*X(1).^2 + gamma2i*X(2).^2;...
                    (1/tildec)*(-omegaB + gamma7i*X(1).^2 + gamma8i*X(2).^2)];


%% Front Turing–Hopf to Superposition

% solve ode
Xinit = [1e-5,1,0,0];

[t,X] = ode45(vecfield,[-Tend,Tend],Xinit);
% figure;
% plot(X(:,1),X(:,2));
% 
% pause();

% plot solution

% interpolate solution to uniform grid
xiGrid = linspace(-Tend,Tend,10000);
rA = interp1(t,X(:,1),xiGrid);
rB = interp1(t,X(:,2),xiGrid);

PhiA = interp1(t,X(:,3),xiGrid);
PhiB = interp1(t,X(:,4),xiGrid);

% figure;
% plot(xiGrid,rA,xiGrid,rB);
% 
% figure;
% plot(xiGrid,PhiA,xiGrid,PhiB);
% 
% pause();

% full solution

eps = 0.2;
c0 = 1;
cp = 0.5;

u = @(x,time) 2*eps*interp1(t,X(:,1),eps^2*(x-c0*time)).*cos(x+interp1(t,X(:,3),eps^2*(x-c0*time))+eps^2*omegaA*time);
v = @(x,time) 2*eps*interp1(t,X(:,2),eps^2*(x-c0*time)).*cos(x-cp*time+interp1(t,X(:,4),eps^2*(x-c0*time))+eps^2*omegaB*time);


xgrid = linspace(0,50,1000);
tgrid = linspace(0,300,1000);

% video
if makevideo
    for ii = 1:numel(tgrid)
        subplot(2,1,1);
        plot(xgrid,u(xgrid,tgrid(ii)));
        title("$u$","Interpreter","latex");
        xlim([0,50]);
        ylim([-0.5,0.5]);
    
        subplot(2,1,2);
        plot(xgrid,v(xgrid,tgrid(ii)));
        title("$v$","Interpreter","latex");
        xlim([0,50]);
        ylim([-0.5,0.5]);
        pause(0.0001)
    end
end



% space-time plot (top-down view)

figure
[Xmesh,T] = meshgrid(linspace(0,50,1000),linspace(0,300,1000));
U = u(Xmesh,T);
V = v(Xmesh,T);

hold on
s = pcolor(Xmesh,T,U);
grid off;
s.EdgeColor="none";
colormap(flipud(gray))
xlabel("$\tilde{\xi}$","FontSize",35,"Interpreter","latex")
ylabel("$t$","FontSize",35,"Interpreter","latex")
title("$u$","FontSize",35,"Interpreter","latex")
xlim([0,50])
ylim([0,300])
pbaspect([4/2 1 1])
set(gca,'FontSize',24)

set(gcf,"Position",[100, 100, 800, 420])

saveas(gcf,"u-TH-to-Superpos.png")

figure

hold on
s = pcolor(Xmesh,T,V);
grid off;
s.EdgeColor="none";
colormap(flipud(gray))
xlabel("$\tilde{\xi}$","FontSize",35,"Interpreter","latex")
ylabel("$t$","FontSize",35,"Interpreter","latex")
title("$v$","FontSize",35,"Interpreter","latex")
xlim([0,50])
ylim([0,300])
pbaspect([4/2 1 1])
set(gca,'FontSize',24)

set(gcf,"Position",[100, 100, 800, 420])

saveas(gcf,"v-TH-to-Superpos.png")


%% Front Turing to Trivial

timeEnd = 350;
Xend = 50;

% solve ode
Xinit = [1-1e-7,0,0,0];

[t,X] = ode45(vecfield,[-Tend,Tend],Xinit);
% figure;
% plot(X(:,1),X(:,2));
% 
% pause();

% plot solution

% interpolate solution to uniform grid
xiGrid = linspace(-Tend,Tend,10000);
rA = interp1(t,X(:,1),xiGrid);
rB = interp1(t,X(:,2),xiGrid);

PhiA = interp1(t,X(:,3),xiGrid);
PhiB = interp1(t,X(:,4),xiGrid);

% figure;
% plot(xiGrid,rA,xiGrid,rB);
% 
% figure;
% plot(xiGrid,PhiA,xiGrid,PhiB);
% 
% pause();

% full solution

eps = 0.2;
c0 = 1;
cp = 0.5;

u = @(x,time) 2*eps*interp1(t,X(:,1),eps^2*(x-c0*time)).*cos(x+interp1(t,X(:,3),eps^2*(x-c0*time))+eps^2*omegaA*time);
v = @(x,time) 2*eps*interp1(t,X(:,2),eps^2*(x-c0*time)).*cos(x-cp*time+interp1(t,X(:,4),eps^2*(x-c0*time))+eps^2*omegaB*time);


xgrid = linspace(0,Xend,1000);
tgrid = linspace(0,timeEnd,1000);

% video
if makevideo
    for ii = 1:numel(tgrid)
        subplot(2,1,1);
        plot(xgrid,u(xgrid,tgrid(ii)));
        title("$u$","Interpreter","latex");
        xlim([0,50]);
        ylim([-0.5,0.5]);
    
        subplot(2,1,2);
        plot(xgrid,v(xgrid,tgrid(ii)));
        title("$v$","Interpreter","latex");
        xlim([0,Xend]);
        ylim([-0.5,0.5]);
        pause(0.0001)
    end
end



% space-time plot (top-down view)

figure
[Xmesh,T] = meshgrid(linspace(0,Xend,1000),linspace(0,timeEnd,1000));
U = u(Xmesh,T);
V = v(Xmesh,T);

hold on
s = pcolor(Xmesh,T,U);
grid off;
s.EdgeColor="none";
colormap(flipud(gray))
xlabel("$\tilde{\xi}$","FontSize",35,"Interpreter","latex")
ylabel("$t$","FontSize",35,"Interpreter","latex")
title("$u$","FontSize",35,"Interpreter","latex")
xlim([0,Xend])
ylim([0,timeEnd])
pbaspect([4/2 1 1])
set(gca,'FontSize',24)

set(gcf,"Position",[100, 100, 800, 420])

saveas(gcf,"u-T-to-trivial.png")

figure

hold on
s = pcolor(Xmesh,T,V);
grid off;
s.EdgeColor="none";
colormap(flipud(gray))
xlabel("$\tilde{\xi}$","FontSize",35,"Interpreter","latex")
ylabel("$t$","FontSize",35,"Interpreter","latex")
title("$v$","FontSize",35,"Interpreter","latex")
xlim([0,Xend])
ylim([0,timeEnd])
pbaspect([4/2 1 1])
set(gca,'FontSize',24)

set(gcf,"Position",[100, 100, 800, 420])

saveas(gcf,"v-T-to-trivial.png")

%% Front Turing–Hopf to Trivial

% change speed and update vectorfield
tildec = -0.5;

% X = (rA, rB, phiA, phiB)
vecfield = @(~,X) [-X(1) + X(1).^3 - gammaA*X(1).*X(2).^2;...
                    (1/tildec)*(X(2) - X(2).^3 + gammaB*X(2).*X(1).^2);...
                    -omegaA + gamma1i*X(1).^2 + gamma2i*X(2).^2;...
                    (1/tildec)*(-omegaB + gamma7i*X(1).^2 + gamma8i*X(2).^2)];

timeEnd = 350;
Xend = 50;

% solve ode
Xinit = [0,1-1e-15,0,0];

[t,X] = ode45(vecfield,[-Tend,Tend],Xinit);
% figure;
% plot(X(:,1),X(:,2));

% plot solution

% interpolate solution to uniform grid
xiGrid = linspace(-Tend,Tend,10000);
rA = interp1(t,X(:,1),xiGrid);
rB = interp1(t,X(:,2),xiGrid);

PhiA = interp1(t,X(:,3),xiGrid);
PhiB = interp1(t,X(:,4),xiGrid);

% figure;
% plot(xiGrid,rA,xiGrid,rB);
% 
% figure;
% plot(xiGrid,PhiA,xiGrid,PhiB);
% 
% pause();

% full solution

eps = 0.2;
c0 = 1;
cp = 0.5;

u = @(x,time) 2*eps*interp1(t,X(:,1),eps^2*(x-c0*time)).*cos(x+interp1(t,X(:,3),eps^2*(x-c0*time))+eps^2*omegaA*time);
v = @(x,time) 2*eps*interp1(t,X(:,2),eps^2*(x-c0*time)).*cos(x-cp*time+interp1(t,X(:,4),eps^2*(x-c0*time))+eps^2*omegaB*time);


xgrid = linspace(0,Xend,1000);
tgrid = linspace(0,timeEnd,1000);

% video
if makevideo
    for ii = 1:numel(tgrid)
        subplot(2,1,1);
        plot(xgrid,u(xgrid,tgrid(ii)));
        title("$u$","Interpreter","latex");
        xlim([0,50]);
        ylim([-0.5,0.5]);
    
        subplot(2,1,2);
        plot(xgrid,v(xgrid,tgrid(ii)));
        title("$v$","Interpreter","latex");
        xlim([0,Xend]);
        ylim([-0.5,0.5]);
        pause(0.0001)
    end
end



% space-time plot (top-down view)

figure
[Xmesh,T] = meshgrid(linspace(0,Xend,1000),linspace(0,timeEnd,1000));
U = u(Xmesh,T);
V = v(Xmesh,T);

hold on
s = pcolor(Xmesh,T,U);
grid off;
s.EdgeColor="none";
colormap(flipud(gray))
xlabel("$\tilde{\xi}$","FontSize",35,"Interpreter","latex")
ylabel("$t$","FontSize",35,"Interpreter","latex")
title("$u$","FontSize",35,"Interpreter","latex")
xlim([0,Xend])
ylim([0,timeEnd])
pbaspect([4/2 1 1])
set(gca,'FontSize',24)

set(gcf,"Position",[100, 100, 800, 420])

saveas(gcf,"u-TH-to-trivial.png")

figure

hold on
s = pcolor(Xmesh,T,V);
grid off;
s.EdgeColor="none";
colormap(flipud(gray))
xlabel("$\tilde{\xi}$","FontSize",35,"Interpreter","latex")
ylabel("$t$","FontSize",35,"Interpreter","latex")
title("$v$","FontSize",35,"Interpreter","latex")
xlim([0,Xend])
ylim([0,timeEnd])
pbaspect([4/2 1 1])
set(gca,'FontSize',24)

set(gcf,"Position",[100, 100, 800, 420])

saveas(gcf,"v-TH-to-trivial.png")


%% Front Turing to Turing–Hopf

% change speed and update vectorfield
tildec = 0.5;
gammaA = 1;
gammaB = 2;

% X = (rA, rB, phiA, phiB)
vecfield = @(~,X) [-X(1) + X(1).^3 - gammaA*X(1).*X(2).^2;...
                    (1/tildec)*(X(2) - X(2).^3 + gammaB*X(2).*X(1).^2);...
                    -omegaA + gamma1i*X(1).^2 + gamma2i*X(2).^2;...
                    (1/tildec)*(-omegaB + gamma7i*X(1).^2 + gamma8i*X(2).^2)];

timeEnd = 350;
Xend = 50;

% solve ode
Xinit = [1,1e-15,0,0];

[t,X] = ode45(vecfield,[-Tend,Tend],Xinit);
% figure;
% plot(X(:,1),X(:,2));

% plot solution

% interpolate solution to uniform grid
xiGrid = linspace(-Tend,Tend,10000);
rA = interp1(t,X(:,1),xiGrid);
rB = interp1(t,X(:,2),xiGrid);

PhiA = interp1(t,X(:,3),xiGrid);
PhiB = interp1(t,X(:,4),xiGrid);

% figure;
% plot(xiGrid,rA,xiGrid,rB);
% 
% figure;
% plot(xiGrid,PhiA,xiGrid,PhiB);
% 
% pause();

% full solution

eps = 0.2;
c0 = 1;
cp = 0.5;

u = @(x,time) 2*eps*interp1(t,X(:,1),eps^2*(x-c0*time)).*cos(x+interp1(t,X(:,3),eps^2*(x-c0*time))+eps^2*omegaA*time);
v = @(x,time) 2*eps*interp1(t,X(:,2),eps^2*(x-c0*time)).*cos(x-cp*time+interp1(t,X(:,4),eps^2*(x-c0*time))+eps^2*omegaB*time);


xgrid = linspace(0,Xend,1000);
tgrid = linspace(0,timeEnd,1000);

% video
if makevideo
    for ii = 1:numel(tgrid)
        subplot(2,1,1);
        plot(xgrid,u(xgrid,tgrid(ii)));
        title("$u$","Interpreter","latex");
        xlim([0,50]);
        ylim([-0.5,0.5]);
    
        subplot(2,1,2);
        plot(xgrid,v(xgrid,tgrid(ii)));
        title("$v$","Interpreter","latex");
        xlim([0,Xend]);
        ylim([-0.5,0.5]);
        pause(0.0001)
    end
end



% space-time plot (top-down view)

figure
[Xmesh,T] = meshgrid(linspace(0,Xend,1000),linspace(0,timeEnd,1000));
U = u(Xmesh,T);
V = v(Xmesh,T);

hold on
s = pcolor(Xmesh,T,U);
grid off;
s.EdgeColor="none";
colormap(flipud(gray))
xlabel("$\tilde{\xi}$","FontSize",35,"Interpreter","latex")
ylabel("$t$","FontSize",35,"Interpreter","latex")
title("$u$","FontSize",35,"Interpreter","latex")
xlim([0,Xend])
ylim([0,timeEnd])
pbaspect([4/2 1 1])
set(gca,'FontSize',24)

set(gcf,"Position",[100, 100, 800, 420])

saveas(gcf,"u-T-to-TH.png")

figure

hold on
s = pcolor(Xmesh,T,V);
grid off;
s.EdgeColor="none";
colormap(flipud(gray))
xlabel("$\tilde{\xi}$","FontSize",35,"Interpreter","latex")
ylabel("$t$","FontSize",35,"Interpreter","latex")
title("$v$","FontSize",35,"Interpreter","latex")
xlim([0,Xend])
ylim([0,timeEnd])
pbaspect([4/2 1 1])
set(gca,'FontSize',24)

set(gcf,"Position",[100, 100, 800, 420])

saveas(gcf,"v-T-to-TH.png")


%% Front Superposition to Turing

% change speed and update vectorfield
tildec = -2;
gammaA = -0.1;
gammaB = 0.1;
Tend = 25;

% X = (rA, rB, phiA, phiB)
vecfield = @(~,X) [-X(1) + X(1).^3 - gammaA*X(1).*X(2).^2;...
                    (1/tildec)*(X(2) - X(2).^3 + gammaB*X(2).*X(1).^2);...
                    -omegaA + gamma1i*X(1).^2 + gamma2i*X(2).^2;...
                    (1/tildec)*(-omegaB + gamma7i*X(1).^2 + gamma8i*X(2).^2)];

timeEnd = 500;
Xend = 50;


% solve ode in reverse time
Xinit = [1,1e-9,0,0];

[t,X] = ode45(@(t,X) -vecfield(t,X),[-Tend,Tend],Xinit);
% figure;
% plot(X(:,1),X(:,2));

% reverse solution
t = -flip(t);
X = flip(X,1);

% plot solution

% interpolate solution to uniform grid
xiGrid = linspace(-Tend,Tend,10000);
rA = interp1(t,X(:,1),xiGrid);
rB = interp1(t,X(:,2),xiGrid);

PhiA = interp1(t,X(:,3),xiGrid);
PhiB = interp1(t,X(:,4),xiGrid);

% figure;
% plot(xiGrid,rA,xiGrid,rB);
% 
% figure;
% plot(xiGrid,PhiA,xiGrid,PhiB);
% 
% pause();

% full solution

eps = 0.2;
c0 = 1;
cp = 0.5;

u = @(x,time) 2*eps*interp1(t,X(:,1),eps^2*(x-c0*time)).*cos(x+interp1(t,X(:,3),eps^2*(x-c0*time))+eps^2*omegaA*time);
v = @(x,time) 2*eps*interp1(t,X(:,2),eps^2*(x-c0*time)).*cos(x-cp*time+interp1(t,X(:,4),eps^2*(x-c0*time))+eps^2*omegaB*time);


xgrid = linspace(0,Xend,1000);
tgrid = linspace(0,timeEnd,1000);

% video
if makevideo
    for ii = 1:numel(tgrid)
        subplot(2,1,1);
        plot(xgrid,u(xgrid,tgrid(ii)));
        title("$u$","Interpreter","latex");
        xlim([0,50]);
        ylim([-0.5,0.5]);
    
        subplot(2,1,2);
        plot(xgrid,v(xgrid,tgrid(ii)));
        title("$v$","Interpreter","latex");
        xlim([0,Xend]);
        ylim([-0.5,0.5]);
        pause(0.0001)
    end
end



% space-time plot (top-down view)

figure
[Xmesh,T] = meshgrid(linspace(0,Xend,1000),linspace(0,timeEnd,1000));
U = u(Xmesh,T);
V = v(Xmesh,T);

hold on
s = pcolor(Xmesh,T,U);
grid off;
s.EdgeColor="none";
colormap(flipud(gray))
xlabel("$\tilde{\xi}$","FontSize",35,"Interpreter","latex")
ylabel("$t$","FontSize",35,"Interpreter","latex")
title("$u$","FontSize",35,"Interpreter","latex")
xlim([0,Xend])
ylim([0,timeEnd])
pbaspect([4/2 1 1])
set(gca,'FontSize',24)

set(gcf,"Position",[100, 100, 800, 420])

saveas(gcf,"u-Superpos-to-T.png")

figure

hold on
s = pcolor(Xmesh,T,V);
grid off;
s.EdgeColor="none";
colormap(flipud(gray))
xlabel("$\tilde{\xi}$","FontSize",35,"Interpreter","latex")
ylabel("$t$","FontSize",35,"Interpreter","latex")
title("$v$","FontSize",35,"Interpreter","latex")
xlim([0,Xend])
ylim([0,timeEnd])
pbaspect([4/2 1 1])
set(gca,'FontSize',24)

set(gcf,"Position",[100, 100, 800, 420])

saveas(gcf,"v-Superpos-to-T.png")