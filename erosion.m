
clear all; clc;
close all;

tic

N = 128;
dt = .1;
etaS = .01;
g = 1.0;
tau = 20;
tauy = 50;
beta = .1;

TotalSteps = 10000/dt - 1;
CalcEvery = 100/dt;
Nsteps = (TotalSteps+1)/CalcEvery;

x = randn(N,1);
x_all = nan(N,TotalSteps);
xlp = randn;
W_all = nan(N,N,Nsteps);
hs_all = nan(N,N,Nsteps);
noise_all = nan(N,N,Nsteps);
B = .5*eye(N);
W = 2*randn(N)/sqrt(N);
r0 = 2*(rand(N,1) - .5);
Ps = zeros(Nsteps,1);
Pa = zeros(Nsteps,1);
Psnoise = zeros(Nsteps,1);
Panoise = zeros(Nsteps,1);
y = randn(N,1);

for i=1:TotalSteps
   
    x_all(:,i) = x;
    r = tanh(g*x);
    xlp = ((-xlp + x/5e-2)/tau)*dt;
    y = y + (r - y)*dt/tauy;
    x = x + (-x + W*r)*dt;
%%% Homeostatic rules
    hs = B - tanh(x - xlp)*tanh(x)';
%     hs = (r0 - r)*r'*W;
%     hs = -beta*W;

    noise = (1*randn(N,N))/sqrt(N);
    W = W + etaS*(noise + hs)*dt;
    
    if ~mod(i-1,CalcEvery)
        W_all(:,:,(i-1)/CalcEvery + 1) = W;
        hs_all(:,:,(i-1)/CalcEvery + 1) = hs/2 + hs'/2;
        noise_all(:,:,(i-1)/CalcEvery + 1) = noise/2 + noise'/2;
        disp((i-1)/TotalSteps);
    end
    
    if i == 3000/dt
        u = randn(N,1)/sqrt(N);
        v = randn(N,1)/sqrt(N);
        W = W + 5*(u*v' - v*u');
%         W = W + 5*(u*u');
    end

end

[Vseq,Dseq] = eigenshuffle(W_all);
[~,I] = sort((imag(Dseq(:,31))),'descend');
figure(101); plot(Dseq(:,end),'o','linewidth',2); axis square; set(gca,'fontsize',18); hold on;
xlabel('Re(\lambda)'); ylabel('Im(\lambda)');

taxis = [1:1:Nsteps]*CalcEvery*dt;
for j=1:1:N
%     figure(100); subplot(2,1,1); plot(taxis,real(Dseq(j,:)),'linewidth',2,'color',[1 1 1]*.8); hold on; 
    figure(100); plot(taxis,imag(Dseq(j,:)),'linewidth',2,'color',[1 1 1]*.8); hold on;
end
plot(taxis,imag(Dseq(I(1),:)),'linewidth',2,'color',[.47 .67 .19]); hold on;
plot(taxis,imag(Dseq(I(end),:)),'linewidth',2,'color',[.47 .67 .19]); hold on;
set(gca,'fontsize',18);ylabel('Im(\lambda)');
xlabel('time'); box off;


toc
