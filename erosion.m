
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
    L = r*y' - y*r';
% %     hs = (r0 - rlp)*ones(1,N);  % additive, activity-dependent rate control
% %     hs = (1 - r'*r/sqrt(N))*W;    % multiplicative, averaged rate control

%%% Homeostatic rules
    hs = B - tanh(x - xlp)*tanh(x)';
%     hs = (r0 - r)*ones(N,1)'*W;
%     hs = (r0 - r)*r'*W;
%     hs = -beta*W;

    noise = (1*randn(N,N))/sqrt(N);
    W = W + etaS*(noise + hs)*dt;
    

    if ~mod(i-1,CalcEvery)
        W_all(:,:,(i-1)/CalcEvery + 1) = W;
        hs_all(:,:,(i-1)/CalcEvery + 1) = hs/2 + hs'/2;
        noise_all(:,:,(i-1)/CalcEvery + 1) = noise/2 + noise'/2;
        S = hs/2 + hs'/2;
        S = S(:);
        Ps((i-1)/CalcEvery + 1) = hs(:)'*S/(norm(S));
        A = hs/2 - hs'/2;
        A = A(:);
        Pa((i-1)/CalcEvery + 1) = hs(:)'*A/(norm(A));
        S = noise/2 + noise'/2;
        S = S(:);
        Psnoise((i-1)/CalcEvery + 1) = noise(:)'*S/(norm(S));
        A = noise/2 - noise'/2;
        A = A(:);
        Panoise((i-1)/CalcEvery + 1) = noise(:)'*A/(norm(A));
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


% [~,Dseq] = eigenshuffle(noise_all);
% for j=1:1:N
%     figure(202); plot(real(Dseq(j,:)),'linewidth',2,'color','k'); hold on; 
% end
% hold on;
% [~,Dseq] = eigenshuffle(hs_all);
% for j=1:1:N
%     figure(202); plot(real(Dseq(j,:)),'linewidth',2,'color','r'); hold on; 
% end

figure; plot(Ps,Pa,'o'); hold on; plot(Psnoise,Panoise,'o'); axis square;
xlabel('S'); ylabel('A');

STD = nan(N,N);
for i=1:N
    for j=1:N
        STD(i,j) = std(W_all(i,j,:));
    end
end
mean(STD(:))

figure; plot(x_all')

S = W/2 + W'/2;
A = W/2 - W'/2;
S = S(:);
A = A(:);
PsW = W(:)'*S/norm(S)
PaW = W(:)'*A/norm(A)


toc