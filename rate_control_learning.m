
%  Simulate rate network with rate-control homeostasis
%%%%  Dynamically learn one memory item

clear all; clc;
close all;

load('data'); 
rng(data.rcseed);

N = 128;
dt = .1;
etaS = .01;
g = 1.0;
tau = 20;
tauy = 100;
Nmem = 1;

TotalSteps = 3000/dt - 1;
CalcEvery = 10/dt;
Nsteps = (TotalSteps+1)/CalcEvery;

x = randn(N,1);
x_all = nan(N,TotalSteps);
xlp = randn;
W_all = nan(N,N,Nsteps);
hs_all = nan(N,N,Nsteps);
noise_all = nan(N,N,Nsteps);
W = 2*randn(N)/sqrt(N);
r0 = 2*(rand(N,1) - .5);
zeta = .01;
y = randn(N,1);

H = sign(randn(N,N)/sqrt(N));
u = H(1:Nmem,:)'/sqrt(N);
v = H(Nmem+1:2*Nmem,:)'/sqrt(N);
input1 = zeros(N,TotalSteps);
input2 = zeros(N,TotalSteps);
base = 500/dt;
inLen = 100/dt;
input = 0;

f = waitbar(0,'Simulating, please wait...');

%%%% Construct input signal for learning
for i = 1:Nmem
    input1(:,base+1:base+inLen) = repmat(u(:,i),1,inLen);
    input2(:,base+1:base+inLen) = repmat(v(:,i),1,inLen);
end

%%%% Evolve network
for i=1:TotalSteps
   
    x_all(:,i) = x;
    r = tanh(g*x);
    xlp = ((-xlp + x/5e-2)/tau)*dt;
    y = y + (r - y)*dt/tauy; 
    input = input + (-zeta*input + randn*input1(:,i) + randn*input2(:,i))*dt;
    
    x = x + (-x + W*r + 10*input)*dt;
    L = r*y' - y*r';
    hs = (r0 - r)*r'*W;

    noise = (1*randn(N,N))/sqrt(N);
    W = W + etaS*(noise + hs + L)*dt;
    
    if ~mod(i-1,CalcEvery)
        W_all(:,:,(i-1)/CalcEvery + 1) = W;
        waitbar((i-1)/TotalSteps,f,'Simulating, please wait...');
    end
end

%%%% Compute eigenspectrum of W over time
[Vseq,Dseq] = eigenshuffle(W_all);
[~,I] = sort((imag(Dseq(:,end))),'descend');
Dsorted = Dseq(I,end);
taxis = [1:1:Nsteps]*CalcEvery*dt;

%%%% Plot spectrum
for j=1:1:N
    figure(100); subplot(2,1,1); plot(taxis,real(Dseq(I(j),:)),'linewidth',2,'color',[1 1 1]*.8); hold on;
    figure(100); subplot(2,1,2); plot(taxis,imag(Dseq(I(j),:)),'linewidth',2,'color',[1 1 1]*.8); hold on;
end
subplot(2,1,2); plot(taxis,imag(Dseq(I(1),:)),'linewidth',2,'color',[.47 .67 .19]); hold on;
plot(taxis,imag(Dseq(I(end),:)),'linewidth',2,'color',[.47 .67 .19]);
xlabel('time'); ylabel('Im(\lambda)'); set(gca,'fontsize',14); 
subplot(2,1,1); plot(taxis,real(Dseq(I(1),:)),'linewidth',2,'color',[.47 .67 .19]); hold on;
plot(taxis,real(Dseq(I(end),:)),'linewidth',2,'color',[.47 .67 .19]); 
set(gca,'fontsize',14);ylabel('Re(\lambda)');

close(f);
