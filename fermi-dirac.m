clear all
close all

%% Q1 - Q2
qmax = 8;
N = 4;
dq = 3; % excited state to plot

combos = combnk(1:qmax,N);
E = sum(combos,2); % E / epsilon
E0 = min(E);

% Plot # of microstates vs. q
h = histogram(E)
xlabel('E_{total}/ \epsilon')
ylabel('# of microstates')
axis([0,30,0,inf])

dist = zeros(qmax,1); %padding a few more zeros on the right

x = find(E == E0 + dq); % find the combos with the right energy
y = combos(x,:);
for j = 1:qmax
    z = find(y==j);
    dist(j) = length(z)/length(x);
end

% Plot probability of a state being occupied vs. q
figure
plot(0:length(dist)+1,[1;dist;0],'*')
xlabel('q')
ylabel('n_{FE}')

hold on

% Fit to FD distribution
a = 0:10;
plot(a,1./(1+exp((a-4.5)/1.5)))

%%
qmax = 20; % number of energy levels
N = 10; % number of electrons
dq = 10; % excited state to plot

combos = combnk(1:qmax,N);
E = sum(combos,2); % E / epsilon
E0 = min(E);

% Plot # of microstates vs. q
h = histogram(E)
xlabel('E_{total}/ \epsilon')
ylabel('# of microstates')
axis([0,100,0,inf])

dist = zeros(qmax,1); %padding a few more zeros on the right

x = find(E == E0 + dq); % find the combos with the right energy
y = combos(x,:);
for j = 1:qmax
    z = find(y==j);
    dist(j) = length(z)/length(x);
end

%% Plot probability of a state being occupied vs. q
figure
plot(0:length(dist)+1,[1;dist;0],'*')
xlabel('q')
ylabel('n_{FE}')
hold on

% Fit to FD distribution
a = 0:25;
plot(a,1./(1+exp((a-10.5)/3)))

%% Entropy
S =  log(h.Values); % S/k
figure
plot(h.BinWidth*(1:length(S)),S)
axis([0,50,0,inf])
xlabel('q')
ylabel('S/k')

% Calculate the slope around dq
m = (S(dq+1)-S(dq-1))/2;
disp(['T = ', num2str(1/m)])