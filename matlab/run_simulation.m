function run_simulation(T, N, K, pEEcr, pEIcr, ampa_ratio, chunksize, sf)

% set core simulation parameters:

if nargin < 1 || isempty(T)
    pars.T = 500; % simulation length in ms
else
    pars.T = T;
end

if nargin < 2 || isempty(N)
    pars.N = 100; % number of neurons (total) - 80% will be made excitatory
else 
    pars.N = N;    
end

if nargin < 3 || isempty(K)
    pars.K = 2;   % number of modules
else
    pars.K = K;
end

if nargin < 4 || isempty(pEEcr)
    pars.ps.pEEcr = 0.05 / (pars.K-1); % cross-cluster 
else
    pars.ps.pEEcr = pEEcr;
end

if nargin < 5 || isempty(pEIcr)
    pars.ps.pEIcr = 0.10 / (pars.K-1);
else
    pars.ps.pEIcr = pEIcr;
end

if nargin < 6 || isempty(ampa_ratio)
    pars.ampa_ratio = 0.3;
else
    pars.ampa_ratio = ampa_ratio;    
end

if nargin < 7 || isempty(chunksize)
    pars.chunksize = 50;
else
    pars.chunksize = chunksize;
end

if nargin < 8 || isempty(sf)
    sf = '';
end

% set more static parameters:

pars.dt = 1;
    
pars.Ne = fix(0.8 * pars.N);
pars.Ni = pars.N - pars.Ne;

pars.syn_params = struct(...
    'Tampa', 5,'Eampa', 0,...
    'Tnmda1', 80,'Tnmda2', 20,'Enmda', 0,...
    'Tgaba', 6,'Egaba', -70);
pars.params_e = struct(...
    'mem', struct('Er',-60,'Et',-50,'C', 100,...
    'k', 3,'a', 0.01,'b', 5,'c',-60,'d', 400),...
    'syn', pars.syn_params);
pars.params_i = struct(...
    'mem', struct('Er',-55,'Et',-40,'C', 20,...
    'k', 1,'a', 0.15,'b', 8,'c',-55,'d', 200),...  
    'syn', pars.syn_params);

pars.ps.pEE = 0.60;
pars.ps.pEI = 0.60;
pars.ps.pIE = 0.40;
pars.ps.pII = 0.20;

pars.rate_e = 0.40;
pars.rate_i = 0.30;

pars.jXE2 = 1.2;
pars.jXE5 = 0.3;
pars.jXI  = 1.5;

%% Run simulation

[spikes_last, setup] = run_network(pars, sf);

%% Assemble data

if strcmp(sf, '')
    if_save = false;
else
    if_save = true;
end

spikes_all = [];
if if_save
    
    for i = 1:floor( (pars.T/pars.dt)/pars.chunksize )
        load([sf, '_ck', num2str(i)])
        disp(['loaded ', num2str(i)])
        spikes_all = [spikes_all; spikes];
    end
    
    spikes = [spikes_all; spikes_last];
    clear spikes_all spikes_last
    setup = rmfield(setup, 'SOL');
    save([sf, '_all'], 'spikes', 'pars', 'setup') 
    
else
    spikes = spikes_last;
    clear spikes_last
end

