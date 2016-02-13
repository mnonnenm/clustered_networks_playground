clearvars;
clc

% set core simulation parameters:

pars.T = 500; % simulation length in ms

pars.N = 100; % number of neurons (total) - 80% will be made excitatory

pars.K = 2;   % number of modules

pars.ps.pEEcr = 0.05 / (pars.K-1); % cross-cluster 

pars.ps.pEIcr = 0.10 / (pars.K-1);


pars.ampa_ratio = 0.3;

pars.chunksize = 50;

sf = '/home/mackelab/Desktop/Projects/Stitching/results/cosyne_poster/gb_net/test';

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
sf = '/home/mackelab/Desktop/Projects/Stitching/results/cosyne_poster/gb_net/cluster/gb_nets';
sf = [sf, '/test_large3'];
if strcmp(sf, '')
    if_save = false;
else
    if_save = true;
end


spikes_all = [];
if if_save
    for i = 1:150 %floor( (pars.T/pars.dt)/pars.chunksize )
        load([sf, '_ck', num2str(i)])
        disp(['loaded ', num2str(i)])
        spikes_all = [spikes_all; spikes];
    end
    
    spikes = [spikes_all; spikes_last];
    clear spikes_all spikes_last

    save([sf, '_all'], 'spikes', 'pars') 
else
    spikes = spikes_last;
    clear spikes_last
end



%%
% close all
% figure
% for i = 1:4
%     subplot(2,2,i)
%     x = cellfun(@(s) {s(:,i)}, setup.SOL); x = [x{:}];
%     plot(setup.t0, x);
% end

figure, hold on
for j = 1:pars.Ne
    flag = spikes(:,j);
    ts = find(flag);
    plot(setup.t0(ts), j * ones(size(ts)),'.k');
end
for j = pars.Ne + (1:pars.Ni)
    flag = spikes(:,j);
    ts = find(flag);
    plot(setup.t0(ts), j * ones(size(ts)),'.r');
end

xlim([setup.t0(1), setup.t0(end)])
ylim([0, pars.N + 1]);
clsz_e = pars.Ne/pars.K;
for k = 1:pars.K     
    set(refline(0, k * clsz_e + 0.5),'Color','b','LineStyle',':')
end
xlabel('time [ms]')
ylabel('neuron id')
hold off

% figure;
% clrs = hsv(5);
% idx = sort( pars.Ne + randsample(pars.Ni, 5) );
% for i = 1:5
%     tmp = setup.SOL{idx(i)};
%     plot(tmp(:,6), 'color', clrs(i,:));
%     hold on
% end
% hold off

