function big_sim(sim_index, seed)

%% set up for storing
dstdir = fullfile(pwd, sprintf('subresults%d', sim_index));
if ~exist(dstdir,'dir')
    mkdir(dstdir)
end

%% initialize random stream
rng(seed);

%%    
syn_params = struct(...
    'Tampa', 5,'Eampa', 0,...
    'Tnmda1', 80,'Tnmda2', 20,'Enmda', 0,...
    'Tgaba', 6,'Egaba', -70);
params_e = struct(...
    'mem', struct('Er',-60,'Et',-50,'C', 100,...
    'k', 3,'a', 0.01,'b', 5,'c',-60,'d', 400),...
    'syn', syn_params);
params_i = struct(...
    'mem', struct('Er',-55,'Et',-40,'C', 20,...
    'k', 1,'a', 0.15,'b', 8,'c',-55,'d', 200),...  
    'syn', syn_params);
    
rate_e = 0.40;
rate_i = 0.30;

M = 1250;
Me = fix(0.8 * M);
Mi = M - Me;
neur = cell(M, 1);
for i = 1:M
    if i <= Me
        neur{i} = make_neuron(params_e);
    else
        neur{i} = make_neuron(params_i);
    end
end
neur = [neur{:}];

%%
N = 6000;
burnin = 100;
t0 = 0:N;
Y0 = arrayfun(@(m) {[-65 + rand(1) * 10, 10 + 20*[1, 1, (m <= Me), (m <= Me)].* rand(1,4), (1 + randn(1))]'}, 1:M);
for i = 1:M
    opts{i} = odeset(...
        'Events',neur(i).events,...
        'RelTol',1e-3,'AbsTol',1e-3,...
        'InitialStep', t0(2) - t0(1),'MaxStep', t0(2) - t0(1));
end

%% Generate external input
spikes_ext_e = sparse(poissrnd(rate_e, N, M));
spikes_ext_i = sparse(poissrnd(rate_i, N, M));

load(fullfile(pwd,'connectivity.mat'))
ampa_ratio = 0.5;
nmda_ratio = 1 - ampa_ratio;

%% Simulate the system
vr = [repmat(params_e.mem.c, Me, 1); repmat(params_i.mem.c, Mi, 1)];
ur = [repmat(params_e.mem.d, Me, 1); repmat(params_i.mem.d, Mi, 1)];

% h = waitbar(0,'Computing....');
Geb = Ge; Gib = Gi;
K=10;
clustersize = Me/K;
for k = 1:K
    idx = (k-1)*clustersize+(1:clustersize);
    Geb(idx,idx) = 0;       
end
spikes = zeros(burnIn, M);
t_on = linspace(0, burnin, K+1); t_on = t_on(2:end);
for i = 1:burnIn
    % update membrane potentials
    tic
    t_1 = i;
    t_2 = i+1;
    for j = 1:M
        [~, Y, TE, YE] = ode23(neur(j).odefun, [t_1, t_2], Y0{j}, opts{j});
        if ~isempty(TE)
            spikes(i,j) = 1;
            Y = YE(end,:)';
            Y(1) = vr(j);
            Y(6) = Y(6) + ur(j);
            neur(j).set_spike(TE);
            [~, Y] = ode23(neur(j).odefun, [TE, t_2], Y, opts{j});
        end
        Y0{j} = Y(end,:);
    end
    toc
    % propagate effect of spikes
    s = spikes(i,:);
    i_e = s(1:Me) * Geb;
    i_i = s(Me+1:end) * Gib;
    for j = 1:M
        % update for next iteration
        Y0{j}(3) = Y0{j}(3) + i_i(j);
        if j <= Me
            Y0{j}(2) = Y0{j}(2) + ampa_ratio * i_e(j);
        else
            Y0{j}(2) = Y0{j}(2) + i_e(j);
        end
        
        if j <= Me            
            Y0{j}(2) = Y0{j}(2) + 1.2 * spikes_ext_e(i,j);
            Y0{j}(5) = Y0{j}(5) + 0.3 * spikes_ext_e(i,j);
            Y0{j}(3) = Y0{j}(3) + 1.5 * spikes_ext_i(i,j);
            % nmda only on excitatory neurons
            Y0{j}(5) = Y0{j}(5) + nmda_ratio * i_e(j);
        end
        % saturate conductances
        Y0{j}(2:end-2) = min(Y0{j}(2:end-2), [100, 100, 50]);
    end
    if mod(i, 10) == 0
        clc
        disp(i)
    end
    if ismember(i, t_on) 
        k = find(i==t_on);
        idx = (k-1)*clustersize+(1:clustersize);
        Geb(idx,idx) = Ge(idx,idx);        
        imagesc(Geb)
        pause;
    end
end

disp('finished "burnin" phase!')
disp('starting main run')

spikes = zeros(1000, M);
local_index = 0;
for i = 1:length(t0)-1
    local_index = local_index + 1;
    % update membrane potentials
    tic
    t_1 = t0(i);
    t_2 = t0(i+1);
    for j = 1:M
        [~, Y, TE, YE] = ode23(neur(j).odefun, [t_1, t_2], Y0{j}, opts{j});
        if ~isempty(TE)
            spikes(local_index,j) = 1;
            Y = YE(end,:)';
            Y(1) = vr(j);
            Y(6) = Y(6) + ur(j);
            neur(j).set_spike(TE);
            [~, Y] = ode23(neur(j).odefun, [TE, t_2], Y, opts{j});
        end
        Y0{j} = Y(end,:);
    end
    toc
    % propagate effect of spikes
    s = spikes(local_index,:);
    i_e = s(1:Me) * Ge;
    i_i = s(Me+1:end) * Gi;
    for j = 1:M
        % update for next iteration
        Y0{j}(3) = Y0{j}(3) + i_i(j);
        if j <= Me
            Y0{j}(2) = Y0{j}(2) + ampa_ratio * i_e(j);
        else
            Y0{j}(2) = Y0{j}(2) + i_e(j);
        end
        
        if j <= Me            
            Y0{j}(2) = Y0{j}(2) + 1.2 * spikes_ext_e(i,j);
            Y0{j}(5) = Y0{j}(5) + 0.3 * spikes_ext_e(i,j);
            Y0{j}(3) = Y0{j}(3) + 1.5 * spikes_ext_i(i,j);
            % nmda only on excitatory neurons
            Y0{j}(5) = Y0{j}(5) + nmda_ratio * i_e(j);
        end
        % saturate conductances
        Y0{j}(2:end-2) = min(Y0{j}(2:end-2), [100, 100, 50]);
    end
    if mod(i, 10) == 0
        clc
        disp(i)
    end
    if mod(i, 1000) == 0
        % save partial data
        spikes = sparse(spikes); %#ok<NASGU>
        index = fix(i / 1000);
        filename = fullfile(dstdir, sprintf('chunk_%02d.mat', index));
        save(filename, 'spikes', 'seed');
        spikes = zeros(1000, M);
        local_index = 0;
    end
end

end