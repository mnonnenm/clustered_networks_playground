clearvars;
close all;

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
t0 = 0:N;
ts = [];
Y0 = arrayfun(@(m) {[-65 + rand(1) * 10, 10 + 20*[1, 1, (m <= Me), (m <= Me)].* rand(1,4), (1 + randn(1))]'}, 1:M);
SOL = cellfun(@(y0) {zeros(length(t0), length(y0))}, Y0);
for i = 1:M
    SOL{i}(1,:) = Y0{i};
end
for i = 1:M
    opts{i} = odeset(...
        'Events',neur(i).events,...
        'RelTol',1e-3,'AbsTol',1e-3,...
        'InitialStep', t0(2) - t0(1),'MaxStep', t0(2) - t0(1));
end
exit = false;
Tend = 1;
spikes = zeros(N+1, M);

%% Generate external input
spikes_ext_e = sparse(poissrnd(rate_e, N, M));
spikes_ext_i = sparse(poissrnd(rate_i, N, M));

%% Define cells' position
pos_e = rand(Me, 2);
pos_i = rand(Mi, 2);

%%
figure; hold on
xy = reshape(pos_e, [], 2, K);
for k = 1:K
    plot(xy(:,1,k), xy(:,2,k),'^');
end
plot(pos_i(:,1), pos_i(:,2),'ok');
pause(0.1);

%% Build connectivity matrix: E -> E 
c_ee = cell(K);
for i = 1:K
    for j = 1:K
        if i == j
            pc = 0.30;
        else
            pc = 0.01;
        end
        c_ee{i,j} = sparse(rand(Me/K, Me/K) < pc); 
    end
end
c_ee = cell2mat(c_ee);

%% Build connectivity matrix: I -> I 
c_ii = zeros(Mi);
for i = 1:Mi-1
    for j = i+1:Mi
        dp = pos_i(i,:) - pos_i(j,:);
        r = sqrt(sum(dp.^2));
        if r < 0.2
            c_ii(i,j) = 1;
            c_ii(j,i) = 1;
        end
    end
end
c_ii = sparse(c_ii);

%% Build connectivity matrix: E -> I, I -> E
dist = 0.2;
for i = 1:Mi
    dx = pos_e(:,1) - pos_i(i,1);
    dy = pos_e(:,2) - pos_i(i,2);
    r = sqrt(dx.^2 + dy.^2);
    % inhibitory to excitatory
    pc_ie = r < dist;
    c_ie(i,:) = sparse(rand(Me, 1) < pc_ie);
    % excitatory to inhibitory
    pc_ei = r < dist;
    c_ei(:,i) = sparse(rand(Me, 1) < pc_ei);
end

%% Combine the matrices
Ge = [K * 125 * c_ee, 250 * c_ei] / M;
Gi = [450 * c_ie, 250 * c_ii] / M;
ampa_ratio = 0.5;
nmda_ratio = 1 - ampa_ratio;

%%
figure
spy([Ge; Gi])

%% Simulate the system
vr = [repmat(params_e.mem.c, Me, 1); repmat(params_i.mem.c, Mi, 1)];
ur = [repmat(params_e.mem.d, Me, 1); repmat(params_i.mem.d, Mi, 1)];


progress = 0;
h = waitbar(0,'Computing....');
for i = 1:length(t0)-1
    % update membrane potentials
    tic
    t_1 = t0(i);
    t_2 = t0(i+1);
    for j = 1:M
        [T, Y, TE, YE] = ode23(neur(j).odefun, [t_1, t_2], Y0{j}, opts{j});
        if ~isempty(TE)
            spikes(i,j) = 1;
            Y = YE(end,:)';
            Y(1) = vr(j);
            Y(6) = Y(6) + ur(j);
            neur(j).set_spike(TE);
            [T, Y] = ode23(neur(j).odefun, [TE, t_2], Y, opts{j});
        end
        SOL{j}(i+1,:) = Y(end,:);
    end
    toc
    % propagate effect of spikes
    s = spikes(i,:);
    i_e = s(1:Me) * Ge;
    i_i = s(Me+1:end) * Gi;
    for j = 1:M
        % update for next iteration
        Y0{j} = SOL{j}(i+1,:);
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
        progress = progress + 1;
        waitbar(double(i) / N, h);
    end
end
toc(t_on)

%%
% close all
% figure
% for i = 1:4
%     subplot(2,2,i)
%     x = cellfun(@(s) {s(:,i)}, SOL); x = [x{:}];
%     plot(t0, x);
% end

figure, hold on
for j = 1:M
    flag = spikes(:,j);
    ts = find(flag);
    plot(t0(ts), j * ones(size(ts)),'.k');
end
xlim([t0(1), t0(end)])
ylim([0, M + 1]);

%%
con = struct(...
    'c_ee', c_ee,'c_ei', c_ei,...
    'c_ie', c_ie,'c_ii', c_ii,...
    'Ge', Ge,'Gi', Gi);
out = struct(...
    'spikes', sparse(spikes),...
    'state', vertcat(Y0{:}));
data = struct;
data.connectivity = con;
data.simulation = out;
data.clustering = struct('num_klusters', K,'num_exc', Me,'num_ihn', Mi);
data.params = struct(...
    'syn', syn_params,...
    'exc', params_e,...
    'ihn', params_i);
data.time = t0;
data.inputs = struct(...
    'exc', sparse(spikes_ext_e),...
    'ihn', sparse(spikes_ext_i));

filename = ['results-', date(),'.mat'];
save(filename,'-struct','data');