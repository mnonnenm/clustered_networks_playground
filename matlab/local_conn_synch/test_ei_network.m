clearvars;
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

M = 100;
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


N = 1000;
t0 = 0:2:N;
ts = [];
Y0 = arrayfun(@(~) {[-65 + rand(1) * 10, 0.05*rand(1,4), 100*randn(1)]'}, 1:M);
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

K = 2;
c_e1e1 = sparse(rand(Me/K, Me/K) < 0.60);
c_e1e2 = sparse(rand(Me/K, Me/K) < 0.05);
c_e2e1 = sparse(rand(Me/K, Me/K) < 0.05);
c_e2e2 = sparse(rand(Me/K, Me/K) < 0.60);
c_ee = [
    [c_e1e1, c_e1e2];
    [c_e2e1, c_e2e2];
    ];

c_ei = sprand(Me, Mi, 0.6);
c_ie = rand(Mi, Me);
c_ii = sprand(Mi, Mi, 0.2);

Ge = [250 * c_ee, 125 * c_ei] / M;
Gi = [900 * c_ie, 100 * c_ii] / M;
ampa_ratio = 0.3;
nmda_ratio = 1 - ampa_ratio;

%%
vr = [repmat(params_e.mem.c, Me, 1); repmat(params_i.mem.c, Mi, 1)];
ur = [repmat(params_e.mem.d, Me, 1); repmat(params_i.mem.d, Mi, 1)];
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
        
        if i == 500
            rate_e = 0.01;
            rate_i = 0.01;
        end
        
        if j <= Me
            spike_e = poissrnd(rate_e);
            
            Y0{j}(2) = Y0{j}(2) + 1.2 * spike_e;
            Y0{j}(5) = Y0{j}(5) + 0.3 * spike_e;
            spike_i = poissrnd(rate_i);
            Y0{j}(3) = Y0{j}(3) + 1.5 * spike_i;
            % nmda only on excitatory neurons
            Y0{j}(5) = Y0{j}(5) + nmda_ratio * i_e(j);
        end
        % saturate conductances
        Y0{j}(2:end-2) = min(Y0{j}(2:end-2), [100, 100, 50]);
    end
end

%%
% close all
figure
for i = 1:4
    subplot(2,2,i)
    x = cellfun(@(s) {s(:,i)}, SOL); x = [x{:}];
    plot(t0, x);
end

figure, hold on
for j = 1:M
    flag = spikes(:,j);
    ts = find(flag);
    plot(t0(ts), j * ones(size(ts)),'.k');
end
xlim([t0(1), t0(end)])
ylim([0, M + 1]);
 
set(refline(0, Me / 2 + 0.5),'Color','r','LineStyle',':')
set(refline(0, Me + 0.5),'Color','r','LineStyle',':')