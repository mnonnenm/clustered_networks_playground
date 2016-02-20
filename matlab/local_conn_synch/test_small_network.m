clearvars;
syn_params = struct(...
    'Tampa', 3.5,'Eampa', 0,...
    'Tnmda1', 80,'Tnmda2', 20,'Enmda', 0,...
    'Tgaba', 5.5,'Egaba', -65);
params = struct(...
    'mem', struct('Er',-65,'Et',-40,'tau', 20,'Vr',-70,'Tref', 3),...
    'syn', syn_params);
rate_e = 0.40;
rate_i = 0.40;

M = 20;
neur = cell(M, 1);
for i = 1:M
    if i < 0.8 * M
        neur{i} = make_neuron(params_e);
    else
        neur{i} = make
end
neur = [neur{:}];
N = 1000;
t0 = 0:N;
ts = [];
Y0 = arrayfun(@(~) {[-50, 0.05*rand(1,4)]'}, 1:M);
SOL = cellfun(@(y0) {zeros(N+1, length(y0))}, Y0);
for i = 1:M
    SOL{i}(1,:) = Y0{i};
end
for i = 1:M
    opts{i} = odeset(...
        'Events',neur(i).events,...
        'RelTol',1e-3,'AbsTol',1e-3,...
        'InitialStep', 1,'MaxStep', 1);
end
exit = false;
Tend = 1;
spikes = zeros(N+1, M);

c_ee = sprand(M, M, 0.3);
g_ee = spfun(@(~) 0.05 * rand(1), c_ee);

vr = repmat(params.mem.Vr, M, 1);
for i = 1:N
    % update membrane potentials
    tic
    for j = 1:M
        [T, Y, TE, YE] = ode45(neur(j).odefun, [t0(i), t0(i+1)], Y0{j}, opts{j});
        if ~isempty(TE)
            spikes(i,j) = 1;
            Y = YE(end,:)';
            Y(1) = vr(j);
            neur(j).set_spike(TE);
            [T, Y] = ode15s(neur(j).odefun, [TE, t0(i+1)], Y, opts{j});
        end
        SOL{j}(i+1,:) = Y(end,:);
    end
    toc
    % propagate effect of spikes
    
    s = spikes(i,:);
    i_ee = s * g_ee;
    for j = 1:M
        % update for next iteration
        Y0{j} = SOL{j}(i+1,:);
        spike_e = poissrnd(rate_e);
        Y0{j}(2) = Y0{j}(2) + 0.020 * spike_e + 0.75 * i_ee(j);
        Y0{j}(5) = Y0{j}(5) + 0.005 * spike_e + 0.25 * i_ee(j);
        spike_i = poissrnd(rate_i);
        Y0{j}(3) = Y0{j}(3) + 0.03*spike_i;
        
        Y0{j}(2:end) = min(Y0{j}(2:end), repmat(0.2, 1, length(Y0{j})-1));
    end
end

%%
close all
V = cellfun(@(s) {s(:,1)}, SOL); V = [V{:}];
plot(V);

figure, hold on
for j = 1:M
    flag = spikes(:,j);
    ts = find(flag);
    plot(ts, j * ones(size(ts)),'.k');
end
xlim([t0(1), t0(end)])
ylim([0, M + 1]);