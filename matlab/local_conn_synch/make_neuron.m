function handles = make_neuron(params)
    %MAKE_ODE Summary of this function goes here
    %   Detailed explanation goes here
    this = struct('tspike', -inf);
    handles = struct('odefun', @f,'events', @spike,...
        'params',@p,'set_spike',@set_t);
    
    function value = p()
        value = params;
    end
    
    function set_t(tspike)
        this.tspike = tspike;
    end
    
    function y1 = f(t, y)
        V = y(1);
        g_ampa = y(2);
        g_gaba = y(3);
        g_nmda = y(4);
        r_nmda = y(5);
        
        y1(1) = params.mem.k * (V - params.mem.Er) *(V - params.mem.Et) +...
            g_ampa * (params.syn.Eampa - V) +...
            g_gaba * (params.syn.Egaba - V) +...
            g_nmda * (params.syn.Enmda - V) * B(V);
        y1(2) = -g_ampa / params.syn.Tampa;
        y1(3) = -g_gaba / params.syn.Tgaba;
        y1(4) = -g_nmda / params.syn.Tnmda1 + r_nmda / params.syn.Tnmda2;
        y1(5) = -r_nmda / params.syn.Tnmda2;
        
        if length(y) == 6
            y1(6) = params.mem.a * (params.mem.b * (y(1) - params.mem.Er) - y(6));
            y1(1) = y1(1) - y(6);
        end
        
        y1(1) = y1(1) / params.mem.C;
        y1 = y1(:);
    end
    
    function x = B(V)
        x = 1./(1 + exp(-0.062 * V) / 3.57);
    end
    
    function [value, isterminal, direction] = spike(~, y)
        value = y(1) + 20;
        isterminal = 1;
        direction = 1;
    end
end

