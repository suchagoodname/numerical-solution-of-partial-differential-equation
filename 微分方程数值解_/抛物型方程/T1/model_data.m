function pde = model_data()
pde = struct('u_initial',@u_initial,'u_left',@u_left,...
    'u_right',@u_right,'f',@f,'time_grid',@time_grid,...
    'space_grid',@space_grid,'a',@a);

    function [T,tau] = time_grid(NT)
        T = linspace(0,1,NT+1);
        tau = 1/NT;
    end
    function [X,h] = space_grid(NS)
        X = linspace(0,1,NS+1)';
        h = 1/NS;
    end
    function u = u_initial(x)
        u = 2*sin(2*pi*x);
    end
    function u = u_left(t)
        u = zeros(size(t));
    end
    function u = u_right(t)
        u = zeros(size(t));
    end
    function f = f(x,t)
        f = zeros(size(x));
    end
    function a = a()
        a = 1/16;
    end
end

