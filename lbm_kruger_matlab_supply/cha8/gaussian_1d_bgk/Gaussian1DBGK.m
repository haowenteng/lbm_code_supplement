classdef Gaussian1DBGK< handle
    properties
        NX = 50; % Domain size
        N = 20; % Time steps
        NOUTPUT = 1; % Output frequency
        xInit = 25; % Initial position
        sigma = 3; % Standard deviation
        g; % Distribution function
        g2; % Post-collision distribution function
        phase; % Phase field
        ux_main = 0.5; % Diffusion parameter
        weights = [2.0/3.0, 1.0/6.0, 1.0/6.0]; % Weights
        cx = [0, 1, -1]; % Lattice velocities
        complement = [1, 3, 2]; % Complement of the lattice velocities  速度对称方向
        omega = 0.3; % Relaxation parameter
    end
    
    methods
        function obj = Gaussian1DBGK()
            obj.g = zeros(obj.NX, 3);
            obj.g2 = zeros(obj.NX, 3);
            obj.phase = zeros(obj.NX, 1);
        end
        
        function init(obj)
            % Phase initialization
            for iX = 1:obj.NX
                if abs(iX -1- obj.xInit) <= 3 * obj.sigma
                    obj.phase(iX) = exp(-0.5 * (iX-1 - obj.xInit)^2 / obj.sigma^2);
                else
                    obj.phase(iX) = 0.0;
                end
            end
            
            % Bulk nodes initialization
            for iX = 1:obj.NX
                for k = 1:3
                    obj.g(iX, k) = obj.weights(k) * obj.phase(iX) * ...
                        (1.0 + 3.0 * obj.cx(k) * obj.ux_main + ...
                        4.5 * (obj.cx(k)^2 - 1.0/3.0) * obj.ux_main^2);
                end
            end
        end
        
        function collide_bulk(obj)
            % Construction of phase field
            for iX = 1:obj.NX
                obj.phase(iX) = sum(obj.g(iX, :));
            end
            
            for iX = 1:obj.NX
                g_plus = zeros(1, 3);
                g_minus = zeros(1, 3);
                geq = zeros(1, 3);
                geq_plus = zeros(1, 3);
                geq_minus = zeros(1, 3);
                
                for k = 1:3
                    g_plus(k) = 0.5 * (obj.g(iX, k) + obj.g(iX, obj.complement(k)));
                    g_minus(k) = 0.5 * (obj.g(iX, k) - obj.g(iX, obj.complement(k)));
                end
                
                for k = 1:3
                    geq(k) = obj.weights(k) * obj.phase(iX) * ...
                        (1.0 + 3.0 * obj.cx(k) * obj.ux_main + ...
                        4.5 * (obj.cx(k)^2 - 1.0/3.0) * obj.ux_main^2);
                end
                
                for k = 1:3
                    geq_plus(k) = 0.5 * (geq(k) + geq(obj.complement(k)));
                    geq_minus(k) = 0.5 * (geq(k) - geq(obj.complement(k)));
                end
                
                omega_minus_phase = obj.omega;
                omega_plus_phase = obj.omega;
                
                for k = 1:3
                    obj.g2(iX, k) = obj.g(iX, k) - omega_plus_phase * (g_plus(k) - geq_plus(k)) - ...
                        omega_minus_phase * (g_minus(k) - geq_minus(k));
                end
            end
        end
        
       function stream(obj)
           for iX = 1:obj.NX
               for iPop = 1:3
                   iX2 = mod(iX - 1 - obj.cx(iPop), obj.NX) + 1;
                    obj.g(iX, iPop) = obj.g2(iX2, iPop);
                end
            end
        end
        
        function writephase(obj, fName)
            fullName = ['gaussian_1d_bgk', fName, '.dat'];
            fout = fopen(fullName, 'w');
            for iX = 1:obj.NX
                fprintf(fout, '%.10f ', obj.phase(iX));
            end
            fprintf(fout, '\n');
            fclose(fout);
        end
        
        function run_simulation(obj)
            obj.init();
            for counter = 0:obj.N
                obj.collide_bulk();
                obj.stream();
                
                if mod(counter, obj.NOUTPUT) == 0
                    disp(['Counter=', num2str(counter)]);
                    filewritephase = sprintf('gaussian_1d_bgk%06d', counter);
                    obj.writephase(filewritephase);
                end
            end
        end
    end
end