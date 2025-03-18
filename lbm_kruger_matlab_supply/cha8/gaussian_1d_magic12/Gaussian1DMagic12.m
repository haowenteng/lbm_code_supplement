% filepath: /d:/BaiduSyncdisk/VScode/Books/lbm_kruger/lbm_principles_practice/code-master/chapter8/gaussian_1d_magic12.m
classdef Gaussian1DMagic12<handle
    properties
        NX = 50
        N = 20
        NOUTPUT = 1
        xInit = 25
        sigma = 3
        g
        g2
        phase
        ux_main = 0.0
        weights = [2.0/3.0, 1.0/6.0, 1.0/6.0]
        cx = [0, 1, -1]
        complement = [1, 3, 2]
        magic = 1.0 / 12.0
        omega = 0.5
    end
    
    methods
        function obj = Gaussian1DMagic12()
            obj.g = zeros(obj.NX, 3);
            obj.g2 = zeros(obj.NX, 3);
            obj.phase = zeros(obj.NX, 1);
            obj.init();
        end
        
        function init(obj)
            for iX = 1:obj.NX
                if abs(iX -1- obj.xInit) <= 3 * obj.sigma
                    obj.phase(iX) = exp(-0.5 * (iX-1 - obj.xInit)^2 / obj.sigma^2);
                else
                    obj.phase(iX) = 0.0;
                end
            end
            
            for iX = 1:obj.NX
                for k = 1:3
                    obj.g(iX, k) = obj.weights(k) * obj.phase(iX) * (1.0 + 3.0 * obj.cx(k) * obj.ux_main + 4.5 * (obj.cx(k)^2 - 1.0 / 3.0) * obj.ux_main^2);
                end
            end
        end
        
        function writephase(obj, fName)
            fullName = ['gaussian_1d_magic12', fName, '.dat'];
            fid = fopen(fullName, 'w');
            for iX = 1:obj.NX
                fprintf(fid, '%.10f ', obj.phase(iX));
            end
            fprintf(fid, '\n');
            fclose(fid);
        end
        
        function collide_bulk(obj)
            for iX = 1:obj.NX
                obj.phase(iX) = 0.0;
                for iPop = 1:3
                    obj.phase(iX) = obj.phase(iX) + obj.g(iX, iPop);
                end
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
                    geq(k) = obj.weights(k) * obj.phase(iX) * (1.0 + 3.0 * obj.cx(k) * obj.ux_main + 4.5 * (obj.cx(k)^2 - 1.0 / 3.0) * obj.ux_main^2);
                end
                
                for k = 1:3
                    geq_plus(k) = 0.5 * (geq(k) + geq(obj.complement(k)));
                    geq_minus(k) = 0.5 * (geq(k) - geq(obj.complement(k)));
                end
                
                omega_minus_phase = obj.omega;
                omega_plus_phase = 1.0 / (obj.magic / (1.0 / obj.omega - 0.5) + 0.5);
                
                for k = 1:3
                    obj.g2(iX, k) = obj.g(iX, k) - omega_plus_phase * (g_plus(k) - geq_plus(k)) - omega_minus_phase * (g_minus(k) - geq_minus(k));
                end
            end
        end
        
        function stream(obj)
            for iX = 1:obj.NX
                for iPop = 1:3
                    iX2 = mod(iX - obj.cx(iPop) - 1, obj.NX) + 1;
                    obj.g(iX, iPop) = obj.g2(iX2, iPop);
                end
            end
        end
        
        function run_simulation(obj)
            obj.init();    
            for counter = 0:obj.N
                obj.collide_bulk();
                obj.stream();
                
                if mod(counter, obj.NOUTPUT) == 0
                    disp(['Counter=', num2str(counter)]);
                    filewritephase = sprintf('gaussian_1d_magic12%06d', counter);
                    obj.writephase(filewritephase);
                end
            end
        end
    end
end