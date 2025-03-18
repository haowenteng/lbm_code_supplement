classdef FilmUniform< handle
    properties
        NY
        NX
        NUM
        NPOP = 9
        N = 50000
        NOUTPUT = 100
        f
        f2
        rho
        ux
        uy
        geometry
        conc_wall = 1.0
        conc_inlet = 0.0
        u0 = 0.05
        omega
        weights = [4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0]
        cx = [0, 1, 0, -1, 0, 1, -1, -1, 1]
        cy = [0, 0, 1, 0, -1, 1, 1, -1, -1]
    end
    
    methods
        function obj = FilmUniform()
            obj.NY = 10 * 160;
            obj.NX = 160;
            obj.NUM = obj.NX * obj.NY;
            obj.omega = 1.0 / (4.0 / 1.8 - (4.0 - 1.0) * 0.5);
            obj.geometry = zeros(1, obj.NUM);
            obj.rho = zeros(1, obj.NUM);
            obj.ux = zeros(1, obj.NUM);
            obj.uy = zeros(1, obj.NUM);
            obj.f = zeros(1, obj.NUM * obj.NPOP);
            obj.f2 = zeros(1, obj.NUM * obj.NPOP);
        end
        
        function initialize_geometry(obj)
            for iY = 1:obj.NY
                for iX = 1:obj.NX
                    counter = (iY - 1) * obj.NX + iX;
                    obj.rho(counter) = 0.0;
                    obj.uy(counter) = obj.u0;
                    obj.ux(counter) = 0.0;
                end
            end
            
            for iX = 2:obj.NX-1
                obj.rho(iX) = obj.conc_inlet;
            end
            
            for iY = 1:obj.NY
                obj.rho((iY - 1) * obj.NX + obj.NX) = obj.conc_wall;
            end
            
            obj.writedensity('conc_initial');
            obj.writevelocityx('ux_initial');
            obj.writevelocityy('uy_initial');
        end
        
        function init(obj)
            for iY = 1:obj.NY
                for iX = 1:obj.NX
                    counter = (iY - 1) * obj.NX + iX;
                    dense_temp = obj.rho(counter);
                    ux_temp = obj.ux(counter);
                    uy_temp = obj.uy(counter);
                    
                    for k = 1:obj.NPOP
                        feq = obj.weights(k) * (dense_temp + 3.0 * dense_temp * (obj.cx(k) * ux_temp + obj.cy(k) * uy_temp) ...
                            + 4.5 * dense_temp * ((obj.cx(k) * obj.cx(k) - 1.0 / 3.0) * ux_temp * ux_temp ...
                            + (obj.cy(k) * obj.cy(k) - 1.0 / 3.0) * uy_temp * uy_temp ...
                            + 2.0 * ux_temp * uy_temp * obj.cx(k) * obj.cy(k)));
                        obj.f((counter-1) * obj.NPOP + k) = feq;
                    end
                end
            end
        end
        
        function collide(obj)
            for iY = 1:obj.NY
                for iX = 1:obj.NX
                    counter = (iY - 1) * obj.NX + iX;
                    obj.rho(counter) = 0.0;
                    
                    offset = (counter-1) * obj.NPOP;
                    sum = 0;
                    for k = 1:obj.NPOP
                        sum = sum + obj.f(offset + k);
                    end
                    
                    obj.rho(counter) = sum;
                    dense_temp = obj.rho(counter);
                    ux_temp = obj.ux(counter);
                    uy_temp = obj.uy(counter);
                    
                    feq = zeros(1, obj.NPOP);
                    for k = 1:obj.NPOP
                        feq(k) = obj.weights(k) * dense_temp * (1.0 + 3.0 * (obj.cx(k) * ux_temp + obj.cy(k) * uy_temp) ...
                            + 4.5 * ((obj.cx(k) * obj.cx(k) - 1.0 / 3.0) * ux_temp * ux_temp ...
                            + 2.0 * obj.cx(k) * obj.cy(k) * ux_temp * uy_temp ...
                            + (obj.cy(k) * obj.cy(k) - 1.0 / 3.0) * uy_temp * uy_temp));
                        obj.f2(offset + k) = obj.f(offset + k) - obj.omega * (obj.f(offset + k) - feq(k));
                    end
                end
            end
        end
        
        function update_bounce_back(obj)
            for iX = 2:obj.NX-1
                offset = (iX-1) * obj.NPOP;
                obj.f2(offset + 3) = -obj.f2((obj.NX + iX) * obj.NPOP + 5) + 2 * obj.weights(3) * obj.conc_inlet;
                obj.f2(offset + 6) = -obj.f2((obj.NX + iX + 1) * obj.NPOP + 8) + 2 * obj.weights(6) * obj.conc_inlet;
                obj.f2(offset + 7) = -obj.f2((obj.NX + iX - 1) * obj.NPOP + 9) + 2 * obj.weights(7) * obj.conc_inlet;
                obj.rho(iX) = obj.conc_inlet;
            end
            
            for iY = 1:obj.NY
                ytop = mod(iY, obj.NY) + 1;
                ybottom = mod(iY - 2, obj.NY) + 1;
                offset = (iY - 1) * obj.NX * obj.NPOP + obj.NX * obj.NPOP;
                
                obj.f2(offset + 4) = -obj.f2((iY - 1) * obj.NX * obj.NPOP + (obj.NX - 1) * obj.NPOP + 2) + 2 * obj.weights(4) * obj.conc_wall;
                obj.f2(offset + 7) = -obj.f2((ytop - 1) * obj.NX * obj.NPOP + (obj.NX - 1) * obj.NPOP + 9) + 2 * obj.weights(7) * obj.conc_wall;
                obj.f2(offset + 8) = -obj.f2((ybottom - 1) * obj.NX * obj.NPOP + (obj.NX - 1) * obj.NPOP + 6) + 2 * obj.weights(8) * obj.conc_wall;
                obj.rho((iY - 1) * obj.NX + obj.NX) = obj.conc_wall;
            end
            
            for iX = 2:obj.NX-1
                offset = (obj.NY - 1) * obj.NX * obj.NPOP + iX * obj.NPOP;
                offset2 = (obj.NY - 2) * obj.NX * obj.NPOP + iX * obj.NPOP;
                for k = 1:obj.NPOP
                    obj.f2(offset + k) = obj.f2(offset2 + k);
                end
            end
            
            for iY = 1:obj.NY
                offset = (iY - 1) * obj.NX * obj.NPOP;
                ytop = mod(iY, obj.NY) + 1;
                ybottom = mod(iY - 2, obj.NY) + 1;
                
                obj.f2(offset + 2) = obj.f2((iY - 1) * obj.NX * obj.NPOP + obj.NPOP + 4);
                obj.f2(offset + 6) = obj.f2((ytop - 1) * obj.NX * obj.NPOP + obj.NPOP + 9);
                obj.f2(offset + 9) = obj.f2((ybottom - 1) * obj.NX * obj.NPOP + obj.NPOP + 7);
                obj.rho((iY - 1) * obj.NX + 1) = obj.rho((iY - 1) * obj.NX + 2);
            end
        end
        
        function stream(obj)
            for iY = 1:obj.NY
                for iX = 1:obj.NX
                    counter = (iY - 1) * obj.NX + iX;
                    for iPop = 1:obj.NPOP
                        iX2 = mod(iX - 1 + obj.cx(iPop), obj.NX) + 1;
                        iY2 = mod(iY - 1 + obj.cy(iPop), obj.NY) + 1;
                        counter2 = (iY2 - 1) * obj.NX + iX2;
                        obj.f((counter2-1) * obj.NPOP + iPop) = obj.f2((counter-1) * obj.NPOP + iPop);
                    end
                end
            end
        end
        
        function writedensity(obj, fname)
            filename = ['film_uniform', fname, '.dat'];
            fid = fopen(filename, 'w');
            for iX = 1:obj.NX
                for iY = 1:obj.NY
                    counter = (iY - 1) * obj.NX + iX;
                    fprintf(fid, '%.10f ', obj.rho(counter));
                end
                fprintf(fid, '\n');
            end
            fclose(fid);
        end
        
        function writevelocityx(obj, fname)
            filename = ['film_uniform', fname, '.dat'];
            fid = fopen(filename, 'w');
            for iX = 1:obj.NX
                for iY = 1:obj.NY
                    counter = (iY - 1) * obj.NX + iX;
                    fprintf(fid, '%.10f ', obj.ux(counter));
                end
                fprintf(fid, '\n');
            end
            fclose(fid);
        end
        
        function writevelocityy(obj, fname)
            filename = ['film_uniform', fname, '.dat'];
            fid = fopen(filename, 'w');
            for iX = 1:obj.NX
                for iY = 1:obj.NY
                    counter = (iY - 1) * obj.NX + iX;
                    fprintf(fid, '%.10f ', obj.uy(counter));
                end
                fprintf(fid, '\n');
            end
            fclose(fid);
        end
        
        function finish_simulation(obj)
            delete(obj.geometry);
            delete(obj.rho);
            delete(obj.ux);
            delete(obj.uy);
            delete(obj.f);
            delete(obj.f2);
        end
        
        function run_simulation(obj)
            obj.initialize_geometry();
            obj.init();
            
            for counter = 0:obj.N
                obj.collide();
                obj.update_bounce_back();
                obj.stream();
                
                if mod(counter, obj.NOUTPUT) == 0
                    disp(['Counter=', num2str(counter)]);
                    filewritedensity = ['film_uniform', num2str(counter, '%07d')];
                    obj.writedensity(filewritedensity);
                end
            end
            
            obj.finish_simulation();
        end
    end
end