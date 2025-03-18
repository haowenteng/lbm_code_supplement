% filepath: /d:/BaiduSyncdisk/VScode/Books/lbm_kruger/lbm_principles_practice/code-master/chapter8/gaussian_2d_trt.m
classdef Gaussian2DTRT<handle
    properties
        NX = 512
        NY = 512
        N = 200
        NOUTPUT = 10
        NPOP = 9
        xInit = 200.0
        yInit = 200.0
        sigma = 10.0
        g
        g2
        phase
        ux_main = 0.1
        uy_main = 0.1
        lambda = 1.0 / 12.0
        omegaminus = 1.95
        omegaplus
        weights = [4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0]
        cx = [0, 1, 0, -1, 0, 1, -1, -1, 1]
        cy = [0, 0, 1, 0, -1, 1, 1, -1, -1]
        complement = [  1 , 4  ,   5   ,  2  ,   3   ,  8 ,    9 ,    6  ,   7]
    end
    
    methods
        function obj = Gaussian2DTRT()
            obj.g = zeros(obj.NY, obj.NX, obj.NPOP);
            obj.g2 = zeros(obj.NY, obj.NX, obj.NPOP);
            obj.phase = zeros(obj.NY, obj.NX);
            obj.omegaplus = 1.0 / (0.5 + obj.lambda * obj.omegaminus / (1.0 - 0.5 * obj.omegaminus));
            obj.init();
        end
        
        function init(obj)
            disp(['Omegaplus: ', num2str(obj.omegaplus), ' Omegaminus: ', num2str(obj.omegaminus)]);
            for iY = 1:obj.NY
                for iX = 1:obj.NX
                    radius = sqrt((iX - obj.xInit)^2 + (iY - obj.yInit)^2);
                    if radius <= 3.0 * obj.sigma
                        obj.phase(iY, iX) = exp(-0.5 * radius^2 / obj.sigma^2);
                    else
                        obj.phase(iY, iX) = 0.0;
                    end
                end
            end
            
            for iY = 1:obj.NY
                for iX = 1:obj.NX
                    for iPOP = 1:obj.NPOP
                        obj.g(iY, iX, iPOP) = obj.weights(iPOP) * obj.phase(iY, iX) * ...
                            (1.0 + 3.0 * obj.cx(iPOP) * obj.ux_main + 3.0 * obj.cy(iPOP) * obj.uy_main + ...
                            4.5 * (obj.cx(iPOP)^2 - 1.0 / 3.0) * obj.ux_main^2 + ...
                            9.0 * obj.cx(iPOP) * obj.cy(iPOP) * obj.ux_main * obj.uy_main + ...
                            4.5 * (obj.cy(iPOP)^2 - 1.0 / 3.0) * obj.uy_main^2);
                    end
                end
            end
        end
        
        function writephase(obj, fName)
            fullName = ['gaussian_2d_trt', fName, '.dat'];
            fid = fopen(fullName, 'w');
            for iY = 1:obj.NY
                for iX = 1:obj.NX
                    fprintf(fid, '%.10f ', obj.phase(iY, iX));
                end
                fprintf(fid, '\n');
            end
            fclose(fid);
        end
        
        function collide_bulk(obj)
            for iY = 1:obj.NY
                for iX = 1:obj.NX
                    obj.phase(iY, iX) = 0.0;
                    for iPop = 1:obj.NPOP
                        obj.phase(iY, iX) = obj.phase(iY, iX) + obj.g(iY, iX, iPop);
                    end
                end
            end
            
            for iY = 1:obj.NY
                for iX = 1:obj.NX
                    gplus = zeros(1, obj.NPOP);
                    gminus = zeros(1, obj.NPOP);
                    geq = zeros(1, obj.NPOP);
                    geqplus = zeros(1, obj.NPOP);
                    geqminus = zeros(1, obj.NPOP);
                    
                    for iPOP = 1:obj.NPOP
                        gplus(iPOP) = 0.5 * (obj.g(iY, iX, iPOP) + obj.g(iY, iX, obj.complement(iPOP)));
                        gminus(iPOP) = 0.5 * (obj.g(iY, iX, iPOP) - obj.g(iY, iX, obj.complement(iPOP)));
                        
                        geq(iPOP) = obj.weights(iPOP) * obj.phase(iY, iX) * ...
                            (1.0 + 3.0 * obj.cx(iPOP) * obj.ux_main + 3.0 * obj.cy(iPOP) * obj.uy_main + ...
                            4.5 * (obj.cx(iPOP)^2 - 1.0 / 3.0) * obj.ux_main^2 + ...
                            9.0 * obj.cx(iPOP) * obj.cy(iPOP) * obj.ux_main * obj.uy_main + ...
                            4.5 * (obj.cy(iPOP)^2 - 1.0 / 3.0) * obj.uy_main^2);
                    end
                    
                    for iPOP = 1:obj.NPOP
                        geqplus(iPOP) = 0.5 * (geq(iPOP) + geq(obj.complement(iPOP)));
                        geqminus(iPOP) = 0.5 * (geq(iPOP) - geq(obj.complement(iPOP)));
                        obj.g2(iY, iX, iPOP) = obj.g(iY, iX, iPOP) - obj.omegaplus * (gplus(iPOP) - geqplus(iPOP)) ...
                            - obj.omegaminus * (gminus(iPOP) - geqminus(iPOP));
                    end
                end
            end
        end
        
        function stream(obj)
            for iY = 1:obj.NY
                for iX = 1:obj.NX
                    for iPOP = 1:obj.NPOP
                        iX2 = mod(iX - obj.cx(iPOP) - 1, obj.NX) + 1;
                        iY2 = mod(iY - obj.cy(iPOP) - 1, obj.NY) + 1;
                        obj.g(iY, iX, iPOP) = obj.g2(iY2, iX2, iPOP);
                    end
                end
            end
        end
        
        function run_simulation(obj)
            for counter = 0:obj.N
                obj.collide_bulk();
                obj.stream();
                
                if mod(counter, obj.NOUTPUT) == 0
                    disp(counter);
                    filewritephase = sprintf('gaussian_2d_trt%06d', counter);
                    obj.writephase(filewritephase);
                end
            end
        end
    end
end