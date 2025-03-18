% Diffusion from Cylinder Without Flow
% D2Q9 Model
% TRT collision model
% two different boundary conditions: anti-bounce-back and Inamuro’s boundary condition
classdef cylinder2   
    properties
        NY  % number of lattice nodes in the y direction
        NX  % number of lattice nodes in the x direction
        NUM % total number of lattice nodes
        radius  % radius of the cylinder
        NPOP = 9    % number of discrete velocities
        N = 64000   % number of iterations
        NOUTPUT = 200   % output frequency
        f   % distribution function
        f2  % distribution function after collision
        rho % density
        ux  % velocity in the x direction
        uy  % velocity in the y direction
        geometry    % geometry of the domain
        conc_wall = 1.0   % concentration of the wall
        bb_nodes    % bounce-back nodes
        dirs    % directions for bounce-back nodes
        omega = 1.0 / (0.5 + 0.5 / 32.0)    % relaxation parameter
        omega_plus  % relaxation parameter for the positive moments
        omega_minus % relaxation parameter for the negative moments
        ce = 1.0 / 3.0  % speed of sound squared
        weights = [4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0]   % weights for the discrete velocities
        weights_trt = [0.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/12.0, 1.0/12.0, 1.0/12.0, 1.0/12.0]   % weights for the discrete velocities in the TRT model
        cx = [0, 1, 0, -1, 0, 1, -1, -1, 1] % x component of the discrete velocities
        cy = [0, 0, 1, 0, -1, 1, 1, -1, -1] % y component of the discrete velocities
        complement = [0, 3, 4, 1, 2, 7, 8, 5, 6]    % complement of the discrete velocities
    end
     
    methods
        function obj = cylinder2()   % constructor
            obj.omega_plus = 2.0 - obj.omega;       % relaxation parameter for the positive moments
            obj.omega_minus = obj.omega;        % relaxation parameter for the negative moments
        end
        
        function obj = initialize_geometry(obj)  % initialize the geometry
            obj.NY = 129;   % number of lattice nodes in the y direction
            obj.NX = 129;   % number of lattice nodes in the x direction
            obj.NUM = obj.NX * obj.NY;  % total number of lattice nodes
            obj.geometry = -ones(obj.NUM, 1);   % geometry of the domain
            obj.rho = zeros(obj.NUM, 1);
            obj.ux = zeros(obj.NUM, 1);
            obj.uy = zeros(obj.NUM, 1);
            obj.radius = 40;
            
            for iX = 1:obj.NX
                for iY = 1:obj.NY
                    counter = (iX - 1) * obj.NY + iY;
                    if (iX - (obj.NX - 1) / 2)^2 + (iY - (obj.NY - 1) / 2)^2 <= obj.radius^2
                        obj.geometry(counter) = 1;
                    end
                end
            end
            obj.writegeometry('geometry_before');
            
            for iX = 1:obj.NX
                for iY = 1:obj.NY
                    counter = (iX - 1) * obj.NY + iY;
                    if obj.geometry(counter) == 1
                        flag = false;
                        for k = 1:obj.NPOP
                            counter2 = (iX + obj.cx(k)) * obj.NY + iY + obj.cy(k);
                            if obj.geometry(counter2) == -1
                                flag = true;
                            end
                        end
                        if flag
                            obj.geometry(counter) = 0;
                        end
                    end
                end
            end
            
            for iX = 1:obj.NX
                for iY = 1:obj.NY
                    counter = (iX - 1) * obj.NY + iY;
                    if obj.geometry(counter) == 0
                        flag = false;
                        for k = 1:obj.NPOP
                            counter2 = (iX + obj.cx(k)) * obj.NY + iY + obj.cy(k);
                            if obj.geometry(counter2) == 1
                                flag = true;
                            end
                        end
                        if ~flag
                            obj.geometry(counter) = -1;
                        end
                    end
                end
            end
            
            obj.writegeometry('geometry_after');
            
            for counter = 1:obj.NUM
                if obj.geometry(counter) == 0
                    obj.rho(counter) = obj.conc_wall;
                    obj.ux(counter) = 0.0;
                    obj.uy(counter) = 0.0;
                    obj.bb_nodes = [obj.bb_nodes; counter];
                elseif obj.geometry(counter) == -1
                    obj.rho(counter) = -1.0;
                    obj.ux(counter) = 0.0;
                    obj.uy(counter) = 0.0;
                else
                    obj.rho(counter) = 0.0;
                    obj.ux(counter) = 0.0;
                    obj.uy(counter) = 0.0;
                end
            end
            
            obj.dirs = cell(length(obj.bb_nodes), 1);
            for counter = 1:length(obj.bb_nodes)
                for k = 2:obj.NPOP
                    counter2 = obj.bb_nodes(counter) + obj.cy(k) * obj.NX + obj.cx(k);
                    if obj.geometry(counter2) == 1
                        obj.dirs{counter} = [obj.dirs{counter}, k];
                    end
                end
            end
        end
        
        function obj = init(obj)    % initialize the distribution function
            obj.f = zeros(obj.NUM, obj.NPOP);
            obj.f2 = zeros(obj.NUM, obj.NPOP);
            
            for iY = 1:obj.NY
                for iX = 1:obj.NX
                    counter = (iX - 1) * obj.NY + iY;
                    dense_temp = obj.rho(counter);
                    ux_temp = obj.ux(counter);
                    uy_temp = obj.uy(counter);
                    
                    for k = 1:obj.NPOP  % loop over discrete velocities
                        feq = obj.weights(k) * (dense_temp + 3.0 * dense_temp * (obj.cx(k) * ux_temp + obj.cy(k) * uy_temp) ...
                            + 4.5 * dense_temp * ((obj.cx(k)^2 - 1.0 / 3.0) * ux_temp^2 ...
                            + (obj.cy(k)^2 - 1.0 / 3.0) * uy_temp^2 + 2.0 * ux_temp * uy_temp * obj.cx(k) * obj.cy(k)));
                        obj.f(counter, k) = feq;
                    end
                end
            end
        end
        
        function obj = collide(obj) % collision step
            for counter = 1:obj.NUM
                if obj.geometry(counter) == 1
                    obj.rho(counter) = sum(obj.f(counter, :));
                    dense_temp = obj.rho(counter);
                    ux_temp = obj.ux(counter);
                    uy_temp = obj.uy(counter);
                    
                    feq_minus = zeros(1, obj.NPOP);
                    feq_plus = zeros(1, obj.NPOP);
                    f_minus = zeros(1, obj.NPOP);
                    f_plus = zeros(1, obj.NPOP);
                    
                    f_plus(1) = obj.f(counter, 1);
                    f_plus(2) = 0.5 * (obj.f(counter, 2) + obj.f(counter, 4));
                    f_plus(3) = 0.5 * (obj.f(counter, 3) + obj.f(counter, 5));
                    f_plus(4) = f_plus(2);
                    f_plus(5) = f_plus(3);
                    f_plus(6) = 0.5 * (obj.f(counter, 6) + obj.f(counter, 8));
                    f_plus(7) = 0.5 * (obj.f(counter, 7) + obj.f(counter, 9));
                    f_plus(8) = f_plus(6);
                    f_plus(9) = f_plus(7);
                    
                    f_minus(1) = 0.0;
                    f_minus(2) = 0.5 * (obj.f(counter, 2) - obj.f(counter, 4));
                    f_minus(3) = 0.5 * (obj.f(counter, 3) - obj.f(counter, 5));
                    f_minus(4) = -f_minus(2);
                    f_minus(5) = -f_minus(3);
                    f_minus(6) = 0.5 * (obj.f(counter, 6) - obj.f(counter, 8));
                    f_minus(7) = 0.5 * (obj.f(counter, 7) - obj.f(counter, 9));
                    f_minus(8) = -f_minus(6);
                    f_minus(9) = -f_minus(7);
                    
                    feq_minus(1) = 0.0;
                    feq_minus(2) = obj.weights_trt(2) * dense_temp * (obj.cx(2) * ux_temp + obj.cy(2) * uy_temp);
                    feq_minus(3) = obj.weights_trt(3) * dense_temp * (obj.cx(3) * ux_temp + obj.cy(3) * uy_temp);
                    feq_minus(4) = -feq_minus(2);
                    feq_minus(5) = -feq_minus(3);
                    feq_minus(6) = obj.weights_trt(6) * dense_temp * (obj.cx(6) * ux_temp + obj.cy(6) * uy_temp);
                    feq_minus(7) = obj.weights_trt(7) * dense_temp * (obj.cx(7) * ux_temp + obj.cy(7) * uy_temp);
                    feq_minus(8) = -feq_minus(6);
                    feq_minus(9) = -feq_minus(7);
                    
                    u_sq = ux_temp^2 + uy_temp^2;
                    
                    feq_plus(2) = obj.weights_trt(2) * dense_temp * (obj.ce + 0.5 * (3.0 * (obj.cx(2) * ux_temp + obj.cy(2) * uy_temp)^2 - u_sq));
                    feq_plus(3) = obj.weights_trt(3) * dense_temp * (obj.ce + 0.5 * (3.0 * (obj.cx(3) * ux_temp + obj.cy(3) * uy_temp)^2 - u_sq));
                    feq_plus(4) = feq_plus(2);
                    feq_plus(5) = feq_plus(3);
                    feq_plus(6) = obj.weights_trt(6) * dense_temp * (obj.ce + 0.5 * (3.0 * (obj.cx(6) * ux_temp + obj.cy(6) * uy_temp)^2 - u_sq));
                    feq_plus(7) = obj.weights_trt(7) * dense_temp * (obj.ce + 0.5 * (3.0 * (obj.cx(7) * ux_temp + obj.cy(7) * uy_temp)^2 - u_sq));
                    feq_plus(8) = feq_plus(6);
                    feq_plus(9) = feq_plus(7);
                    feq_plus(1) = dense_temp - 2.0 * (feq_plus(2) + feq_plus(3) + feq_plus(6) + feq_plus(7));
                    
                    for k = 1:obj.NPOP
                        obj.f2(counter, k) = obj.f(counter, k) - obj.omega_plus * (f_plus(k) - feq_plus(k)) - obj.omega_minus * (f_minus(k) - feq_minus(k));
                    end
                end
            end
        end
        
        function obj = update_bounce_back(obj)  % bounce-back nodes
            for counter = 1:length(obj.bb_nodes)
                for k = 1:length(obj.dirs{counter})
                    dir = obj.dirs{counter}(k);
                    counter2 = obj.bb_nodes(counter) + obj.cy(dir) * obj.NX + obj.cx(dir);
                    obj.f2(obj.bb_nodes(counter), dir) = -obj.f2(counter2, obj.complement(dir)) + 2 * obj.weights(dir) * obj.conc_wall;
                end
            end
        end
        
        function obj = stream(obj)  % streaming step
            for counter = 1:obj.NUM
                if obj.geometry(counter) == 1
                    for iPop = 1:obj.NPOP
                        counter2 = counter - obj.cy(iPop) * obj.NX - obj.cx(iPop);
                        obj.f(counter, iPop) = obj.f2(counter2, iPop);
                    end
                end
            end
        end
        
        function writedensity(obj, fname)   % write the density
            filename = ['cylinder_', fname, '.dat'];
            fout = fopen(filename, 'w');
            for iX = 1:obj.NX
                for iY = 1:obj.NY
                    counter = (iX - 1) * obj.NY + iY;
                    fprintf(fout, '%.10f ', obj.rho(counter));
                end
                fprintf(fout, '\n');
            end
            fclose(fout);
        end
        
        function writegeometry(obj, fname)    % write the geometry
            filename = ['cylinder_', fname, '.dat'];
            fout = fopen(filename, 'w');
            for iX = 1:obj.NX
                for iY = 1:obj.NY
                    counter = (iX - 1) * obj.NY + iY;
                    fprintf(fout, '%.10f ', obj.geometry(counter));
                end
                fprintf(fout, '\n');
            end
            fclose(fout);
        end
        
        function finish_simulation(obj)   % finish the simulation
            clear obj.geometry obj.rho obj.ux obj.uy obj.f obj.f2 obj.dirs;
        end
        
        function run_simulation(obj)    % run the simulation
            obj = obj.initialize_geometry();
            obj = obj.init();
            
            for counter = 1:obj.N
                obj = obj.collide();
                obj = obj.update_bounce_back();
                obj = obj.stream();
                
                if mod(counter, obj.NOUTPUT) == 0
                    disp(['Counter=', num2str(counter)]);
                    filewritedensity = sprintf('cylinder%07d', counter);
                    obj.writedensity(filewritedensity);
                end
            end
            
            obj.finish_simulation();
        end
    end
end