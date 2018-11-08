classdef RobotObj < handle
    %wrapper of robot object
    % under development
    % Yu Zhao, Aug04, 2017
    
    % public properties
    properties
        ModelName; % model name of of robot
    end
    % hidden but public properties
    properties (Hidden)
        r;  % struct contains robot parameters
        rtb; % robotic toolbox object
        trans; % transmission matrix from motor angle to joint angle
        trans_rtb; % transmission matrix in robotic toolbox convention
    end
    
    methods
        % constructor
        function obj = RobotObj(model)
            switch model
                case 'LRMate'
                    obj = obj.buildLRMate;
                otherwise
                    error('Robot model not supported.');
            end
            obj.ModelName=model;
        end
        
        % build LRMate model, only kinematic, self edited parameters
        function obj = buildLRMate(obj)
            r.name = 'LRMate';
            r.manufacturer = 'FANUC Corporation';
            r.n = 6;

            % DH parameters
            r.Theta0  = [0,-pi/2,0,0,0,0];                                       % Joint variable offset (rad)
            r.D       = [0,0,0,-0.42,0,-0.08];                              % Joint extension (m)
            r.A       = [0.05,0.44,0.035,0,0,0];                              % Joint offset (m)
            r.Alpha   = [-pi/2,pi,-pi/2,pi/2,-pi/2,pi];                           % Joint twist (rad)
            r.Qlim    = [-185 185; -60 76; -132.04 90; ...
                        -210 210; -120 120; -205 205]'*pi/180;                  % Joint angle limit (rad)
            
            % some common poses
            r.qz  = [0 0 0 0 0 0];         % zero angles, L shaped pose
            r.qh  = [0 0 0 0 -pi/2 0];     % home pose, J5 down
                        
            obj.r = r;
            
            % define each link
            for i = 1:6
                L(i) = Revolute('d', r.D(i), 'a', r.A(i), 'alpha', r.Alpha(i), ...
                    'offset', r.Theta0(i));    
                    % the Coulomb friction is disabled as multiplied by 0. The friction
                    % parameters and modeling seem not quite fit with the experimental
                    % data and thus is disabled here.
            end  

            % construct full robot 
            obj. rtb = SerialLink(L, 'name', r.name, 'manufacturer', r.manufacturer);
            obj.rtb.tool = eye(4);
            obj.rtb.tool(3,4) = 0.027+0.0157+0.027 + 0.1714;%0.14605;
            
            obj.trans = eye(6);
            
            obj.trans_rtb = eye(6);
            % preferred joint range
            %obj.r.pos_lb = [-170, -100, -72, -190, -125, -360] + 2; %[deg]     
            %obj.r.pos_ub = [170,  145,  240, 190, 125, 360] - 2; %[deg]
            obj.r.pos_lb = [-60, -40, -10, -30, -60, -80] + 2; %[deg]     
            obj.r.pos_ub = [60,  30,  30, 30, 60, 80] - 2; %[deg]
        end
        
        
    end
    
end

