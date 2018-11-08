classdef RobVizLRMate < RobotObj
    %ROBVIZ visualization class for robot
    % @Yu Zhao
    
    properties (Hidden)
        si;
        parts;
        h_rob;
        h_marker;
        sta;
        xyz_target;
        jnt_pos;
        jnt_vel;
        Vm;
        Am;
        dt;
        delta;
        bound_c=[0.2595;0.6345;0.330];% center of boundary ball
        bound_rout=0.95;% outer radius of boundary ball
        bound_rin=0.35;% inner radius of boundary ball
        scale=0.05;
        v_threshold=1/100;% velocity threshold, lower cmd suppressed
        T_B2W;% base transfrom from world frame
        x_bound=[30,80]/100;
        y_bound=[-20,40]/100;
        z_bound=[-0,50]/100;
        enb = false;
        udp_enable=false;
        home_xyz;
    end
    
    methods
        function obj = RobVizLRMate()
            obj = obj@RobotObj('LRMate');
            obj.loadCAD;
            obj.si.DH(:,1) = obj.r.Alpha(:);
            obj.si.DH(:,2) = obj.r.A(:);
            obj.si.DH(:,3) = obj.r.D(:);
            obj.si.DH(:,4) = obj.r.Theta0(:);
            obj.h_rob=[];
            obj.h_marker=[];
            obj.sta = true;
            obj.jnt_pos=obj.r.qh(:);% initial position
            obj.jnt_vel=[0 0 0 0 0 0].';
            
            BaseScrewPos = [0, 0, 0; 165, 0, 0; 165, -165, 0; 0, -165, 0]/1000;   % position coordinates of robot base screws in design space frame
            BaseCenter = mean(BaseScrewPos, 1);                             % center position of robot base
            Base2J1O = [0, 0, 0.330]';
            J1Origin = BaseCenter' + Base2J1O;        % robot J1 origin in design space frame
            RPY_DS2B = [0, 0, 0];                                 % rotation angle RPY from design space frame to robot base (J1) frame
            T_B2DS = [[rpy2r(RPY_DS2B, 'zyx') J1Origin]; 0 0 0 1];   % transformation from robot base (J1) frame to design space
            RPY_W2DS = [0, 0, 0];
            WOrigin = [ -0.127, -0.717, 0];      % [-5, -21.25-7, 0] inch define world frame origin as the center of the middle stage
            T_DS2W = [rpy2r(RPY_W2DS, 'zyx') * [eye(3) -WOrigin']; 0 0 0 1];   % transformation from the design space frame to world frame
            obj.T_B2W = T_DS2W * T_B2DS;                     % transformation from robot base (J1) frame to world frame

            [jntT]=kfwd_rob_LRMate_full(obj.jnt_pos,obj.si.DH,obj.rtb.tool,obj.T_B2W);
            obj.xyz_target=jntT(1:3,4,end);
            obj.home_xyz=obj.xyz_target;
        end
        
        function loadCAD(obj)
            cad=load('./accOpt/Utils/CAD_model/LRMate200iD7L_20180211.mat');
            obj.parts(1).data.faces=cad.base{1}.f;
            obj.parts(1).data.vertices=cad.base{1}.v;
            obj.parts(1).id=1;
            obj.parts(1).color=cad.base{1}.color;
            for i=1:length(cad.link)
                obj.parts(i+1).data.faces=cad.link{i}.f;
                obj.parts(i+1).data.vertices=cad.link{i}.v;
                obj.parts(i+1).id=i+1;
                obj.parts(i+1).color=cad.link{i}.color;
            end
            obj.parts(end+1).data.faces=cad.payload{1}.f;
            obj.parts(end).data.vertices=cad.payload{1}.v;
            obj.parts(end).id=length(obj.parts)-1;
            obj.parts(end).color=cad.payload{1}.color;
        end
        
        function anim(obj,Vm,Am,dt)
            obj.dt=dt;
            obj.Vm=Vm;
            obj.Am=Am;
            obj.v_threshold=Am*dt*0.5;
            obj.enb=true;
            while(obj.sta)
                obj.stateUp;
            end
        end
        
        function stateUp(obj)
            [jntT]=kfwd_rob_LRMate_full(obj.jnt_pos,obj.si.DH,obj.rtb.tool,obj.T_B2W);
            pos_= jntT(1:3,4,end)+obj.delta;
            J = obj.jacobian();
            vel_=J*obj.jnt_vel;
            vel_=vel_(1:3);
            
            global x_rec xdot_rec q_rec qdot_rec
            x_rec=[x_rec,pos_(:)];
            xdot_rec=[xdot_rec,vel_(:)];
            q_rec=[q_rec,obj.jnt_pos(:)];
            qdot_rec=[qdot_rec,obj.jnt_vel(:)];
                
            [~,vel_]=OTG(obj.xyz_target,pos_,vel_,obj.Vm,obj.Am,obj.dt);
            if norm(vel_,inf) < obj.v_threshold
                vel_=vel_*0;
            end
            % follow xyz motion
            %obj.jnt_vel = pinv(J(1:3,:))*vel_(:);% min qdot move
            %obj.jnt_vel = obj.minWMov(J(1:3,:),vel_(:),diag(1./obj.r.G).^2);
            W=diag(1./(obj.r.pos_ub-obj.r.pos_lb).^2);
            obj.jnt_vel = obj.minWMov(J(1:3,:),vel_(:),W);
            
            % follow xyz motion, keep two angles
            %obj.jnt_vel = obj.minWMov(J(1:5,:),[vel_(:);0;0],diag(1./obj.r.G));
            
            obj.posUpdate();
            obj.draw();
            %pause(obj.dt*0.5);
        end
        
        function posUpdate(obj)
            obj.jnt_pos=obj.jnt_pos+...
                        obj.dt*obj.jnt_vel(:);
        end
        
        % move with joint weights, if redunent
        function qdot = minWMov(obj,J,xdot,W)
            s = pinv(J)*xdot;% basis joint vel
            N = null(J);% qdot=s+N*v
            if isempty(N)
                qdot = s;
            else
                v=-pinv(N.'*W*N)*N.'*W*s;
                qdot = s+N*v;
            end
        end
        
        function J = jacobian(obj)
            A=obj.si.DH(:,2);
            D=obj.si.DH(:,3);
            alpha=obj.si.DH(:,1);
            R_tool=obj.rtb.tool(1:3,1:3);
            T_tool=obj.rtb.tool(1:3,4);
            % Kinematics
            q=obj.jnt_pos;
            for i=1:6
                q(i)=q(i)+obj.si.DH(i,4);
            end
            % manully loop for each joint
            w=nan(3,6);
            r_0=nan(3,7);
            J=nan(6,6);% jacobian
            TCP_T=obj.T_B2W;
            for i=1:6
                w(:,i)=TCP_T(1:3,3);
                r_0(:,i)=TCP_T(1:3,4);
                TCP_T=TCP_T*...
                    [cos(q(i)) -sin(q(i))*cos(alpha(i))  sin(q(i))*sin(alpha(i)) A(i)*cos(q(i));...
                    sin(q(i))  cos(q(i))*cos(alpha(i)) -cos(q(i))*sin(alpha(i)) A(i)*sin(q(i));...
                    0            sin(alpha(i))                cos(alpha(i))            D(i);...
                    0                0                       0               1];
            end
            r_0(:,7)=TCP_T(1:3,4);
            TCP_T=TCP_T*... % up to now, forward kinematics, S0 to S6
                [R_tool,T_tool;0 0 0 1]; % S6 to tool
            TCP=TCP_T(1:3,4);%orig_abs;% TCP point in the world frame
            for i=1:6
                J(:,i)=[cross(r_0(:,i)-TCP,w(:,i));w(:,i)];
            end
        end
        
        function draw(obj)
            [jntT]=kfwd_rob_LRMate_full(obj.jnt_pos,obj.si.DH,obj.rtb.tool,obj.T_B2W);
            for j=1:length(obj.parts)
                set(obj.h_rob(j),'vertices',bsxfun(@plus,...
                    obj.parts(j).data.vertices*...
                    jntT(1:3,1:3,obj.parts(j).id).',...
                    jntT(1:3,4,obj.parts(j).id).'));
            end
            drawnow;
        end
        
        function InitPlot(obj)
            % axis limit
            xlimit=[0,1];
            ylimit=[0,1];
            zlimit=[0,1];
            hold on;
            viewang=[128,15];
            [jntT]=kfwd_rob_LRMate_full(obj.jnt_pos(:),...
                                     obj.si.DH,...
                                     obj.rtb.tool,obj.T_B2W);
            partsnum=numel(obj.parts);
            obj.h_rob=nan(partsnum,1);
            for j=1:partsnum
                obj.h_rob(j)=patch('faces',obj.parts(j).data.faces,...
                    'vertices',bsxfun(@plus,obj.parts(j).data.vertices*...
                    jntT(1:3,1:3,obj.parts(j).id).',...
                    jntT(1:3,4,obj.parts(j).id).'));
                set(obj.h_rob(j),'facecolor',obj.parts(j).color);
                set(obj.h_rob(j),'edgecolor','none');
                set(obj.h_rob(j),'Clipping','off');
            end
            [X,Y,Z]=sphere();
            obj.h_marker=mesh(obj.scale*X+jntT(1,4,end),...
                              obj.scale*Y+jntT(2,4,end),...
                              obj.scale*Z+jntT(3,4,end));
            set(obj.h_marker,'facecolor','c');
            set(obj.h_marker,'edgecolor','none');
            set(obj.h_marker,'Clipping','off');
            patch('faces',[1,2,3,4],'vertices',...
                [1,1,0;1,-1,0;-1,-1,0;-1,1,0]*1.5,...
                'facealpha',0.2);
            % set lighting, view, etc
            axis equal;
            view(viewang);
            xlim(xlimit);ylim(ylimit);zlim(zlimit);
            axis vis3d;
            box on;grid on;
            set(gca,'boxstyle','full');
            lighting flat;
            camlight right;
            set(gca,'ClippingStyle','rectangle');
            set(gca,'Clipping','off');
            axis off;
            set(gca,'Position',[-0 -0, 1 1]);
            % init 
            obj.moveit3();
            % modification
            obj.delta=0*obj.xyz_target;
            temp=[mean(mean(get(obj.h_marker,'XData')));...
                  mean(mean(get(obj.h_marker,'YData')));...
                  mean(mean(get(obj.h_marker,'ZData')))];
            dis=obj.xyz_target-temp;
            
            set(obj.h_marker,'XData',get(obj.h_marker,'XData') + dis(1));
            set(obj.h_marker,'YData',get(obj.h_marker,'YData') + dis(2));
            set(obj.h_marker,'ZData',get(obj.h_marker,'ZData') + dis(3));
            
            %posold=obj.xyz_target;
            %obj.xyz_target=[mean(mean(get(obj.h_marker,'XData')));...
            %                mean(mean(get(obj.h_marker,'YData')));...
            %                mean(mean(get(obj.h_marker,'ZData')))];
            %obj.delta=obj.xyz_target-posold;
            
            %set(gcf,'color','w');% change background color
            % plot bouncary limit
            
        end
        
        function moveit3(obj)
            % Unpack gui object
            gui = get(gcf,'UserData');
            % Make a fresh figure window
            set(obj.h_marker,'ButtonDownFcn',@obj.startmovit);
            set(gcf,'KeyPressFcn',@obj.keypre);
            set(gcf,'KeyReleaseFcn',@obj.keyrel);
            % Store gui object
            set(gcf,'UserData',gui);
        end
        
        function startmovit(obj,src,evnt)
            % Unpack gui object
            gui = get(gcf,'UserData');
            % Remove mouse pointer
            set(gcf,'PointerShapeCData',nan(16,16));
            set(gcf,'Pointer','custom');
            % Set callbacks
            gui.currenthandle = src;
            thisfig = gcbf();
            set(thisfig,'WindowButtonMotionFcn',@obj.movit);
            set(thisfig,'WindowButtonUpFcn',@obj.stopmovit);
            % Store starting point of the object
            gui.startpoint = get(gca,'CurrentPoint');
            set(gui.currenthandle,'UserData',{get(gui.currenthandle,'XData') get(gui.currenthandle,'YData') ...
                get(gui.currenthandle,'ZData')});
            % Store gui object
            set(gcf,'UserData',gui);
        end
        
        function keypre(obj,src,evnt)
            switch evnt.Key
                case 'k'
                    obj.sta=false;
                    obj.enb=false;
                case 'q'
                    obj.UDPenable();
            end
        end
        
        function sendUDP(obj)
            judp('send',25000,'192.168.1.100',typecast(...
                obj.xyz_target(:).'-obj.home_xyz(:).','int8'));
        end
        function UDPenable(obj)   
            obj.udp_enable = ~obj.udp_enable;
            disp(['udp status: ', num2str(obj.udp_enable)]);
        end
        
        function keyrel(obj,src,evnt)
            obj.sta=true;
        end
        
        function movit(obj,src,evnt)
            % Unpack gui object
            gui = get(gcf,'UserData');
            try
                if isequal(gui.startpoint,[])
                    return
                end
            catch
            end
            % Do "smart" positioning of the object, relative to starting point...
            pos = get(gca,'CurrentPoint')-gui.startpoint;
            XYData = get(gui.currenthandle,'UserData');
            offsetX=mean(mean(XYData{1}))-obj.delta(1);
            offsetY=mean(mean(XYData{2}))-obj.delta(2);
            offsetZ=mean(mean(XYData{3}))-obj.delta(3);
            % ball boundary
            posVec=[offsetX;offsetY;offsetZ]...
                   +[pos(1,1);pos(1,2);pos(1,3)]...
                   -obj.bound_c(:);% vector to center of boundary ball
            posR=norm(posVec);
            posU=posVec/norm(posVec);
            posR=obj.boundVal(posR,obj.bound_rin,obj.bound_rout);% limit target pos
            posVec=obj.bound_c(:)+posR*posU;
            % zlim
            if posVec(3)<obj.bound_c(3)+obj.z_bound(1)
                posVec(3)=obj.bound_c(3)+obj.z_bound(1);
            elseif posVec(3)>obj.bound_c(3)+obj.z_bound(2)
                posVec(3)=obj.bound_c(3)+obj.z_bound(2);
            end
            % ylim
            if posVec(2)<obj.bound_c(2)+obj.y_bound(1)
                posVec(2)=obj.bound_c(2)+obj.y_bound(1);
            elseif posVec(2)>obj.bound_c(2)+obj.y_bound(2)
                posVec(2)=obj.bound_c(2)+obj.y_bound(2);
            end
            % xlim
            if posVec(1)<obj.bound_c(1)+obj.x_bound(1)
                posVec(1)=obj.bound_c(1)+obj.x_bound(1);
            elseif posVec(1)>obj.bound_c(1)+obj.x_bound(2)
                posVec(1)=obj.bound_c(1)+obj.x_bound(2);
            end
            pos(1,1)=posVec(1)-offsetX;
            pos(1,2)=posVec(2)-offsetY;
            pos(1,3)=posVec(3)-offsetZ;
            % box boundary
            %pos(1,1)=obj.boundVal(pos(1,1),...
            %                      obj.Xboundary(1)-offsetX,...
            %                      obj.Xboundary(2)-offsetX);
            %pos(1,2)=obj.boundVal(pos(1,2),...
            %                      obj.Yboundary(1)-offsetY,...
            %                      obj.Yboundary(2)-offsetY);
            %pos(1,3)=obj.boundVal(pos(1,3),...
            %                      obj.Zboundary(1)-offsetZ,...
            %                      obj.Zboundary(2)-offsetZ);
            if ~obj.enb
            else
            set(gui.currenthandle,'XData',XYData{1} + pos(1,1));
            set(gui.currenthandle,'YData',XYData{2} + pos(1,2));
            set(gui.currenthandle,'ZData',XYData{3} + pos(1,3));
            % update target
            obj.xyz_target=[mean(mean(get(obj.h_marker,'XData')));...
                            mean(mean(get(obj.h_marker,'YData')));...
                            mean(mean(get(obj.h_marker,'ZData')))];
            if obj.udp_enable
                obj.sendUDP();
            end
            % draw
            drawnow;
            % update robot pos
            obj.stateUp();
            % Store gui object
            set(gcf,'UserData',gui);
            end
        end
        
        function stopmovit(obj,src,evnt)
            % Clean up the evidence ...
            thisfig = gcbf();
            gui = get(gcf,'UserData');
            set(gcf,'Pointer','arrow');
            set(thisfig,'WindowButtonUpFcn','');
            set(thisfig,'WindowButtonMotionFcn','');
            drawnow;
            set(gui.currenthandle,'UserData','');
            set(gcf,'UserData',[]);
        end
        
        function plotBound(obj)
            alpha=0.1;
            facecolor='none';
            edgecolor='k';
            n1=20;
            [t1,t2]=meshgrid(linspace(-pi,pi,n1),linspace(0,pi/2,n1));
            X=obj.bound_c(1)+obj.bound_rout*cos(t2).*cos(t1);
            Y=obj.bound_c(2)+obj.bound_rout*cos(t2).*sin(t1);
            Z=obj.bound_c(3)+obj.bound_rout*sin(t2);
            mesh(X,Y,Z,'facealpha',alpha,...
                'facecolor',facecolor,'edgecolor',edgecolor);
            n2=20;
            [t1,t2]=meshgrid(linspace(-pi,pi,n2),linspace(0,pi/2,n2));
            X=obj.bound_c(1)+obj.bound_rin*cos(t2).*cos(t1);
            Y=obj.bound_c(2)+obj.bound_rin*cos(t2).*sin(t1);
            Z=obj.bound_c(3)+obj.bound_rin*sin(t2);
            mesh(X,Y,Z,'facealpha',alpha,...
                'facecolor',facecolor,'edgecolor',edgecolor);
            n3=20;
            [r,t]=meshgrid(linspace(obj.bound_rin,obj.bound_rout,n3),...
                linspace(-pi,pi,n3));
            X=obj.bound_c(1)+r.*cos(t);
            Y=obj.bound_c(2)+r.*sin(t);
            Z=obj.bound_c(3)*ones(n3);
            mesh(X,Y,Z,'facealpha',alpha,...
                'facecolor',facecolor,'edgecolor',edgecolor);
            faces=[5 6 7 8;2 6 7 3;3 7 8 4;4 8 5 1;1 5 6 2;1 2 3 4];
            vertices=[obj.x_bound(2) obj.y_bound(1) obj.z_bound(2);...
                obj.x_bound(2) obj.y_bound(2) obj.z_bound(2);...
                obj.x_bound(1) obj.y_bound(2) obj.z_bound(2);...
                obj.x_bound(1) obj.y_bound(1) obj.z_bound(2);...
                obj.x_bound(2) obj.y_bound(1) obj.z_bound(1);...
                obj.x_bound(2) obj.y_bound(2) obj.z_bound(1);...
                obj.x_bound(1) obj.y_bound(2) obj.z_bound(1);...
                obj.x_bound(1) obj.y_bound(1) obj.z_bound(1)];
            for i=1:3
                vertices(:,i)=vertices(:,i)+obj.bound_c(i);
            end            
            patch('faces',faces,'vertices',vertices,...
                'facealpha',0.3,...
                'facecolor','none','edgecolor','k');
        end
        
    end
    
    methods(Static)
        
        function val = boundVal(val,lb,ub)
            if val>ub
                val=ub;
            elseif val<lb
                val=lb;
            end
        end
        
    end
end

