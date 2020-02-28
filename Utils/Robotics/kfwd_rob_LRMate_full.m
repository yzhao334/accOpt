%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Analytical forward kinematics of LRMate robot     %
%                   Yu Zhao, 2018/02/11                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Take one group of theta one time

% Inputs:   theta       joint position in rad unit
%           DH          DH parameter of robot, [m], [rad]
%           tool        tool frame relative to link 6 frame
%           base        base frame transform
% Outputs:  jntT         transform of each joint
function [ jntT ] = kfwd_rob_LRMate_full( theta, DH, tool, base )
%% DH parameters
alpha=DH(:,1);
A=DH(:,2);
D=DH(:,3);
offset=DH(:,4);
%%
jntT = zeros(4,4,8);
jntT(:,:,1)=base;
q = theta(:) + offset;
TCP_T=jntT(:,:,1);
for i=1:6
    TCP_T=TCP_T*...
        [cos(q(i)) -sin(q(i))*cos(alpha(i))  sin(q(i))*sin(alpha(i)) A(i)*cos(q(i));...
         sin(q(i))  cos(q(i))*cos(alpha(i)) -cos(q(i))*sin(alpha(i)) A(i)*sin(q(i));...
          0            sin(alpha(i))                cos(alpha(i))            D(i);...
          0                0                       0               1];
    jntT(:,:,i+1) = TCP_T;
end
TCP_T=TCP_T*... % up to now, forward kinematics, S0 to S6
    double(tool); % S6 to tool
jntT(:,:,end) = TCP_T;

end

