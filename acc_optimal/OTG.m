function [ pos_, vel_, acc_ ] = OTG( pd, pos, vel, Vm, Am, dt )
%Online trajectory generation with bounded vel and acc bounds
% @Yu Zhao
%   pd:     target position, default reach with 0 vel
%   pos:    current position
%   vel:    current velocity
%   Vm:     maximum velocity
%   Am:     maximum acceleration
%   dt:     sample time
dim = length(pd);
negDelta=pd(:)-pos(:);
ret = directionGeneral( negDelta );
% isotropy goalchanged axis limit
if length(Vm)==length(ret(:,1))
    VmI=min(abs(Vm(:))./abs(ret(:,1)));
elseif length(Vm)==1
    VmI=min(abs(Vm*ones(size(ret(:,1))))./abs(ret(:,1)));
end
if length(Am)==length(ret(:,1))
    AmI=min(abs(Am(:))./abs(ret(:,1)));
elseif length(Am)==1
    AmI=min(abs(Am*ones(size(ret(:,1))))./abs(ret(:,1)));
end
% otg code
delta_temp=pos(:)-pd(:);
tempv=zeros(dim,1);
tempa=zeros(dim,1);
for jj=1:dim
    x=[dot(delta_temp,ret(:,jj)),...
       dot(vel(:),ret(:,jj))];
    a=acc_control_optimal(x,VmI,...
                          AmI,dt);
    x_next=acc_dyn(x,a,dt);
    tempv(jj)=x_next(2);
    tempa(jj)=a;
end
vel_=ret*tempv;
acc_=ret*tempa;
pos_=pos(:)+vel(:)*dt+acc_(:)*dt^2/2;

function  [ ret ] = directionGeneral( Delta)
%get n orthonormal vectors from given displacement Delta
ret = RoughDirectionSet( Delta(:) );
n=size(ret,1);
ret(:,1)=ret(:,1)/norm(ret(:,1));
for ii=2:n
    vk=ret(:,ii);
    for jj=1:ii-1
        vk=vk-dot(ret(:,jj),vk)/norm(ret(:,jj))^2*ret(:,jj);
    end
    ret(:,ii)=vk/norm(vk);
end


function [ directionSet ] = RoughDirectionSet( Delta )
%from given delta generate rough directionset
%   directionSet: linear independent directions, start from resized Delta
v0 = Delta(:);
[~, ind] = max(abs(v0));
directionSet=zeros(length(v0));
directionSet(:,1)=v0;
% modify
if norm(v0,inf)==0
    directionSet(ind,1)=1;
end
for ii=ind+1:length(v0)
    directionSet(ii,ii+1-ind)=1;
end
m=length(v0)+1-ind;
for ii=1:ind-1
    directionSet(ii,ii+m)=1;
end

function [ a ] = acc_control_optimal( x, Vm, U, T)
yn=x(1);
ynd=x(2);
zn=1/(T*U)*(yn/T+ynd/2);
znd=ynd/(T*U);
m=floor((1+sqrt(1+8*abs(zn)))/2);
sigman=znd+zn/m+(m-1)/2*sign(zn);
a=-U*sat_fcn(sigman)*...
    (1+sign(...
    ynd*sign(sigman)+Vm-T*U...
    ))/2;

function [ y ]=sat_fcn(x)
y=min(max(x,-1),1);

function [ x_next ] = acc_dyn( x, a, dt )
ek=x(1);
vk=x(2);
ek_next=ek+vk*dt+1/2*a*dt^2;
vk_next=vk+a*dt;
x_next=[ek_next;vk_next];

