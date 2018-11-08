function frameplot( R, P, length )
% plot frame
R=R*length;
hold on;
% h=nan(3,1);
% xaxis
quiver3(P(1),P(2),P(3),R(1,1),R(2,1),R(3,1),'r');
% yaxis
quiver3(P(1),P(2),P(3),R(1,2),R(2,2),R(3,2),'g');
% zaxis 
quiver3(P(1),P(2),P(3),R(1,3),R(2,3),R(3,3),'b');


end

