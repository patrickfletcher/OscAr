function DY = slopeY(t,Y,boundaryMethod)
% central difference approximation for the slope of Y that handles unequal
% spacing.

%TODO: 
% - add options for forward/backward differences, different orders, low-noise versions
% - add options for different boundary strategies


% centDiff=@(xp1,xm1) (xp1-xm1) / (hp1+h);

if ~exist('boundaryMethod','var')
    boundaryMethod='order1';
end

DY=zeros(size(Y));

%compute slope only at inner points
% for i=2:length(Y)-1
%     h=t(i)-t(i-1);
%     hp1=t(i+1)-t(i);
%     DY(i,:)=(Y(i+1,:)-Y(i-1,:))/(hp1+h);
% end

%vectorized:
h=t(2:length(Y)-1)-t(1:length(Y)-2);
hp1=t(3:length(Y))-t(2:length(Y)-1);
DY(2:length(Y)-1,:)=(Y(3:length(Y),:)-Y(1:length(Y)-2,:))./(hp1+h);

%use difference that works at ends
switch boundaryMethod
    case 'order1'
    DY(1,:)=(Y(2,:)-Y(1,:))/(t(2)-t(1));
    DY(end,:)=(Y(end,:)-Y(end-1,:))/(t(end)-t(end-1));

    case 'clamp'
    %clamp the endpoints to the second values
    DY(1,:)=DY(2,:); 
    DY(end,:)=DY(end-1,:);
end

end