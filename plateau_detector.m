function [F,Fdist,points,fcns]=plateau_detector(t, X, varargin)
%PLATEAU_DETECTOR A plateau detector for oscillating timeseries. Uses Shmitt
% trigger concept - upward transition threshold >= downward transition
% threshold. The timeseries is divided into periods measured from upward
% transition to upward transition.
%
% If results are plotted, use left/right arrows to switch between columns
% of X.
%
% Required Inputs:
%  t - vector of timepoints
%  X - vector matching t, or matrix with columns matching t
%
% Optional input:
%  f - fraction of amplitude for threshold:
%     scalar {default, f=0.5}, upward=downward threshold,
%     vector, f=[upFraction,downFraction] with upFraction>=downFraction.
%
% optional name-value pair arguments:
%  'Interpolate', {true}/false
%  'ThresholdPercentiles', default=[0,100] (threshold is fraction of this range of amplitudes)
%  'Plot',true/{false} (use to control plotting with output arguments)
%
% Outputs:
%  features - struct containg feature distributions for each column of X
%  points - struct containing special points, for plotting


%TODO: no inputs => output list of feature names, fcn handles (compute_features+plot), etc?

%TODO: optionally set period marker to other than "up" points? eg. periodMarker='min';

%TODO: make sure initialization of upstate toggle is robust

%TODO: add slope thresholds as option (good for voltage?)
%TODO: noise-robust method for numerical slope

%TODO: option to find maxslope only between up/down. Otherwise anywhere in period. only weird because slope not zero at
%up transition; min-min period avoids this.

%TODO: scalar feature summary output

%TODO: per-period feature idea
%       - phase of special point within current period (eg. phase of max V1 relative to period detected in V2)


if isvector(X)
    X=X(:);
end

[thrPtiles,fracUp,fracDown,fracUpDx,fracDownDx,minAmp,threshMethod,doInterp,doPlot,figID,dokeypress]=parseArgs(t, X, varargin{:});

nX=size(X,2); %number of traces

DX=slopeY(t,X);

globalXmin=prctile(X,thrPtiles(1),1);
globalXmax=prctile(X,thrPtiles(2),1);
globalXamp=globalXmax-globalXmin;

globalDXmin=prctile(DX,thrPtiles(1),1);
globalDXmax=prctile(DX,thrPtiles(2),1);
% globalDXamp=globalDXmax-globalDXmin;

%independent thresholds for each trace
switch lower(threshMethod)
    case 'fracamp'
        %fraction of global amplitude method
        thrUp=globalXmin + fracUp*globalXamp;
        if fracDown>0
            thrDown=globalXmin + fracDown*globalXamp;
        else
            thrDown=thrUp;
        end
        
    case 'rawvals'
        %specific raw values as thresholds
        thrUp=fracUp;
        thrDown=fracDown;
        
    case 'meanfrac'
        %mean +/- fraction of amp
        meanX=mean(X,1);
        thrUp=meanX + fracUp*globalXamp;
        thrDown=meanX - fracDown*globalXamp;
        
    case 'medfrac'
        %median +/- fraction of amp
        medX=median(X,1);
        thrUp=medX + fracUp*globalXamp;
        thrDown=medX - fracDown*globalXamp;
        
    case 'medraw'
        %median +/- raw values
        medX=median(X,1);
        thrUp=medX + fracUp;
        thrDown=medX - fracDown;
end
thrUpDx=fracUpDx*globalDXmax;
if fracDownDx>0
    thrDownDx=fracDownDx*globalDXmin;
else
    thrDownDx=globalDXmin;
end

pt=struct('ix',[],'t',[],'x',[],'dx',[]);

points=repmat(struct('period',pt,'up',pt,'down',pt,...
    'max',pt,'min',pt,'dxmax',pt,'dxmin',pt),1,nX);

sufficientGlobalAmp=globalXamp>minAmp;

% isUp=ones(1,nX); %forces first point to be downward transition
% isUp=zeros(1,nX); %forces first point to be upward transition
isUp=X(1,:)>thrUp;
% isUp=X(1,:)>=thrUp | (X(1,:)<=thrUp&X(1,:)>=thrDown&DX(1,:)<0);

for i=2:length(t)
    
    
    %check for downward-transition and upward-transition
    %     isDownTransition=isUp & X(i,:)<=thrDown;
    %     isUpTransition=~isUp & X(i,:)>=thrUp;
    isDownTransition=isUp & X(i,:)<=thrDown & DX(i,:)>=thrDownDx;
    isUpTransition=~isUp & X(i,:)>=thrUp & DX(i,:)>=thrUpDx;
    
    isDownTransition(~sufficientGlobalAmp)=0;
    isUpTransition(~sufficientGlobalAmp)=0;
    
    %if started up, then went below thrUp but not below thrDown, set isUp=0
    firstDown=isDownTransition & arrayfun(@(x)isempty(x.up.ix),points);
    isUp(firstDown)=0;
    
    if any(isDownTransition)
        idx=find(isDownTransition);
        for j=idx
            %index
            points(j).down.ix(end+1)=i;
            
            if doInterp
                %interpolate t
                thisX=thrDown(j);
                thisDX=thrDownDx(j);
                thisT=interp1(X(i-1:i,j),t(i-1:i),thisX);
                if isnan(thisT)
                    thisT=interp1(DX(i-1:i,j),t(i-1:i),thisDX);
                end
                %                     thisT=(thrDown(j)-X(i-1,j))*(t(i)-t(i-1))/(X(i,j)-X(i-1,j))+t(i-1);
            else
                thisT=t(i);
                thisX=X(i,j);
            end
            
            points(j).down.t(end+1)=thisT;
            points(j).down.x(end+1)=thisX;
        end
        isUp(isDownTransition)=0;
    end
    
    if any(isUpTransition)
        idx=find(isUpTransition);
        for j=idx
            points(j).up.ix(end+1)=i;
            
            if doInterp
                %interpolate t
                thisX=thrUp(j);
                thisDX=thrUpDx(j);
                thisT=interp1(X(i-1:i,j),t(i-1:i),thisX);
                if isnan(thisT)
                    thisT=interp1(DX(i-1:i,j),t(i-1:i),thisDX);
                end
                %                     thisT=(thrUp(j)-X(i-1,j))*(t(i)-t(i-1))/(X(i,j)-X(i-1,j))+t(i-1);
            else
                thisT=t(i);
                thisX=X(i,j);
            end
            
            points(j).up.t(end+1)=thisT;
            points(j).up.x(end+1)=thisX;
        end
        isUp(isUpTransition)=1;
    end
    
end

%clean up period maker points
for i=1:nX
    points(i).period.ix=points(i).up.ix;
    points(i).period.t=points(i).up.t;
    points(i).period.x=points(i).up.x;
    
    if ~isempty(points(i).down.t) && ~isempty(points(i).up.t)
        % trim away down transitions outside first/last up transitions
        if points(i).down.t(1)<points(i).up.t(1)
            points(i).down.ix=points(i).down.ix(2:end);
            points(i).down.t=points(i).down.t(2:end);
            points(i).down.x=points(i).down.x(2:end);
        end
    end
    if ~isempty(points(i).down.t) && ~isempty(points(i).up.t)
        if points(i).down.t(end)>points(i).up.t(end)
            points(i).down.ix=points(i).down.ix(1:end-1);
            points(i).down.t=points(i).down.t(1:end-1);
            points(i).down.x=points(i).down.x(1:end-1);
        end
    end
    
    
end


%compute the extra points and features - per-period distributions
[F,Fdist,points]=compute_features(t,X,points,DX);



fcns.compute_features=@compute_features;
fcns.plot_data=@plot_data; %simple plot fcn for one trace of interest
fcns.plot_interactive=@plot_interactive; %adds keypressfcn to switch traces

%plot to show performance
if nargout==0||doPlot==1
    if dokeypress
        plot_interactive(figID,t,X,points)
    else
        if isempty(figID)
            gcf; %new figure if none available, otherwise use current fig
        else
            figure(figID); %bring figID into focus
        end
        plot_data(t,X,points,1)
    end
end

end

%up/down times are retained, all other points/features computed from those
function [F,Fdist,points]=compute_features(t,X,points,DX)

if ~exist('DX','var')
    DX=slopeY(t,X);
end

nX=size(X,2);

Fdist=repmat( struct('period',0,'apd',0,'spd',0,'pf',0,'f1',0,'f2',0,'f3',0,'f4',0,'amp',0,...
    'baseline',0,'peaks',0,'maxslope',0,'minslope',0) ,1,nX); %'range',[],'thrUp',[],'thrDown',[],'pthresh',[]

for i=1:nX
    
    %find the x value at up and down marker times (up is period marker)
    %     points(i).up.x=interp1(t,X(:,i),points(i).up.t);
    points(i).period.x=points(i).up.x;
    %     points(i).down.x=interp1(t,X(:,i),points(i).down.t);
    
    nPer=length(points(i).up.t)-1;
    
    if nPer>=1
        
        for j=1:nPer
            ix=find(t>=points(i).period.t(j) & t<=points(i).period.t(j+1));
            tt=t(ix);
            xx=X(ix,i);
            dx=DX(ix,i);
            
            [xmax,imax]=max(xx);
            points(i).max.ix(j)=ix(1)+imax-1;
            points(i).max.t(j)=tt(imax);
            points(i).max.x(j)=xmax;
            
            [xmin,imin]=min(xx);
            points(i).min.ix(j)=ix(1)+imin-1;
            points(i).min.t(j)=tt(imin);
            points(i).min.x(j)=xmin;
            
            [dxmin,idxmin]=min(dx);
            points(i).dxmin.ix(j)=ix(1)+idxmin-1;
            points(i).dxmin.t(j)=tt(idxmin);
            points(i).dxmin.x(j)=xx(idxmin);
            points(i).dxmin.dx(j)=dxmin;
            
            %compute maxslope only in the active phase? Not if finding max/min in a different trace than periods
            ix=find(t>=points(i).up.t(j) & t<=points(i).down.t(j));
            xx=X(ix,i);
            dx=DX(ix,i);
            [dxmax,idxmax]=max(dx);
            points(i).dxmax.ix(j)=ix(1)+idxmax-1;
            points(i).dxmax.t(j)=tt(idxmax);
            points(i).dxmax.x(j)=xx(idxmax);
            points(i).dxmax.dx(j)=dxmax;
        end
        
        %these are scalar features:
        %     Fdist(i).range=globalXamp(i);
        %     Fdist(i).thrUp=thrUp(i);
        %     Fdist(i).thrDown=thrDown(i);
        %     Fdist(i).pthresh=(thrUp(i)+thrDown(i))/2;
        
        Fdist(i).period=diff(points(i).up.t);
        Fdist(i).apd=points(i).down.t-points(i).up.t(1:end-1); %active phase duration: we stripped any extra down transitions
        Fdist(i).spd=Fdist(i).period-Fdist(i).apd;
        Fdist(i).pf=Fdist(i).apd./Fdist(i).period;
        Fdist(i).f1=(points(i).max.t-points(i).up.t(1:end-1))./Fdist(i).period; %up2max time
        Fdist(i).f2=(points(i).down.t-points(i).max.t)./Fdist(i).period; %max2down time
        Fdist(i).f3=(points(i).min.t-points(i).down.t)./Fdist(i).period; %down2min time
        Fdist(i).f4=(points(i).up.t(2:end)-points(i).min.t)./Fdist(i).period; %min2up time
        Fdist(i).baseline=points(i).min.x;
        Fdist(i).peaks=points(i).max.x;
        Fdist(i).amp=points(i).max.x-points(i).min.x;
        Fdist(i).maxslope=points(i).dxmax.dx;
        Fdist(i).minslope=points(i).dxmin.dx;
        
    else
        %had less than two up transitions: can't compute features.
        if ~isempty(points(i).min.x)
            Fdist(i).baseline=points(i).min.x;
        end
        if ~isempty(points(i).max.x)
            Fdist(i).peaks=points(i).max.x;
        end
        if ~isempty(points(i).max.x) && ~isempty(points(i).min.x)
            Fdist(i).amp=points(i).max.x-points(i).min.x;
        end
        if ~isempty(points(i).dxmax.dx)
            Fdist(i).maxslope=points(i).dxmax.dx;
        end
        if ~isempty(points(i).dxmin.dx)
            Fdist(i).minslope=points(i).dxmin.dx;
        end
    end
end

%compute summary statistics of feature distributions
fnames=fieldnames(Fdist);

for i=1:nX
    for j=1:length(fnames)
        F(i).([fnames{j},'Mean'])=mean(Fdist(i).(fnames{j}));
        %         F(i).([fnames{j},'Stdev'])=std(Fdist(i).(fnames{j}));
    end
end

end

function plot_data(t,X,points,tix)

plot(t,X(:,tix),'k-')
hold on

plot(points(tix).up.t,points(tix).up.x,'bs')
plot(points(tix).down.t,points(tix).down.x,'bo')
plot(points(tix).max.t,points(tix).max.x,'rv')
plot(points(tix).min.t,points(tix).min.x,'r^')
plot(points(tix).dxmax.t,points(tix).dxmax.x,'g>')
plot(points(tix).dxmin.t,points(tix).dxmin.x,'g<')
% plot(xlim(),thrUp(tix)*[1,1],'r--')
% if thrUp~=thrDown
%     plot(xlim(),thrDown(tix)*[1,1],'b--')
% end
plot(points(tix).period.t,points(tix).period.x,'bs')
hold off

xlabel('t')
ylabel('x')
axis tight
YLIM=ylim();
ylim([YLIM(1)-0.05*abs(YLIM(1)),YLIM(2)+0.05*abs(YLIM(2))]);
end

function plot_interactive(figID,t,X,points)
tix=1;
if isempty(figID)
    figID=gcf; %new figure if none available, otherwise current fig
else
    figID=figure(figID);
end
figID.KeyPressFcn=@keypressFcn;
figID.UserData.tix=tix; %store tix (trace to plot) in userdata
figID.UserData.t=t;
figID.UserData.X=X;
figID.UserData.points=points;

plot_data(t,X,points,tix)
end

function keypressFcn(src,event)
tix=src.UserData.tix;
t=src.UserData.t;
X=src.UserData.X; nX=size(X,2);
points=src.UserData.points;
switch(event.Key)
    case {'leftarrow'}
        if tix>1
            tix=tix-1;
        end
    case {'rightarrow'}
        if tix<nX
            tix=tix+1;
        end
end
plot_data(t,X,points,tix)
src.UserData.tix=tix;
end

function [thrPtiles,fracUp,fracDown,fracUpDx,fracDownDx,minAmp,threshMethod,doInterp,doPlot,figID,dokeypress]=parseArgs(t, X, varargin)

%default parameters
defaultF=[0.5,0.4]; %near halfmax
defaultFslope=[0,0];
defaultthrPtiles=[0,100];
defaultminAmp=0;
defaultThresholdMethod='fracamp';
doInterp=0;
doPlot=false;
doKeypress=true;
figID=[];

p=inputParser;
validX=@(x) isreal(x) && size(x,1)==length(t); %traces are columns of X
validSwitch=@(x) isscalar(x) && (isnumeric(x)||islogical(x));
addRequired(p,'t',@(x) isreal(x));
addRequired(p,'X',validX);
addOptional(p,'f',defaultF);  %somehow adding the validation function here messes things up
addParameter(p,'ThresholdPercentiles',defaultthrPtiles,@(x) isreal(x) && numel(x)==2);
addParameter(p,'MinimumAmplitude',defaultminAmp,@isreal);
addParameter(p,'ThresholdMethod',defaultThresholdMethod);
addParameter(p,'Interpolate',doInterp,validSwitch);
addParameter(p,'Plot',doPlot,validSwitch);
addParameter(p,'Keypress',doKeypress,validSwitch);
addParameter(p,'FigureID',figID);

parse(p,t,X,varargin{:});

f=p.Results.f;
thrPtiles=p.Results.ThresholdPercentiles;
minAmp=p.Results.MinimumAmplitude;
threshMethod=p.Results.ThresholdMethod;
doInterp=p.Results.Interpolate;
doPlot=p.Results.Plot;
figID=p.Results.FigureID;
dokeypress=p.Results.Keypress;

fracUpDx=defaultFslope(1);
fracDownDx=defaultFslope(2);
if isscalar(f)
    fracUp=f;
    fracDown=f;
elseif length(f)==2
    fracUp=f(1);
    fracDown=f(2);
elseif length(f)==4
    fracUp=f(1);
    fracDown=f(2);
    fracUpDx=f(3);
    fracDownDx=f(4);
else
    error('Invalid value for fractional threshold')
end

%     %snap to fracUp?
%     if fracDown>fracUp
%         fracDown=fracUp;
%     end

end