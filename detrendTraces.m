function [Xdetrend,XTrend]=detrendTraces(t,X,method,methodparam,doPlot)
% DETRENDTRACES detrend the columns of X using various methods

% TODO: finish error checking, documentation

% TODO: no inputs - return cell array of possible methods with their possible params
%  {{method},{methodpar}}

nX=size(X,2);

dt=mode(diff(t));
fs=1/dt; 

if ~exist('doPlot','var')
    doPlot=0;
end

steepness=0.5;
stopAtten=90;

XTrend=zeros(size(X));
switch lower(method)
    case {'none'}
%         XTrend=zeros

    case {'linear'}
        
        for i=1:nX
            XTrend(:,i)=polyval(polyfit(t,X(:,i),1),t);
        end
        
    case {'robustlin'}
        
        for i=1:nX
            brob = robustfit(t,X(:,i));
            XTrend(:,i)=brob(1)+brob(2)*t;
        end

    case {'poly','polynomial'}
        
        %methodparam=degree of polynomial
        if ~exist('methodparam','var')||isempty(methodparam)
            degree=1;
        else
            degree=methodparam;
        end
        
        for i=1:nX
            XTrend(:,i)=polyval(polyfit(t,X(:,i),degree),t);
        end

    case {'ptile'}
        %use a moving window prctile
        %methodparam=window width
        if ~exist('methodparam','var')||isempty(methodparam)
            error('ptile trendline requires a window duration (in time units) and a prctile');
        else
            wwidth=methodparam{1};
            ptile=methodparam{2};
        end
        if wwidth>(t(end)-t(1))
            wwidth=(t(end)-t(1));
            warning(['windowsize shrunk to time interval, ', num2str(wwidth) ])
        end
        wsz=round(wwidth/dt);
        wsz=max(wsz,1);
        wsz2=ceil(wsz/2);
        for i=wsz2:length(t)-wsz2
            ix=i-wsz2+1:i+wsz2;
            XTrend(i,:)=prctile(X(ix,:),ptile,1);
        end
        XTrend(1:wsz2-1,:)=repmat(XTrend(wsz2,:),wsz2-1,1);
        XTrend(end-wsz2+1:end,:)=repmat(XTrend(end-wsz2,:),wsz2,1);
        
    case {'lowpass'}
        
        %methodparam=lowpass cutoff frequency
        if ~exist('methodparam','var')||isempty(methodparam)
            error('lowpass trendline requires a cutoff frequency');
        else
            fpass=methodparam;
        end
        
%         XTrend=lowpass(X,fpass,fs,'ImpulseResponse','iir');
        XTrend=lowpass(X,fpass,fs,'ImpulseResponse','iir','Steepness',steepness,'StopbandAttenuation',stopAtten);
    
    case {'highpass'}
        
        %methodparam=highpass cutoff frequency
        if ~exist('methodparam','var')||isempty(methodparam)
            error('highpass trend removal requires a cutoff frequency');
        else
            fpass=methodparam;
        end
        
        Xdetrend=highpass(X,fpass,fs,'ImpulseResponse','iir','Steepness',steepness,'StopbandAttenuation',stopAtten);
        XTrend=X-Xdetrend;

    case {'movmean','movmedian','gaussian','lowess','loess','rlowess','rloess','sgolay'}
        
        %methodparam=window width
        if ~exist('methodparam','var')||isempty(methodparam)
            error('smoothdata trendline requires a window duration (in time units)');
        else
            wwidth=methodparam;
        end
        
        if wwidth>(t(end)-t(1))
            wwidth=(t(end)-t(1));
            warning(['windowsize shrunk to time interval, ', num2str(wwidth) ])
        end
        
        wsz=round(wwidth/dt);
        wsz=max(wsz,1);
        XTrend=smoothdata(X,method,wsz);
end

Xdetrend=X-XTrend;



%plot to show result
if nargout==0 || doPlot==1
    
nX=size(X,2);
tix=1;
figure('Name','Detrend Traces','KeyPressFcn',@keypressFcn);
plotData()
    
end


%nested functions can see variables in caller's scope
    function plotData()
        subplot(2,1,1)
        plot(t,X(:,tix),t,XTrend(:,tix))
        grid on
        xlabel('Time')
%         ylabel('raw')
        axis tight
        
        subplot(2,1,2)
        plot(t,Xdetrend(:,tix))
        grid on
        xlabel('Time')
        ylabel('detrend')
        axis tight
        
    end

    function keypressFcn(src,event)
        switch(event.Key)
            case {'leftarrow'}
                if tix>1
                    tix=tix-1;
                    plotData()
                end
            case {'rightarrow'}
                if tix<nX
                    tix=tix+1;
                    plotData()
                end
        end
        
    end

end