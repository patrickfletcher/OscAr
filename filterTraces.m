function [XFilt,tfilt]=filterTraces(t,X,method,methodparam,doPlot,opts)

arguments
    t
    X
    method
    methodparam
    doPlot=0
    opts.steepness=0.85;
    opts.stopAtten=60;
    opts.doPlot = 0
end

% TODO: no inputs - return cell array of possible methods with their possible params
%  {{method},{methodpar}}

nX=size(X,2);

dt=mode(diff(t));
fs=1/dt; 

steepness=opts.steepness;
stopAtten=opts.stopAtten;

tfilt=t;
switch lower(method)
    case {'none'}
        XFilt=X;
        

    case {'movmean','movmedian','gaussian','lowess','loess','rlowess','rloess','sgolay'}
        
        %methodparam=window width
        if ~exist('methodparam','var')
            error('smoothdata requires a window duration (in time units)');
        else
            wwidth=methodparam;
        end
        
        wsz=round(wwidth/dt);
        wsz=max(wsz,1);
        XFilt=smoothdata(X,method,wsz);

    case {'lowpass'}
        
        %methodparam=lowpass cutoff frequency
        if ~exist('methodparam','var')
            error('lowpass requires a cutoff frequency');
        else
            fpass=methodparam;
        end
        
%         XFilt=lowpass(X,fpass,fs);
%         XFilt=lowpass(X,fpass,fs,'ImpulseResponse','iir');
        XFilt=lowpass(X,fpass,fs,'ImpulseResponse','auto','Steepness',steepness,'StopbandAttenuation',stopAtten);
       
    case {'highpass'}
        
        %methodparam=lowpass cutoff frequency
        if ~exist('methodparam','var')
            error('lowpass requires a cutoff frequency');
        else
            fpass=methodparam;
        end
        
%         XFilt=highpass(X,fpass,fs);
%         XFilt=highpass(X,fpass,fs,'ImpulseResponse','iir');
        XFilt=highpass(X,fpass,fs,'ImpulseResponse','auto','Steepness',steepness,'StopbandAttenuation',stopAtten);
    
    case {'bandpass'}
        %TODO: edge effect corrections
        
        %methodparam=lowpass cutoff frequency
        if ~exist('methodparam','var')
            error('bandpass requires a low and high frequency cutoff');
        else
            fpass=methodparam;
        end
        
        XFilt=bandpass(X,fpass,fs);
        % XFilt=bandpass(X,fpass,fs,'ImpulseResponse','auto');
        % XFilt=bandpass(X,fpass,fs,'ImpulseResponse','auto','Steepness',steepness,'StopbandAttenuation',stopAtten);
        
    otherwise
        error(['unknown method: ' method]);
end


%plot to show result
if nargout==0 || doPlot==1
    
tix=1;
figure('Name','Filter Traces','KeyPressFcn',@keypressFcn);
plotData()
    
end


%nested functions can see variables in caller's scope
    function plotData()
        
        plot(t,X(:,tix))
        hold on
        plot(tfilt,XFilt(:,tix),'linewidth',1.5)
        hold off
        grid on
        xlabel('Time')
        ylabel('filtered')
        axis tight
        
    end

    function keypressFcn(~,event)
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