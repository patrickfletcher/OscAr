function XNorm=normalizeTraces(t, X, method, methodparam, doPlot, opts)
% NORMALIZETRACES centers and or scales X values by a variety of methods, columnwise

% TODO: no inputs - return cell array of possible methods with their possible params
%  {{method},{methodpar}}

% TODO: update plotting method to not rely on nested functions and be more efficient

arguments
    t
    X
    method
    methodparam
    doPlot=0
    opts.doPlot = 0
end

nX=size(X,2);

switch method
    case {'none'}
        XNorm=X;
        
    case {'zscore'}
        XNorm=(X-mean(X,1))./std(X,[],1);
        
    case {'unit'}
        XNorm=(X-min(X,[],1))./(max(X,[],1)-min(X,[],1));
        if exist('methodparam','var')&&length(methodparam)==2
            intrvl=methodparam;
            XNorm=intrvl(1)+diff(intrvl)*XNorm;
        end
        
    case {'ptile'} %methodparam=[lo,hi] ptile
        p=prctile(X,methodparam,1);
        XNorm=(X-p(1,:))./abs(p(2,:)-p(1,:));
        % XNorm(XNorm<0)=0;
        % XNorm(XNorm>1)=1;
        
    case {'devmean'}
        XNorm=(X-mean(X,1))./mean(X,1);
        
    case {'devmean2'}
        XNorm=(X-mean(X,1))./mean(X,1).^2;
        
    case {'devmeanpow'}
        XNorm=(X-mean(X,1))./mean(X,1).^methodparam;
        
    case {'devmedian'}
        XNorm=(X-median(X,1))./median(X,1);
        
    case {'devmedianpow'}
        XNorm=(X-median(X,1))./median(X,1).^methodparam;
        
    case {'devtrend','trend'}
        
        %no errorchecking for now..
        trendmethod=methodparam{1};
        trendpar=methodparam{2};
        
        [~,XTrend]=detrendTraces(t,X,trendmethod,trendpar);
        XNorm=(X-XTrend)./XTrend;
        
    case {'devtrend2','trend2'}
        
        %no errorchecking for now..
        trendmethod=methodparam{1};
        trendpar=methodparam{2};
        
        [~,XTrend]=detrendTraces(t,X,trendmethod,trendpar);
        XNorm=(X-XTrend)./XTrend.^2;
        
    case {'center'}
        switch methodparam
            case {'mean'}
                XNorm=X-mean(X,1);
                
            case {'median'}
                XNorm=X-median(X,1);
                
            case {'max'}
                XNorm=X-max(X,[],1);
                
            case {'min'}
                XNorm=X-min(X,[],1);
                
            otherwise
                if isnumeric(methodparam) &&isscalar(methodparam)
                    k=methodparam;
                else
                    k=mean(X,1); %empty param defaults to mean
                end
                XNorm=X-k;
        end
        
    case {'scale'}
        switch methodparam
            case {'mean'}
                XNorm=X./mean(X,1);
                
            case {'stdev'}
                XNorm=X./std(X,[],1);
                
            case {'median'}
                XNorm=X./median(X,1);
                
            case {'max'}
                XNorm=X./max(X,[],1);
                
            case {'min'}
                XNorm=X./min(X,[],1);
                
            case {'range'}
                XNorm=X./(max(X,[],1)-min(X,[],1));
                
            case {'iqr'}
                XNorm=X./(prctile(X,25,1)-prctile(X,75,1));
                
            otherwise
                if isnumeric(methodparam) &&isscalar(methodparam)
                    k=methodparam;
                else
                    k=0;
                end
                XNorm=X/k;
        end
        
    case {'wptile','windowptile'}
        %use a moving window prctile
        %methodparam=window width
        if ~exist('methodparam','var')||isempty(methodparam)
            error('smoothdata trendline requires a window duration (in time units)');
        else
            wwidth=methodparam{1};
            ptile=methodparam{2}; %[lo,hi]
        end
        dt=mode(diff(t));
        wsz=round(wwidth/dt);
        wsz=max(wsz,1);
        wsz2=ceil(wsz/2);
        
        Xlo=zeros(size(X));
        Xhi=zeros(size(X));
        for i=wsz2:length(t)-wsz2
            ix=i-wsz2+1:i+wsz2;
            Xlo(i,:)=prctile(X(ix,:),ptile(1),1);
            Xhi(i,:)=prctile(X(ix,:),ptile(2),1);
        end
        Xlo(1:wsz2-1,:)=repmat(Xlo(wsz2,:),wsz2-1,1);
        Xlo(end-wsz2+1:end,:)=repmat(Xlo(end-wsz2,:),wsz2,1);
        Xhi(1:wsz2-1,:)=repmat(Xhi(wsz2,:),wsz2-1,1);
        Xhi(end-wsz2+1:end,:)=repmat(Xhi(end-wsz2,:),wsz2,1);
        
%         Xlo=smoothdata(Xlo,1,'gaussian',round(3/dt));
%         Xhi=smoothdata(Xhi,1,'gaussian',round(3/dt));
        
        XNorm=(X-Xlo)./abs(Xhi-Xlo);
%         XNorm(XNorm>1)=1;
%         XNorm(XNorm<0)=0;
        
    otherwise
        error(['unknown method: ' method]);
end


%plot to show result
if nargout==0 || doPlot==1
    tix=1;
    fh=figure('name','Normalize Traces','KeyPressFcn',@keypressFcn);
    XLIM = [min(t), max(t)];
    YLIM = [min(X(:)), max(X(:))];
    YLIMn = [min(XNorm(:)), max(XNorm(:))];
    plotData()
    
end


%nested functions can see variables in caller's scope
    function plotData()
        subplot(2,1,1)
        plot(t,X(:,tix))
        switch method
            case {'trend'}
                hold on
                plot(t,XTrend(:,tix))
                hold off
            case {'wptile'}
                hold on
                plot(t,Xlo(:,tix))
                plot(t,Xhi(:,tix))
                hold off
        end
        grid on
        xlabel('Time')
        %         ylabel('raw')
        xlim(XLIM)
        ylim(YLIM)
        title(tix)
        line(t,X,-1*ones(size(t)),'Color',0.7*[1,1,1])

        
        subplot(2,1,2)
        plot(t,XNorm(:,tix))
        grid on
        xlabel('Time')
        ylabel('normalized')
        xlim(XLIM)
        ylim(YLIMn)
        line(t,XNorm,-1*ones(size(t)),'Color',0.7*[1,1,1])
        
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

