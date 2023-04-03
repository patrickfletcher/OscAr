function [P1,f,Pmax,fmax]=powerSpectrumOscar(X,fs,n,doPlot)
%if x=matrix, columns are timeseries

%actually this is a periodogram

%make vector column
if any(size(X)==1)
    X=X(:);
end

if ~exist('doPlot','var')
    doPlot=0;
end

if ~exist('n','var')||isempty(n)
    n=4096;
end

L=size(X,1); %assume column is timeseries

f=fs*(0:n/2-1)/n; f=f(:);

Y=fft(X,n);

% phaseY = unwrap(angle(Y));

%Power of each freq components 
% P2=Y.*conj(Y)/(n*L);
P2=Y.*conj(Y)/(fs*L);
% P2=abs(Y).^2/(n*L);
P1=P2(1:n/2,:);
P1(2:end-1,:)=2*P1(2:end-1,:);

[Pmax,ix]=max(P1,[],1);
fmax=f(ix);


%plot to show result
if nargout==0 || doPlot==1
    
nX=size(X,2);
tix=1;
figure('KeyPressFcn',@keypressFcn);
plotData()
    
end


%nested functions can see variables in caller's scope
    function plotData()
        
        plot(f,pow2db(P1(:,tix)),'LineWidth',1); 
        % set(gca,'yscale','log')
        xlim([0,1.5])
        title('One Sided Power Spectral Density');       
        xlabel('frequency')         
        ylabel('power');
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