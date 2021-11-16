classdef Experiment < handle & matlab.mixin.Copyable
    %class to represent a single experiment containing >=1 timeseries that share the same time sample points. For
    %example, a fluorescence videomicroscopy experiment with multiple ROIs in a field of view. 

    
    %TODO: Force user to load data into Matlab? data=[t,X];
    
    %TODO: make parameterDialog a regular figure with "go" button, like
    %selectFeaturesPopup
    
    
    %TODO: add normalization to parameterDialog
    
    %TODO: get valid options for preprocessing functions exported from those
    %functions to display in dropdown lists
    
    %TODO: forced order, but flags to choose which steps (if any) to apply?
    %     - alt: arbitrary list of operations, with fully qualified options.
    %TODO: Customizeable pipeline? eg. array processingStep objects to define
    %steps (name+operation+method+params), and pointer to which one
    %gets periods measured and which one features get measured. pages
    %of X correspond to each step: eg. X(:,:,1)=raw, X(:,:,2)=norm, ..., X(:,:,end)=per
    
    
    %TODO: plotting alternatives - plot all timeseries onces, and use an
    %update function to set the visibility on/off (also would support
    %option to plot out-of-focus traces as light gray)
    
    %TODO: updatePlots needs to be axis-aware? store plottype & extra info
    %in axis userdata not figure? store axes handles instead of figure
    %handles? => single figure GUI 
    
    %TODO: plot mean+/- stdev for per-trace features, when these are distrubtions for the trace. or boxplot. or violin.
    %if so, remove the stdev from features to plot list.
    
    
    %TODO: finalize private/public properties and methods for full OOP
        
    
    %TODO: support 3D array - 3rd dim is observable id (one [nT x nX] page per observable)
    % this requires supporting different pipeline options for each
    % observable and allows computation of new higher order features,
    % such as phase relationship
    % --> actually, will likely determine period markers using only one
    % observable, then compute features relative to those markers in
    % other traces?
    
    properties 
        name=''
        date=''
        sex=''
        condition='' %eg. WT vs KO
        filename=''
        fullfilename=''
        notes={}
        
        %segment is a subinterval of the trace, the basic unit within which
        %periodicity will be detected and features will be computed
        nS
        segment=struct('name','','endpoints',[],'ix',[],...
            'points',struct(),'features_periods',struct(),'features_trace',struct())
        
        t
        dt
        fs %1/dt
        f
        psd
        psdFilt
        
        X
        nX
        includeT %place to store trace-wise include vector
        
        traceID %from original file
        groupT
        groupID
        nG
        includeG %groupwise include vector
        
        groupMode=false
        
        featureFcn %not used?
        fnames_periods
        fnames_trace
        
        resultsTrace=table
        resultsPeriod=table
        
        
        %private (?):
        N
        group
        include
        Xnorm
        Xtrend
        Xdetrend
        Xfilt
        
    end
    
    %store parameters used
    properties

%support variable length parameter lists, eg. for name-value pairs
        normMethod='none'
        normParam=[]
        normPerSegment=false
        
        trendMethod='none'
        trendParam=[]
        detrendPerSegment=false
        flattenMethod='none'
        flattenPerSegment=false
        
        filterMethod='none'
        filterParam=[]
        filterPerSegment=false
        
        averageMethod='arithmetic'
        averageParam=[]
        
        featureMethod='none'
        featureParam={}
        featureExtras={}
        
        interpMethod='linear'
        
        Tab=table;
        resultsFile='';
    end
    
    properties (SetAccess=protected)
        %interactive plots:
        fig_handles=matlab.ui.Figure.empty() %array of figure handles spawned by this object - use to synchronize trace in focus?
%         gobj_array=gobjects()
        tix %in focus
        active_fig %index into fig_handles to maintain its focus upon updates
        featurePlotType='per-period'
        xfeature='segment'
        yfeature='period'
        
        resfig=matlab.ui.Figure.empty %figure to show the results table
        featureSelectDlg
        paramDlg 
    end
    
    
    methods
        %constructor
        % Experiment(): requests file
        % Experiment(filename): looks for given file
        % Experiment(data): uses numeric matrix data=[t,X]
        % TODO: Experiment(___,Name,Value) : name-value pairs to set parameters, metadata, etc.
        function expt=Experiment(varargin)
            
            %parse inputs
            loadFile=true;
            fullfilename=[];
            filename=[];
            if ~isempty(varargin) && ( ischar(varargin{1}) || isstring(varargin{1}) )
                [path,filename,ext]=fileparts(varargin{1});
                fullfilename=[path,filesep,filename,ext];
            elseif ~isempty(varargin) && isnumeric(varargin{1}) && ~isempty(varargin{1})
                data=varargin{1};
                loadFile=false;
            end
            
            %load a new file
            if loadFile
                if isempty(fullfilename) || ~exist(fullfilename,'file')
                    [filename,path]=uigetfile({'*.xls*'});
                    fullfilename=[path,filename];
                end
                
                [num,txt,raw]=xlsread(fullfilename);
                %extract data
%                 ixs=find(cellfun(@(x)(isnumeric(x)&&~isnan(x)),raw(:,1)),1,'first');
%                 data=cell2mat(raw(ixs:end,1:size(num,2)));
%                 ixs=find(cellfun(@(x)(isnumeric(x)&&~isnan(x)),raw(:,1)),1,'first');
                data=num;
                
                %look for special metadata entries
                headerNames=txt(:,1);
                keyword={'name','date','sex','condition'};
                for i=1:length(keyword)
                    r=find(strcmpi(headerNames,keyword{i}));
                    if ~isempty(r)
                        expt.(keyword{i})=txt{r,2};
                    end
                end
                
                r=find(strcmpi(headerNames,'group'));
                if ~isempty(r)
                    expt.setGroup(cell2mat(raw(r,2:size(num,2))))
                end
                
            end
            
            expt.filename=filename;
            expt.fullfilename=fullfilename;
                
            time=data(:,1);
            time=time-time(1); %shift time to start at zero
            Xraw=data(:,2:end);
            
            DT=mode(diff(time)); 

            %interpolate any missing data using neighboring time points
            %TODO: make this optional?
            
            timeIsNan=isnan(time); %this would be a problem!
            if any(timeIsNan)
                expt.notes=[expt.notes,{'missing value in time column'}];
            end
            
            skip=any(abs(diff(time)-DT)>2*DT);
            if any(skip)
                expt.notes=[expt.notes,{'large timestep detected'}];
            end
            
            rowIsNan=any(isnan(Xraw),2);
            if any(rowIsNan)
                expt.notes=[expt.notes,{'some rows have NaN'}];
                r=find(rowIsNan);

                if r(1)==1
                    time=time(2:end,1);
                    Xraw=Xraw(2:end,:);
                    r=r(2:end)-1;
                end
                for ii=1:length(r)
                    for jj=1:size(Xraw,2)
                        iix=[r(ii)-1,r(ii)+1];
                        x=time(iix,1);
                        v=Xraw(iix,jj);
                        xq=time(r(ii),1);
                        Xraw(r(ii),jj)=interp1(x,v,xq,expt.interpMethod);
                    end
                end
            end
%             expt.t=time;
%             expt.X=Xraw;
            
            %interpolate on a fixed grid with given DT
            %TODO: DT as parameter (for same DT across expts)
            % - alt, make this a method with DT parameter
            
            ti=0:DT:time(end); ti=ti';
            XI=zeros(length(ti),size(Xraw,2));
            for j=1:size(Xraw,2)
                XI(:,j)=interp1(time,Xraw(:,j),ti,'pchip');
            end
            
            expt.t=ti;
            expt.X=XI;
            
            expt.nX=size(expt.X,2);
            expt.traceID=1:expt.nX;
            
            expt.dt=DT;
            expt.fs=1/expt.dt;
            
            %default segment = whole trace
            expt.defineSegments({'1'},[ti(1),ti(end)]);
            expt.includeT=true(1,expt.nX);
            
            expt.setGroup(ones(1,expt.nX));
            
            expt.tix=1;
             
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup helper functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function defineSegments(expt, names, startTimes, stopTimes)
            %expects names as string array, times as vector
            expt.nS=length(names);
            
            startTimes=startTimes(:);
                
            if ~exist('stopTimes','var')||isempty(stopTimes)
                startTimes(end+1)=expt.t(end); %append final time
                endpoints = [startTimes(1:expt.nS),startTimes(2:expt.nS+1)];
            else
                stopTimes=stopTimes(:);
                if length(startTimes)==length(stopTimes)+1
                    stopTimes(end+1,1)=expt.t(end);
                elseif length(startTimes)~=length(stopTimes)
                    error(['Start times and Stop times have incorrect length'])
                end
                endpoints=[startTimes,stopTimes];
            end
            
            for i=1:expt.nS
                expt.segment(i).name=names(i);
                expt.segment(i).endpoints=endpoints(i,:);
                expt.segment(i).ix=expt.t>=endpoints(i,1) & expt.t<=endpoints(i,2);
            end
        end
        
        function setGroup(expt,groupvar)
            groupvar=groupvar(:)';%row vector
            [groupvar,ixs]=sort(groupvar);
            expt.traceID=expt.traceID(ixs);
            expt.X=expt.X(:,ixs);
            expt.includeT=expt.includeT(ixs);
            
            expt.groupT=groupvar;
            expt.groupID=unique(groupvar);
            expt.nG=length(expt.groupID);
            expt.includeG=true(1,expt.nG); %could add second input to set this at same time
            
            %update the internal current group/include/N
            if expt.groupMode
                expt.group=expt.groupID;
                expt.include=expt.includeG;
                expt.N=expt.nG;
            else
                expt.group=expt.groupT;
                expt.include=expt.includeT;
                expt.N=expt.nX;
            end
        end
        
        function excludeTraces(expt,excludeIx)
            %TODO: groupMode sensitive?
            expt.includeT(ismember(expt.traceID,excludeIx))=false;
            if expt.groupMode
                expt.include=expt.includeG;
            else
                expt.include=expt.includeT;
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Preprocessing functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function setNormalizeParameters(expt,method,methodPar,normPerSegment)
            arguments
                expt Experiment
                method = 'none'
                methodPar = []
                normPerSegment = 0
            end
            expt.normMethod=method;
            expt.normParam=methodPar;
            expt.normPerSegment=normPerSegment;
        end
        
        function setDetrendParameters(expt,method,methodPar,detrendPerSegment,flattenMethod,flattenPerSegment)
            arguments
                expt Experiment
                method = 'none'
                methodPar = []
                detrendPerSegment = 0
                flattenMethod = 'none'
                flattenPerSegment = 0
            end
            expt.trendMethod=method;
            expt.trendParam=methodPar;
            expt.detrendPerSegment=detrendPerSegment;
            expt.flattenMethod=flattenMethod;
            expt.flattenPerSegment=flattenPerSegment;
        end
        
        function setFilterParameters(expt,method,methodPar, filterPerSegment)
            arguments
                expt Experiment
                method = 'none'
                methodPar = []
                filterPerSegment = 0
            end
            expt.filterMethod=method;
            expt.filterParam=methodPar;
            expt.filterPerSegment=filterPerSegment;
        end
        
        
        function preprocess(expt)
            expt.normalize();
            expt.detrend();
            if expt.groupMode
                expt.averageTraces();
            end
            expt.filter();
        end
        
        
        function normalize(expt,method,methodPar,normPerSegment)
            arguments
                expt Experiment
                method = 'none'
                methodPar = []
                normPerSegment = 0
            end
            if nargin>1
                expt.normMethod=method;
                expt.normParam=methodPar;
                expt.normPerSegment=normPerSegment;
            end
            if expt.normPerSegment
                for i=1:expt.nS
                    thisIx=expt.segment(i).ix;
                    xnorm=normalizeTraces(expt.t(thisIx),expt.X(thisIx,:),expt.normMethod,expt.normParam);
                    expt.Xnorm(thisIx,:)=xnorm;
                end
                segix=any([expt.segment(:).ix],2);
                expt.Xnorm(~segix,:)=nan;
                
            else
                expt.Xnorm=normalizeTraces(expt.t,expt.X,expt.normMethod,expt.normParam);
            end
            
        end
        
        
        %detrend supports per-segment: for piecewise linear detrend.
        function detrend(expt,method,methodPar,detrendPerSegment,flattenMethod,flattenPerSegment)
            arguments
                expt Experiment
                method = 'none'
                methodPar = []
                detrendPerSegment = 0
                flattenMethod = 'none'
                flattenPerSegment = 0
            end
            if nargin>1
                expt.trendMethod=method;
                expt.trendParam=methodPar;
                expt.detrendPerSegment=detrendPerSegment;
                expt.flattenMethod=flattenMethod;
                expt.flattenPerSegment=flattenPerSegment;
            end
            
            if isempty(expt.Xnorm)
                expt.normalize();
            end
            
            [Xdt_full,Xt_full]=detrendTraces(expt.t,expt.Xnorm,expt.trendMethod,expt.trendParam);
            full_offset = getFlattenOffset(expt.Xnorm, Xt_full, expt.flattenMethod);
            if expt.detrendPerSegment
                last=0;
                for i=1:expt.nS
                    thisIx=expt.segment(i).ix;
                    [Xdt,Xt]=detrendTraces(expt.t(thisIx),expt.Xnorm(thisIx,:),expt.trendMethod,expt.trendParam);
                    Xt_full(thisIx,:)=Xt;
                    
                    % flatten logic: all branches set seg_offset.
                    if expt.flattenPerSegment
                        % Trace will be discontinuous
                        seg_offset = getFlattenOffset(expt.Xnorm(thisIx,:), Xt, expt.flattenMethod);
                        
                    else
                        %use full trace to get flatten value: only
                        %norm-based flatten supported. Trace will be continuous
                        if i==1
                            seg_offset=full_offset-Xdt(1,:);
                        else
                            seg_offset=last-Xdt(1,:);
                        end
                    end
                    Xdt = Xdt + seg_offset;
                    
                    last=Xdt(end,:);
                    Xdt_full(thisIx,:)=Xdt;
                end
                
                segix=any([expt.segment(:).ix],2);
                Xt_full(~segix,:)=nan;
                Xdt_full(~segix,:)=nan;
            else
                % flatten as one segment (expt.flattenPerSegment=1 has no effect) 
                % Trace will be continuous
                Xdt_full = Xdt_full + full_offset;
            end
            expt.Xtrend=Xt_full;
            expt.Xdetrend=Xdt_full;
            
            
        end
        
        
        function averageTraces(expt) 
            %TODO: handle edge case - all traces for a group excluded 
            
            XDTg=zeros(length(expt.t),expt.nG);
            for j=1:expt.nG  
                XDTg(:,j)=mean(expt.Xdetrend(:,expt.groupT==expt.groupID(j) & expt.includeT),2);  
            end
%             emptygroup=isnan(XDTg(1,:));
%             XDTg(:,emptygroup)=[]; %remove any groups with no members
%             expt.groupID(emptygroup)=[];
            
            expt.Xdetrend=XDTg; %groupmode overwrites detrend... TODO: keep both for plotting

%             segix=any([expt.segment(:).ix],2);
%             expt.Xdetrend(~segix,:)=nan;
        end
        
        
        function filter(expt, method, methodPar, filterPerSegment)
            arguments
                expt Experiment
                method = 'none'
                methodPar = []
                filterPerSegment = 0
            end
            if nargin>1
                expt.filterMethod=method;
                expt.filterParam=methodPar;
                expt.filterPerSegment=filterPerSegment;
            end
            
            if isempty(expt.Xnorm)
                expt.normalize();
            end
            if isempty(expt.Xdetrend)
                expt.detrend();
            end
            
            if expt.filterPerSegment
                for i=1:expt.nS
                    thisIx=expt.segment(i).ix;
                    xfilt=filterTraces(expt.t(thisIx),expt.Xdetrend(thisIx,:),expt.filterMethod,expt.filterParam);
                    expt.Xfilt(thisIx,:)=xfilt;
                end
                segix=any([expt.segment(:).ix],2);
                expt.Xfilt(~segix,:)=nan;
            else
                expt.Xfilt=filterTraces(expt.t,expt.Xdetrend,expt.filterMethod,expt.filterParam);
            end
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Analysis functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TODO expand this - detect_periods return only points/fcns, then compute_features should apply fcn to relevant
        %traces.
        
        function setFeatureParameters(expt,method,methodpar,varargin)
            expt.featureMethod=method;
            expt.featureParam=methodpar;
            expt.featureExtras=varargin;
        end
        
        function compute_features(expt,method,methodpar,varargin)
            if nargin>1
            expt.featureMethod=method;
            expt.featureParam=methodpar;
            expt.featureExtras=varargin;
            end
            if ~strcmp(expt.featureMethod,'none')
            expt.detect_periods();
            expt.periodogram();
            expt.fnames_periods=fieldnames(expt.segment(1).features_periods)';
            expt.fnames_trace=fieldnames(expt.segment(1).features_trace)';
            expt.buildResultsTable();
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function plotTrace(expt,whichPlot,tix,showPts,doInteractive)
            %dispatch function for plotting expt data
            
            %TODO: support cell array for whichplot
                   
            if exist('tix','var')&&~isempty(tix)
                expt.tix=tix;
            end
            
            if ~exist('showPts','var')||isempty(showPts)
                showPts=true;
            end
            
            if ~exist('doInteractive','var')||isempty(doInteractive)
                doInteractive=true;
            end
            
            if ~iscell(whichPlot)
                whichPlot={whichPlot};
            end
            
            %interactive figure: set up callbacks
            if doInteractive
                figID=gcf;
                figID.Name=strcat('Traces: ',expt.name);
                figID.NumberTitle='off';
                figID.KeyPressFcn=@expt.commonKeypress;
%                 figID.Interruptible='off';
                figID.UserData={'trace',whichPlot,showPts};
                figID.CloseRequestFcn=@expt.figureCloseFcn;
                if ~ismember(figID,expt.fig_handles)
                    expt.fig_handles(end+1)=figID; %register the new fig with expt
                end
            end
            
            expt.active_fig=gcf;
            nPlots=length(whichPlot);
%             tiledlayout(expt.active_fig, nPlots,1,'TileSpacing','compact','Padding','compact');
            
            %set up the objects
            for j=1:nPlots
%                 ax(j)=nexttile(j);
                ax(j)=subplot(nPlots,1,j);
%                 uistack(ax(j),'bottom')
                ax(j).XTick=[];
%                 ax(j).YLabel.String=whichPlot;
                
                for i=1:expt.nX
                    line(ax(j),nan,nan,'color',0.85*[1,1,1],'tag','xall_line');
                end
                line(ax(j),nan,nan,'color',0.1*[1,1,1],'tag','x_line');
                line(ax(j),nan,nan,'color',[0.3,0.7,1],'tag','x2_line')
            
                for i=1:expt.nS
                    line(ax(j),nan,nan,'color','r','marker','^','linestyle','none','tag','min_line')
                    line(ax(j),nan,nan,'color','r','marker','v','linestyle','none','tag','max_line')
                    line(ax(j),nan,nan,'color','g','marker','<','linestyle','none','tag','dxmin_line')
                    line(ax(j),nan,nan,'color','g','marker','>','linestyle','none','tag','dxmax_line')
                    line(ax(j),nan,nan,'color','b','marker','d','linestyle','none','tag','up_line')
                    line(ax(j),nan,nan,'color','b','marker','o','linestyle','none','tag','down_line')
                    line(ax(j),nan,nan,'color','b','marker','s','linestyle','none','tag','period_line')
                end
                for i=1:expt.nS-1
                    line(ax(j),nan,nan,'color',[0,0.75,0],'tag','seg_start_line')
                    line(ax(j),nan,nan,'color',[0.75,0,0],'tag','seg_stop_line')
                end
            end
            ax(end).XTickMode='auto';
            ax(end).XLabel.String='t';
            linkaxes(ax,'x')
            
%             expt.active_fig.Children=expt.active_fig.Children(nPlots:-1:1);
            
            expt.updatePlots('trace')
            
%             uistack(
        end
        
        function plotPeriodogram(expt,tix,doInteractive)
            %overlay psd for each segment
            if ~isempty(expt.psd)
            
            if exist('tix','var')&&~isempty(tix)
                expt.tix=tix;
            end
            
            if ~exist('doInteractive','var')||isempty(doInteractive)
                doInteractive=true;
            end
                      
            if doInteractive
                figID=gcf;
                figID.Name=strcat('Periodogram: ',expt.filename);
                figID.NumberTitle='off';
                if ~ismember(figID,expt.fig_handles)
                    expt.fig_handles(end+1)=figID; %register the new fig with expt
                end
                figID.KeyPressFcn=@expt.commonKeypress;
%                 figID.Interruptible='off';
                figID.UserData={'psd'};
                figID.CloseRequestFcn=@expt.figureCloseFcn;
            end
            
            expt.plot_psd();
            expt.active_fig=gcf;
            end
        end
        
        function plotFeatures(expt,xname,yname,tix,doInteractive)
            %TODO: how to manage plotting options?
                
%             if ~isempty(expt.segment(1).features)
            
            if exist('tix','var')&&~isempty(tix)
                expt.tix=tix;
            end
            
            if ~exist('doInteractive','var')||isempty(doInteractive)
                doInteractive=true;
            end
            
            %this part is a kludge
            if exist('xname','var')&&~isempty(xname)
                xname=validatestring(xname,['t','segment',expt.fnames_periods,expt.fnames_trace]);
            else
                xname=expt.xfeature;
            end
            if exist('yname','var')&&~isempty(yname)
                yname=validatestring(yname,[expt.fnames_periods,expt.fnames_trace]);
            else
                yname=expt.yfeature;
            end
            if xname=="segment"
                if ismember(yname,expt.fnames_periods)
                    featPlotType='per-period';
                else
                    featPlotType='per-trace';
                end
            else
                xp=ismember(xname,['t',expt.fnames_periods]); %valid for period, otherwise trace
                yp=ismember(yname,expt.fnames_periods);
                if xp&&yp
                    featPlotType='per-period';
                elseif ~xp&&~yp
                    featPlotType='per-trace';
                else
                    error('Incompatible feature types: mixed trace/period')
                end
            end
            
%             expt.xfeature=xname;
%             expt.yfeature=yname;
            
            %interactive figure: set up callbacks
            if doInteractive
                figID=gcf;
                figID.Name=strcat('Features: ',expt.filename);
                figID.NumberTitle='off';
                figID.KeyPressFcn=@expt.commonKeypress;
%                 figID.Interruptible='off';
                figID.UserData={'feature',xname,yname,featPlotType};
                figID.CloseRequestFcn=@expt.figureCloseFcn;
                if ~ismember(figID,expt.fig_handles)
                    expt.fig_handles(end+1)=figID; %register the new fig with expt
                end
            end
            
            switch featPlotType
                case 'per-period'
                    expt.plot_features_periods();
                case 'per-trace'
                    expt.plot_features_trace();
            end
            expt.active_fig=gcf;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % saving/displaying results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function displayResults(expt,createNewFig)
            arguments
                expt Experiment
                createNewFig=true 
            end
            
            %simple table view of results to cross reference with traces
            %examined by keypress
            %
            % keypress here to switch segment?
            %  TODO: include=true/false column, select traces for result
            %  export?
            
            
            if isempty(expt.resultsTrace)
%                 warning('results table is empty, nothing to do')
                return
            end

            figID=[];
            if ~createNewFig %check if fig exists
                for i=1:length(expt.fig_handles)
                    if expt.fig_handles(i).UserData{1}=="results"
                        figID=expt.fig_handles(i);
                        break
                    end
                end
            end
            if isempty(figID)
                figID=gcf; %uses current figure, or creates one if no figs
                figID.Name=strcat('Result Table: ',expt.filename);
                figID.NumberTitle='off';
                figID.KeyPressFcn=@expt.commonKeypress;
                figID.UserData={'results'};
    %             figID.Interruptible='off';
                figID.CloseRequestFcn=@expt.figureCloseFcn;
                if ~ismember(figID,expt.fig_handles)
                    expt.fig_handles(end+1)=figID; %register the new fig with expt
                end
            end
            figure(figID);
            
            colnames=[{'Trace'};{'Group'};{'Segment'};{'Include'};fieldnames(expt.segment(1).features_trace)];
%             colnames=[{'Trace'};{'Segment'};{'Include'};fieldnames(expt.segment(1).features_trace)];
%             colformat=repmat({'numeric'},1,length(colnames));
%             coledit=false(1,length(colnames));
%             colformat(3)={'logical'}; 
%             coledit(3)=true; 
            %needs a callback function to link edited values to expt.include
            
            uit=uitable(figID,'Data',expt.resultsTrace{:,:},'tag','oscar_resultsTable');
            uit.ColumnName=colnames;
%             uit.ColumnFormat=colformat;
%             uit.ColumnEditable=coledit;
            uit.RowName=[];
            uit.Units='normalized';
            uit.Position=[0.025,0.025,0.95,0.95];
            
            figure(figID) %bring focus back to figure (for keypressfcn)
            figure(expt.active_fig) %bring focus back to active figure
        end
        
%         function displayResultsKeypress
%         end

        function destinationExperiment=copyProcessedData(expt, destinationExperiment)
            
            destinationExperiment.normMethod=expt.normMethod;
            destinationExperiment.normParam=expt.normParam;
            destinationExperiment.normPerSegment=expt.normPerSegment;
            destinationExperiment.trendMethod=expt.trendMethod;
            destinationExperiment.trendParam=expt.trendParam;
            destinationExperiment.detrendPerSegment=expt.detrendPerSegment;
            destinationExperiment.flattenMethod=expt.flattenMethod;
            destinationExperiment.flattenPerSegment=expt.flattenPerSegment;
            destinationExperiment.filterMethod=expt.filterMethod;
            destinationExperiment.filterParam=expt.filterParam;
            destinationExperiment.filterPerSegment=expt.filterPerSegment;
            destinationExperiment.featureMethod=expt.featureMethod;
            destinationExperiment.featureParam=expt.featureParam;
            destinationExperiment.featureExtras=expt.featureExtras;
            
            destinationExperiment.averageMethod=expt.averageMethod;
            destinationExperiment.averageParam=expt.averageParam;
            
            destinationExperiment.Xnorm=expt.Xnorm;
            destinationExperiment.Xtrend=expt.Xtrend;
            destinationExperiment.Xdetrend=expt.Xdetrend;
            destinationExperiment.Xfilt=expt.Xfilt;
            destinationExperiment.psd=expt.psd;
            destinationExperiment.psdFilt=expt.psdFilt;
            destinationExperiment.segment=expt.segment;
            destinationExperiment.fnames_periods=expt.fnames_periods;
            destinationExperiment.fnames_trace=expt.fnames_trace;
            destinationExperiment.resultsPeriod=expt.resultsPeriod;
            destinationExperiment.resultsTrace=expt.resultsTrace;
            
            destinationExperiment.include=expt.include;
            destinationExperiment.includeT=expt.includeT;
            destinationExperiment.includeG=expt.includeG;
            destinationExperiment.group=expt.group;
            destinationExperiment.groupMode=expt.groupMode;
        end

        function setSaveArgs(expt, Tab, resultsFile)
            expt.Tab=Tab;
            expt.resultsFile=resultsFile;
        end
        
        function saveResults(expt)
            normPar=expt.normParam; if isempty(normPar), normPar="none"; end
            trendPar=expt.trendParam; if isempty(trendPar), trendPar="none"; end
            filterPar=expt.filterParam; if isempty(filterPar), filterPar="none"; end
            
            thisInfo=expt.Tab;
            thisInfo.normMethod=repmat(string(expt.normMethod),expt.nS*expt.nX,1);
            thisInfo.normParam=repmat(normPar,expt.nS*expt.nX,1);
            thisInfo.normPerSegment=repmat(expt.normPerSegment,expt.nS*expt.nX,1);
            thisInfo.trendMethod=repmat(string(expt.trendMethod),expt.nS*expt.nX,1);
            thisInfo.trendParam=repmat(trendPar,expt.nS*expt.nX,1);
            thisInfo.detrendPerSegment=repmat(expt.detrendPerSegment,expt.nS*expt.nX,1);
            thisInfo.filterMethod=repmat(string(expt.filterMethod),expt.nS*expt.nX,1);
            thisInfo.filterParam=repmat(filterPar,expt.nS*expt.nX,1);
            thisInfo.filterPerSegment=repmat(expt.filterPerSegment,expt.nS*expt.nX,1);
            thisInfo.featureMethod=repmat(string(expt.featureMethod),expt.nS*expt.nX,1);
            thisInfo.featureParam=repmat(expt.featureParam,expt.nS*expt.nX,1);
            thisInfo.minAmp=repmat(expt.featureExtras{1}.MinimumAmplitude,expt.nS*expt.nX,1);

            thisResult=expt.resultsTrace(:,3:end);

            segNames=[expt.segment(:).name]';
            thisSegNames=segNames(thisResult.Segment);
            thisInfo.SegmentName=thisSegNames;

            thisInfo=[thisInfo,thisResult];

            if ~exist(expt.resultsFile,'file')
                resultsTab=thisInfo; %first time, make a new table
            else
                resultsTab=load(expt.resultsFile);  %need a way to give this var a name???
                FN=fieldnames(resultsTab);
                resultsTab=resultsTab.(FN{1});
                rows=ismember(resultsTab(:,1:2),thisInfo(:,1:2),'rows');
                if any(rows) 
                    resultsTab(rows,:)=thisInfo; %replace existing rows
                else
                    resultsTab=[resultsTab;thisInfo]; %append new rows
                end
            end
            save(expt.resultsFile,'resultsTab')
            disp("Saved expt.resultsTrace to " + expt.resultsFile)
        end

        
        function writeToExcel(expt,outfilename,doDistributions)
            %write a header region with file/experiment info
            %default: trace level (means/stdevs) features
            
            %trace features: column 1=traceID, col2=segment, then features
            
            %option to instead do period-wise distributions - how to handle
            %jagged arrays? one segment per sheet?
            
            if isempty(expt.resultsTrace)
                warning('results table is empty, nothing to do')
                return
            end
            
            
            %header?
            
            %auto-filename...?
            if ~exist('outfilename','var')||isempty(outfilename)
                [path,file,fext]=fileparts(expt.fullfilename);
                outfilename=fullfile(path,[file,'_results',fext]);
                
                %suggest file name and request confirmation:
                filter='*.xls*';
                title='Enter the output file name';
                [file,path]  = uiputfile(filter,title,outfilename);
                if isequal(file,0) || isequal(path,0)
                    return
                else
                    outfilename=fullfile(path,file);
                end
            end
            
            writetable(expt.resultsTrace,outfilename);
        end
        
        function buildResultsTable(expt)
            ID=repmat((1:expt.N)',expt.nS,1);
            INC=repmat(expt.include',expt.nS,1);
            GRP=repmat(expt.group',expt.nS,1);
            SEG=[];TAB=table;
            for i=1:expt.nS
                SEG=[SEG;i*ones(expt.N,1)];
                TAB=[TAB;struct2table(expt.segment(i).features_trace)];
            end
            expt.resultsTrace=table();
            expt.resultsTrace.ID=ID;
            expt.resultsTrace.Group=GRP;
            expt.resultsTrace.Segment=SEG;
            expt.resultsTrace.Include=INC;
            expt.resultsTrace=[expt.resultsTrace,TAB];
        end
        
        
%         function set.tix(expt,newval)
%             expt.tix=newval;
%             expt.updatePlots();
%         end
        
        function clearFigs(expt)
            %need better way to keep track of all figs registered?
            expt.resfig=matlab.ui.Figure.empty;
            expt.active_fig=matlab.ui.Figure.empty;
            expt.fig_handles=matlab.ui.Figure.empty;
            expt.featureSelectDlg=matlab.ui.Figure.empty;
        end

    end
    
    methods (Access=private)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Analysis functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function detect_periods(expt)
%             expt.featureMethod=method;
%             expt.featureParam=methodpar;
%             expt.featureExtras=varargin{:};
            
            %TODO: per-segment toggle

            extras=expt.featureExtras;
            for i=1:expt.nS
                thisIx=expt.segment(i).ix;
                
%                 ix=find(thisIx,5,'first'); %exclude k first time points
%                 thisIx(ix)=false;
                
                tt=expt.t(thisIx);
                xx=expt.Xfilt(thisIx,:);
                switch expt.featureMethod
                    case {'peaks'}
                        delta=expt.featureParam;
%                         if length(varargin)>1, extras=varargin(2:end); end
                        [F, Fdist, points, fcns]=peak_detector(tt,xx,delta,extras{:});

                    case {'threshold'}
                        frac=expt.featureParam;
%                         if length(varargin)>1, extras=varargin(2:end); end
                        [F, Fdist, points, fcns]=plateau_detector(tt,xx,frac,extras{:});
                end
                
                expt.segment(i).points=points;
                expt.segment(i).features_periods=Fdist;
                expt.segment(i).features_trace=F;
            end
            expt.featureFcn=fcns.compute_features;
        end
        
        function periodogram(expt)
            
            nsamp=2048;
            peakdelta=0.05;
            
            expt.psd=zeros([nsamp/2,expt.N,expt.nS]);
            expt.psdFilt=zeros([nsamp/2,expt.N,expt.nS]);
            
            for i=1:expt.nS
                thisIx=expt.segment(i).ix;
                DT=expt.Xdetrend(thisIx,:);
                DT=DT-mean(DT,1);
                [PSD1,F,Pmax1,fmax1]=powerSpectrumOscar(DT, expt.fs, nsamp);
                Tpsd1=1./fmax1;
                expt.psd(:,:,i)=PSD1;
                
                DTF=expt.Xfilt(thisIx,:);
                DTF=DTF-mean(DTF,1);
                [PSD2,~,Pmax2,fmax2]=powerSpectrumOscar(DTF,expt.fs, nsamp);
                Tpsd2=1./fmax2;
                expt.psdFilt(:,:,i)=PSD2;
                
                PSD=PSD1;
                Pmax=Pmax1;
                fmax=fmax1;
                Tpsd=Tpsd1;
%                 PSD=PSD2;
%                 Pmax=Pmax2;
%                 fmax=fmax2;
%                 Tpsd=Tpsd2;
                
                PSD(1,:)=Pmax;
                [~,~,points]=peak_detector(F, PSD, peakdelta);
                m2=zeros(1,length(points));
                for j=1:length(points)
                    othermax=points(j).max.x~=Pmax(j);
                    if any(othermax)
                        m2(j)=max(points(j).max.x(othermax));
                    end
                end
                Rp21=m2./Pmax; %second peak power to first peak power ratio
                
                Tpsd=num2cell(Tpsd);
                fmax=num2cell(fmax);
                Pmax=num2cell(Pmax);
                Rp21=num2cell(Rp21);
                
                [expt.segment(i).features_trace.Tpsd]=Tpsd{:};
                [expt.segment(i).features_trace.fmax]=fmax{:};
                [expt.segment(i).features_trace.Pmax]=Pmax{:};
                [expt.segment(i).features_trace.Rp21]=Rp21{:};
                
            end
            expt.f=F; %frquency vector
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function plot_t(expt,whichPlot,showPts)
            
            ax=gca;
%             reset(ax)
%             axes(ax)
%             delete(findobj(ax,'tag','oscar_line'));
                
            points_option='period';
            x2=[];
            thisG=expt.group(expt.tix);
            ix2plot=expt.group==thisG;
            gixlist=find(ix2plot);
            gix=find(gixlist==expt.tix);
            
            segix=any([expt.segment(:).ix],2);
            segstart=arrayfun(@(x) find(x.ix,1,'first'),expt.segment);
            segstart(1)=[];
            tt=expt.t;
            switch whichPlot
                case {'r','raw'}
%                     if ~expt.groupMode
                    xall=expt.X(:,ix2plot);
                    x=expt.X(:,expt.tix);
%                     end
                    
                case {'n','norm','normalize','normalized'}
                    xall=expt.Xnorm(:,ix2plot);
                    x=expt.Xnorm(:,expt.tix);
                    if expt.trendMethod~="none"
                        x2=expt.Xtrend(:,expt.tix);
                    end
%                     if ~expt.groupMode
%                     x=expt.Xnorm(:,ix2plot);
%                     x2=expt.Xtrend(:,expt.tix);
%                     end
                    
                case {'d','detrend','detrended'}
                    xall=expt.Xdetrend(:,ix2plot);
                    x=expt.Xdetrend(:,expt.tix);
                    if expt.filterMethod~="none"
                        x2=expt.Xfilt(:,expt.tix);
                    end
                    
                case {'f','filt','filter','filtered'}
                    xall=expt.Xfilt(:,ix2plot);
                    x=expt.Xfilt(:,expt.tix);
                    points_option='all';
                    
                otherwise
                    error([whichPlot, ' is not a supported trace type for plotting']);
            end

%                 hxall=line(tt,xall,zeros(size(tt)),'color',0.85*[1,1,1],'tag','xall_line');
%                 hx=line(tt(segix),x(segix),ones(size(tt(segix))),'color',0.1*[1,1,1],'tag','x_line');
%                 hx2=line(tt(segix),x2(segix),2*ones(size(tt(segix))),'color',[0.3,0.7,1],'tag','x2_line')
%                 line(pts.min.t,pts.min.x,'color','r','marker','^','linestyle','none','tag','min_line')
%                 line(pts.max.t,pts.max.x,'color','r','marker','v','linestyle','none','tag','max_line')
%                 line(pts.dxmin.t,pts.dxmin.x,'color','g','marker','<','linestyle','none','tag','dxmin_line')
%                 line(pts.dxmax.t,pts.dxmax.x,'color','g','marker','>','linestyle','none','tag','dxmax_line')
%                 line(pts.up.t,pts.up.x,'color','b','marker','d','linestyle','none','tag','up_line')
%                 line(pts.down.t,pts.down.x,'color','b','marker','o','linestyle','none','tag','down_line')
%                 line(pts.period.t,xPer,'color','b','marker','s','linestyle','none','tag','period_line')

            if ~isempty(x)
                
                xall(segstart,:)=nan;
                x(segstart)=nan;
                
                hxall=findobj(ax.Children,'tag','xall_line');
                for i=1:size(xall,2);
                    set(hxall(i),'XData',tt,'YData',xall(:,i))
                end
                
                hx=findobj(ax.Children,'tag','x_line');
                set(hx,'XData',tt(segix),'YData',x(segix))
                
%                 hxall=line(tt,xall,zeros(size(tt)),'color',0.85*[1,1,1],'tag','xall_line');
%                 hx=line(tt(segix),x(segix),ones(size(tt(segix))),'color',0.1*[1,1,1],'tag','x_line');

                if ~isempty(x2) %&& ~all(x2==0)
                    x2(segstart)=nan;
                    hx2=findobj(ax.Children,'tag','x2_line');
                    set(hx2,'XData',tt(segix),'YData',x2(segix))
                    
%                     hx2=line(tt(segix),x2(segix),2*ones(size(tt(segix))),'color',[0.3,0.7,1],'tag','x2_line')
                end
                
                    
                %slower, but benefit of encapsulating special point plots in
                %detector code
    %             plotTrace=false;
    %             for i=1:expt.nS
    %                 tt=expt.t(expt.segment(i).ix);
    %                 xx=x(expt.segment(i).ix);
    %                 pts=expt.segment(i).points(expt.tix);
    %                 pts.period.x=interp1(tt,xx,pts.period.t); %interp x value of period markers
    %                 
    %                 expt.segment(i).plot(tt,xx,pts,1,plotTrace,points_option);
    %             end

                if showPts
                    for i=1:expt.nS
                        tt=expt.t(expt.segment(i).ix);
                        xx=x(expt.segment(i).ix);
                        if isfield(expt.segment(i), 'points') && isfield(expt.segment(i).points, 'period')
                            
                            pts=expt.segment(i).points(expt.tix);
                            hmin=findobj(ax.Children,'tag','min_line');
                            hmax=findobj(ax.Children,'tag','max_line');
                            hdxmin=findobj(ax.Children,'tag','dxmin_line');
                            hdxmax=findobj(ax.Children,'tag','dxmax_line');
                            hup=findobj(ax.Children,'tag','up_line');
                            hdown=findobj(ax.Children,'tag','down_line');
                            hperiod=findobj(ax.Children,'tag','period_line');
                            if length(pts.period.t)>1
                                if points_option=="all"
                                    set(hmin(i),'XData',pts.min.t,'YData',pts.min.x)
                                    set(hmax(i),'XData',pts.max.t,'YData',pts.max.x)
                                    set(hdxmin(i),'XData',pts.dxmin.t,'YData',pts.dxmin.x)
                                    set(hdxmax(i),'XData',pts.dxmax.t,'YData',pts.dxmax.x)
                                    set(hup(i),'XData',pts.up.t,'YData',pts.up.x)
                                    set(hdown(i),'XData',pts.down.t,'YData',pts.down.x)
                                    
%                                     line(pts.min.t,pts.min.x,'color','r','marker','^','linestyle','none','tag','min_line')
%                                     line(pts.max.t,pts.max.x,'color','r','marker','v','linestyle','none','tag','max_line')
%                                     line(pts.dxmin.t,pts.dxmin.x,'color','g','marker','<','linestyle','none','tag','dxmin_line')
%                                     line(pts.dxmax.t,pts.dxmax.x,'color','g','marker','>','linestyle','none','tag','dxmax_line')
%                                     line(pts.up.t,pts.up.x,'color','b','marker','d','linestyle','none','tag','up_line')
%                                     line(pts.down.t,pts.down.x,'color','b','marker','o','linestyle','none','tag','down_line')
                                    
            %                         tupdwn=[pts.up.t(1:length(pts.down.t)); pts.down.t];
            %                         ythresh=[1;1]*expt.segment(i).features(expt.tix).pthresh;
            %                         plot(ax,tupdwn,ythresh,'b-')
                                else
                                    xPer=interp1(tt,xx,pts.period.t); %this is necessary except for filt
                                    set(hperiod(i),'XData',pts.period.t,'YData',xPer)
                                    
%                                     line(pts.period.t,xPer,'color','b','marker','s','linestyle','none','tag','period_line')
                                    
                                end
                            else
                                set(hmin(i),'XData',nan,'YData',nan)
                                set(hmax(i),'XData',nan,'YData',nan)
                                set(hdxmin(i),'XData',nan,'YData',nan)
                                set(hdxmax(i),'XData',nan,'YData',nan)
                                set(hup(i),'XData',nan,'YData',nan)
                                set(hdown(i),'XData',nan,'YData',nan)
                                set(hperiod(i),'XData',nan,'YData',nan)
                            end
                        end
                    end
                end
%                 
                axis(ax,'tight')
                ax.YLabel.String="X_{"+whichPlot+"}"+num2str(expt.tix);
                
                endpoints=[expt.segment(:).endpoints];
                starts=endpoints(1:2:end);
                stops=endpoints(2:2:end);
                hseg_start=findobj(ax.Children,'tag','seg_start_line');
                hseg_stop=findobj(ax.Children,'tag','seg_stop_line');
                for i=1:expt.nS-1
%                     set(hseg_start(i),'YData',ylim(ax))
%                     set(hseg_stop(i),'YData',ylim(ax))
                    set(hseg_start(i),'XData',[1,1]*starts(i+1),'YData',ylim(ax))
                    set(hseg_stop(i),'XData',[1,1]*stops(i),'YData',ylim(ax))
                end
%                     line(ax(j),[1,1]*starts(i+1),ylim,'color',[0,0.75,0],'tag','seg_start_line')
%                     line(ax(j),[1,1]*stops(i),ylim,'color',[0.75,0,0],'tag','seg_stop_line')
                
                
            end
        end
        
        function plot_psd(expt)
            %BUG - per-axis isn't resetting when more than one segment
%             ax=findobj(gcf,'type','axes');
%             ax=flipud(ax);
%             reset(ax)
            delete(findobj(gcf,'tag','oscar_line'));
            
            for i=1:expt.nS
                ax(i)=subplot(expt.nS,1,i);
%                 axes(ax(i))
%                 psd=expt.psd(:,expt.tix,i);
                psd=[expt.psd(:,expt.tix,i),expt.psdFilt(:,expt.tix,i)];
                line(ax(i),expt.f,psd,'tag','oscar_line'); 
                title(ax(i),expt.segment(i).name)
                axis(ax(i),'tight')
            end     
            xlabel('frequency')
            ylabel('power');
            legend(ax(end),{'detrended','filtered'},'AutoUpdate','off')
            
            maxFilt=max([expt.segment(1).features_trace(expt.tix).fmax, ...
                         expt.segment(2).features_trace(expt.tix).fmax]);

            for i=1:expt.nS
                x=expt.segment(i).features_trace(expt.tix).fmax;
                y=expt.segment(i).features_trace(expt.tix).Pmax;
                line(ax(i),x,y,'color','r','marker','v','tag','oscar_line'); 
                xlim(ax(i),[0,maxFilt*3]);
            end     
        end
        
        function plot_features_periods(expt)
            
            ax=gca;
            reset(ax)
            delete(findobj(ax,'tag','oscar_line'));
            
            fig=gcf;
            xfeat=fig.UserData{2};
            yfeat=fig.UserData{3};
            
            xJitter=0.05; %make this a property of the class?
           
            XX=cell(1,expt.nS);
            YY=cell(1,expt.nS);
            for i=1:expt.nS
                        
                yall=[];
                yallix=[];
                includeix=[];
                for e=1:expt.nX
                    thise=expt.segment(i).features_periods(e).(yfeat);
                    yall=[yall,thise];
                    yallix=[yallix,e*ones(size(thise))];
                    includeix=[includeix,expt.include(e)*ones(size(thise))];
                end
                INCIX{i}=logical(includeix);
%                 y=expt.segment(i).features_periods(expt.tix).(yfeat);
                y=yall(yallix==expt.tix);
                if isempty(y)
                    return
                end
                
                switch expt.xfeature
                    case {'t'}
                        %TODO: option for marker location?
                         xall=([expt.segment(i).points(:).period.t(1:end-1)]+[expt.segment(i).points(:).period.t(2:end)])/2;
%                          x=xall(yallix==expt.tix);
%                          x=(expt.segment(i).points(expt.tix).period.t(1:end-1)+expt.segment(i).points(expt.tix).period.t(2:end))/2; %midpoint of period
%                         x=expt.segment(i).points(expt.tix).max.t;
%                         x=expt.segment(i).points(expt.tix).min.t(1:end-1); %first min
%                         x=expt.segment(i).points(expt.tix).min.t(2:end); %second min
%                         x=expt.segment(i).points(expt.tix).down.t;
%                         x=expt.segment(i).points(expt.tix).up.t;
%                         x=(expt.segment(i).points(expt.tix).down.t+expt.segment(i).points(expt.tix).period.t(2:end))/2; %midpoint of silent post active
%                         x=(expt.segment(i).points(expt.tix).up.t+expt.segment(i).points(expt.tix).period.t(1:end-1))/2; %midpoint of silent pre active   
                        
                    case {'segment'}
                        %TODO: option for violin plot, boxplot?
                        rng(42)
                        xall=i+xJitter*randn(1,length(yall));
%                         x=i+xJitter*randn(1,length(y));
%                         colorBySegment=true;

                    case expt.fnames_periods
                        xall=[expt.segment(i).features_periods(:).(xfeat)];
%                         x=expt.segment(i).features_periods(expt.tix).(xfeat);
%                         colorBySegment=true;

                    otherwise
                        error(['Invalid name for x-axis: ' xfeat])
                end
                x=xall(yallix==expt.tix);
                
                if isempty(x)
                    y=[];
                end
                
                XX{i}=x;
                YY{i}=y;
                XXall{i}=xall;
                YYall{i}=yall;

            end
                          
            
            hl=matlab.graphics.chart.primitive.Line.empty(expt.nS,0);
            for i=1:expt.nS
                col=ax.ColorOrder(mod(i-1,length(ax.ColorOrder(:,1)))+1,:);
                if ~isempty(XXall{i})
                    line(XXall{i}(INCIX{i}),YYall{i}(INCIX{i}),'marker','o','linestyle','none',...
                        'color',0.75*[1,1,1],'tag','oscar_line');
                    line(XXall{i}(~INCIX{i}),YYall{i}(~INCIX{i}),'marker','x','linestyle','none',...
                        'color',0.75*[1,1,1],'tag','oscar_line');
                end
                if ~isempty(XX{i})
                    if (expt.include(expt.tix))
                    hl(i,1)=line(XX{i},YY{i},'marker','o','linestyle','none',...
                        'color',col,'tag','oscar_line');
                    else
                    hl(i,1)=line(XX{i},YY{i},'marker','x','linestyle','none',...
                        'color',col,'tag','oscar_line');
                    end
                end
            end
%             colorBySegment=false;  %TODO: always color by segment?
%             if colorBySegment
%                 legend(hl,[expt.segment(:).name],'AutoUpdate','off') %TODO: option to suppress legend?
%             else
%                 legend('off')
%                 for i=1:expt.nS
%                     if ~isempty(XX{i})
%                         hl(i).Color='k';
%                     end
%                 end
%             end
            
            axis tight
            
            switch xfeat
                    case {'t'}
                        for i=1:expt.nS
                            if strcmpi(class(hl(i)),'matlab.graphics.chart.primitive.Line')
                            hl(i).LineStyle='-';
                            end
                        end
                        hold on
                        for i=2:expt.nS
                            line(expt.segment(i).endpoints(1)*[1,1],ylim(),'color','g','tag','oscar_line')
                        end
                        hold off
                        xlim([expt.t(1),expt.t(end)])
                        
                    case {'segment'} %TODO: do jitter here?
                        xticks(1:expt.nS)
                        xticklabels({expt.segment.name})
                        xlim([0.5,expt.nS+0.5])
            end
            
            xlabel(xfeat)
            ylabel(yfeat)
            
        end
        
        function plot_features_trace(expt)
            
            ax=gca;
            reset(ax)
            delete(findobj(ax,'tag','oscar_line'));
            
            fig=gcf;
            xfeat=fig.UserData{2};
            yfeat=fig.UserData{3};
            
            HX=[];
            HY=[];
            XX=[];
            YY=[];
            
            for i=1:expt.nS
                
                switch xfeat
                        
                    case {'segment'}
                        %plot all trace data as line seg, highlight point expt.tix
                        xx=i*ones(1,expt.N);
                        yy=[expt.segment(i).features_trace.(yfeat)];
                        HX(end+1,1)=i;
                        HY(end+1,1)=yy(expt.tix);
                        

                    case expt.fnames_trace
                        xx=[expt.segment(i).features_trace.(xfeat)];
                        yy=[expt.segment(i).features_trace.(yfeat)];
                        HX(end+1,1)=xx(expt.tix);
                        HY(end+1,1)=yy(expt.tix);

                    otherwise
                        error(['Invalid name for x-axis: ' xfeat])
                end
                
                XX(end+1,:)=xx;
                YY(end+1,:)=yy;

            end
            
            if expt.N==1 %dummy trace to make coloring by segment work when only one trace
                XX=[XX,nan(size(XX))];
                YY=[YY,nan(size(YY))];
            end
            
            %line segment below markers:
            if expt.nS>1
                line(XX,YY,'color',[0.5,0.5,0.5],'linestyle','-','tag','oscar_line');
                line(HX,HY,'color','k','linestyle','-','linewidth',1.5,'tag','oscar_line');
            end
            
            %markers
            
            ax.ColorOrderIndex=1;
            hmAll=line(XX(:,expt.include)',YY(:,expt.include)','marker','o',...
                'linestyle','none','tag','oscar_line');
            
            hmAllex=line(XX(:,~expt.include)',YY(:,~expt.include)','marker','x',...
                'linestyle','none','tag','oscar_line');
            
            axis tight
            
            for i=1:expt.nS
                line(HX(i),HY(i),'marker','o','linestyle','none',...
                    'color',hmAll(i).Color,'markerfacecolor',hmAll(i).Color,'tag','oscar_line');
            end
            
%             for i=1:expt.nS
%                 ax.ColorOrderIndex=i;
%                 hm=line(XX(i),YY(i),'marker','o','linestyle','none','tag','oscar_line');
%                 line(HX(i),HY(i),'marker','o','linestyle','none',...
%                     'color',hm.Color,'markerfacecolor',hm.Color,'tag','oscar_line');
%                 hmAll(i)=hm;
%             end
%             axis tight
            
            if expt.nS>1
                legend(hmAll,[expt.segment(:).name],'AutoUpdate','off')
            end
            
            if xfeat=="segment"
                xticks(1:expt.nS)
%                 xticklabels({expt.segment.name})
                xlim([0.5,expt.nS+0.5])
            end
            
            if (xfeat=="periodMean"||xfeat=="Tpsd") && (yfeat=="periodMean"||yfeat=="Tpsd")
                line(xlim,xlim,'tag','oscar_line')
            end
            
            xlabel(xfeat)
            ylabel(yfeat)
            
        end
        
        function updatePlots(expt,type)
            
            if ~exist('type','var')
                type='all';
            end
            
            if ~isempty(expt.fig_handles)
                
            for i=1:length(expt.fig_handles)
                thisfig=expt.fig_handles(i);
                figData=thisfig.UserData;
                plotType=figData{1};
                
                if isequal(plotType,type)||type=="all"
                    figure(thisfig);
                    switch plotType
                        case 'trace'
                            whichPlot=figData{2};
                            showPts=figData{3};

                            nPlots=length(whichPlot);
                            ax=thisfig.Children;
                            for j=1:nPlots
                                axes(ax(nPlots+1-j));
                                expt.plot_t(whichPlot{j},showPts)
                            end

                        case 'psd'
                            expt.plot_psd()

                        case 'feature'
                            featPlotType=figData{4};
                            switch featPlotType
                                case 'per-period'
                                    expt.plot_features_periods();
                                case 'per-trace'
                                    expt.plot_features_trace();
                            end
                        case 'results'
                            
                    end
                end
            end
            
            figure(expt.active_fig);
            end
            
        end
        
        function commonKeypress(expt,src,event)
            expt.active_fig=src;
            switch(event.Key)
                case {'leftarrow'}
                    %previous trace
                    if expt.tix>1
                        expt.tix=expt.tix-1;
                        expt.updatePlots()
                    end
                    
                case {'rightarrow'}
                    %next trace
                    if expt.tix<expt.N
                        expt.tix=expt.tix+1;
                        expt.updatePlots()
                    end
                    
                case 'e'
                    %set enpoints of intervals
                    
                case 'f'
                    %select features to plot
                    expt.selectFeaturesPopup()
                    
                case 'g'
                    %pre-toggle actions
                    if expt.groupMode
                        expt.includeG=expt.include; %store include list
%                         expt.ixG=expt.tix; %store last group ix
                        expt.tix=find(expt.groupT==expt.tix,1,'first'); %determine new tix: first trace from currently in-focus group
                    else
                        expt.includeT=expt.include;
%                         expt.ixT=expt.tix;
                        expt.tix=expt.groupT(expt.tix); %determine new tix: group of currently in-focus trace 
                    end
                    
                    %toggle grouped mode
                    expt.groupMode=~expt.groupMode;
                    
                    %post-toggle actions
                    if expt.groupMode
                        expt.group=expt.groupID;
                        expt.include=expt.includeG;
                        expt.N=expt.nG;
                    else
                        expt.group=expt.groupT;
                        expt.include=expt.includeT;
                        expt.N=expt.nX;
                    end
                    
                    expt.preprocess()
                    expt.compute_features()
                    expt.updatePlots()
                    expt.displayResults(0)
                    
                case 'i'
                    %toggle include
                    expt.include(expt.tix)=~expt.include(expt.tix);
                    
                    expt.updatePlots('feature')
                    expt.buildResultsTable()
                    expt.displayResults(0)
                    
                case 'p'
                    %parameter dialog
                    expt.parameterDialog()
                    %needs full pipeline helper function
                    
                case 'q'
                    %quit - close all windows
                    expt.quit();
                    
                case 'r'
                    %re-run analysis
                    expt.preprocess()
                    expt.compute_features()
                    expt.updatePlots()
                    
                case 's'
                    %save results
%                     expt.writeToExcel()
                    expt.saveResults()
                    
                    
            end
        end
        
        function figureCloseFcn(expt,~,~)
            thisfig=gcf;
%             expt.gobj_array(ismember(expt.gobj_array, thisfig.Children))=[];
            expt.fig_handles(ismember(expt.fig_handles,thisfig))=[];
            delete(thisfig)
        end
        
        function resfigCloseFcn(expt,~,~)
            thisfig=gcf;
            if isequal(expt.resfig,thisfig) 
                expt.resfig=[];
                delete(thisfig)
            end
        end
        
        function selectFeaturesPopup(expt)
            
            expt.featureSelectDlg=figure('Name','Select features to plot','NumberTitle','off','MenuBar','none');
            expt.featureSelectDlg.KeyPressFcn=@expt.commonKeypress;
            expt.featureSelectDlg.Position(3:4)=[300,150];
            uibg=uibuttongroup(expt.featureSelectDlg,'Position',[0.1,0.75,0.8,0.2],'SelectionChangedFcn',@changeRadio);
%             uibg=uibuttongroup(figh,'Position',[0,80,300,40],'SelectionChangedFcn',@changeRadio);
            uirb1=uicontrol(uibg,'Style','radiobutton','String','per-period');
            uirb1.Units='normalized'; uirb1.Position=[0,0,0.5,1];
            uirb2=uicontrol(uibg,'Style','radiobutton','String','per-trace');
            uirb2.Units='normalized'; uirb2.Position=[0.5,0,0.5,1];
            if ~strcmp(expt.featurePlotType,uirb1.String)
                uirb2.Value=true;
            end
            
            switch expt.featurePlotType
                case 'per-period'
                    selections={[{'t'},{'segment'},expt.fnames_periods],expt.fnames_periods};
                case 'per-trace'
                    selections={[{'t'},{'segment'},expt.fnames_trace],expt.fnames_trace};
            end
            
            uit=uitable(expt.featureSelectDlg,'Units','Normalized','Position',[0.1,0.3,0.8,0.4]);
            uit.RowName=[];
            uit.ColumnWidth={118,118};
            uit.ColumnEditable=true;
            uit.ColumnName={'x feature','y feature'};
            uit.ColumnFormat=selections;
            uit.Data={expt.xfeature,expt.yfeature};
            uit.CellEditCallback=@cellEdit;
            
            if expt.featurePlotType=="per-trace"
                uibg.SelectedObject=uirb2;
            end
            
            uicontrol('Units','Normalized','Position',[0.2,0.05,0.6,0.2],'String','Plot',...
              'Callback',@launchFeatureUpdate);
          
            function launchFeatureUpdate(~,~)
                
                for i=1:length(expt.fig_handles)
                    if strcmp(expt.fig_handles(i).UserData,'feature')
                        expt.active_fig=expt.fig_handles(i);
                    end
                end
                
                expt.updatePlots('feature')
            end
            
            function changeRadio(~,cbdata)
                switch cbdata.NewValue.String
                    case 'per-period'
                        valid_features={[{'t'},{'segment'},expt.fnames_periods],expt.fnames_periods};
                        uit.ColumnFormat=valid_features;
                        if ~ismember(uit.Data(1),valid_features{1})
                            expt.xfeature=valid_features{1}{1};
                            uit.Data(1)=valid_features{1}(1);
                        end
                        if ~ismember(uit.Data(2),valid_features{2})
                            expt.yfeature=valid_features{2}{1};
                            uit.Data(2)=valid_features{2}(1);
                        end
                            
                    case 'per-trace'
                        valid_features={[{'segment'},expt.fnames_trace],expt.fnames_trace};
                        uit.ColumnFormat=valid_features;
                        if ~ismember(uit.Data(1),valid_features{1})
                            expt.xfeature=valid_features{1}{1};
                            uit.Data(1)=valid_features{1}(1);
                        end
                        if ~ismember(uit.Data(2),valid_features{2})
                            expt.yfeature=valid_features{2}{1};
                            uit.Data(2)=valid_features{2}(1);
                        end
                end
                expt.featurePlotType=cbdata.NewValue.String;
            end
            
            function cellEdit(~,cbdata)
                if cbdata.Indices(2)==1
                    expt.xfeature=cbdata.EditData;
                elseif cbdata.Indices(2)==2
                    expt.yfeature=cbdata.EditData;
                end
            end

        end
        
       
        function parameterDialog(expt)
            prompt = {'Detrend:','Filter:','Threshold:'};
            title = 'Enter parameter values';
            dims = [1,50];
            oldtrend=expt.trendParam;
            oldfilt=expt.filterParam;
            oldfeat=expt.featureParam;
            definput = {num2str(oldtrend), num2str(oldfilt), num2str(oldfeat,'%g, %g')};
            opts.WindowStyle='normal';
            answer = inputdlg(prompt,title,dims,definput,opts);
            if isempty(answer)
                return
            else
                doUpdate=false;
                newtrend = str2double(answer{1});
                if newtrend~=oldtrend
                    expt.trendParam=newtrend;
                    doUpdate=true;
                end
                
                newfilt = str2double(answer{2});
                if newfilt~=oldfilt
                    expt.filterParam=newfilt;
                    doUpdate=true;
                end
                
                newfeat = str2num(answer{3});
                if any(newfeat~=oldfeat)
                    expt.featureParam=newfeat;
                    doUpdate=true;
                end
                
                if doUpdate==true
                    expt.preprocess();
                    expt.compute_features();
                    expt.updatePlots()
                    expt.displayResults(0);
                end
            end
        end
        
        function quit(expt)
            if ~isempty(expt.fig_handles)
                close(expt.fig_handles)
            end
            if ~isempty(expt.resfig)
                close(expt.resfig)
            end
            if ~isempty(expt.featureSelectDlg)
                close(expt.featureSelectDlg)
            end
        end
    end
end

        
function offset = getFlattenOffset(Xnorm, Xt, flattenMethod)

    X=Xnorm;
    doNeg=false;
    if endsWith(flattenMethod,"Trend")
        flattenMethod=strrep(flattenMethod,'Trend','');
        X=Xt;
    elseif endsWith(flattenMethod,"Detrend")
        flattenMethod=strrep(flattenMethod,'Detrend','');
        X=Xnorm-Xt; 
        doNeg=true;%reverese sign makes this work...
    end
        
    offset=0;
    switch flattenMethod
        case 'none'
            offset=0;
        case 'first'
            offset=X(1,:);
        case 'min'
            offset=min(X,[],1);
        case 'max'
            offset=max(X,[],1);
        case 'mean'
            offset=mean(X,1);
        case 'median'
            offset=median(X,1);
        otherwise
            error(['Unknown flatten method: ', flattenMethod])
    end
    
    if doNeg
        offset=-offset;
    end
    
end