function app=ns_explorer(dataFile,varargin)
%NS_EXPLORER Interactive ST054 velocity/pressure/vorticity/residual explorer.
%   ns_explorer                  Open bundled Q2/M3 with continuous time.
%   ns_explorer('grid.mat')       Import PR665 [time,x,y,z] velocity arrays.
%   a=ns_explorer(...,'Visible','off') is useful for a GUI smoke test.
% Display controls NEVER modify frozen candidate coefficients or viscosity.
% MATLAB R2021a+; base MATLAB only. Python is NOT required to use the app.
if nargin<1||isempty(dataFile),dataFile=fullfile(fileparts(mfilename('fullpath')),'data','st054_models.mat');end
p=inputParser;addParameter(p,'Visible','on');parse(p,varargin{:});
if ~isfile(dataFile),error('ns:Data','Missing MAT file: %s. Extract the complete package or run export_st054.py.',dataFile);end
D=load(dataFile);
if ~isfield(D,'schema')||~strcmp(D.schema,'ns_matlab_spectral_v1')
    D=ns_import_grid(dataFile);
end
if ~iscell(D.models),D.models=num2cell(D.models);end
models=D.models;idx=1;m=models{idx};
if ~isfield(m,'mode')||~ismember(m.mode,{'spectral','sampled'}),error('ns:Data','Unsupported model schema.');end
if isfield(m,'pde_validated')&&m.pde_validated,error('ns:Claims','This viewer does not accept unverified promoted claims.');end
state=struct('t',(m.tmin+m.tmax)/2,'slice',0,'radius',1,'seedz',0,'count',200,...
    'length',5,'iso',1.5,'alpha',.20,'extent',2,'cut',2,'colorScale',1);
V=[];cacheTime=NaN;cacheN=0;busy=false;lastDraw=tic;closed=false;controlLabels=struct();
fig=uifigure('Name','NS field explorer | frozen research candidates','Position',[40 40 1440 900],...
    'Visible',p.Results.Visible,'Color',[.96 .97 .98]);
fig.CloseRequestFcn=@closeApp;
g=uigridlayout(fig,[4 3]);g.RowHeight={50,'1x',86,40};g.ColumnWidth={300,'1.2x','1x'};
g.Padding=[12 10 12 10];g.RowSpacing=8;
head=uilabel(g,'Text','NS FIELD EXPLORER  |  frozen field, not an accepted NS solution',...
    'FontSize',18,'FontWeight','bold');head.Layout.Row=1;head.Layout.Column=[1 3];
panel=uipanel(g,'Title','Field / display controls','Scrollable','on');panel.Layout.Row=2;panel.Layout.Column=1;
cg=uigridlayout(panel,[18 2]);cg.ColumnWidth={124,'1x'};cg.RowHeight=repmat({32},1,18);cg.Padding=[8 8 8 8];cg.RowSpacing=7;
items=cellfun(@(q)q.id,models,'UniformOutput',false);
label('Candidate',1);candidate=uidropdown(cg,'Items',items,'Value',m.id,'ValueChangedFcn',@changeCandidate);place(candidate,1,2);
label('3-D view',2);mode=uidropdown(cg,'Items',{'Streamlines','Vorticity surface','Velocity arrows','Streamlines + surface'},'Value','Streamlines','ValueChangedFcn',@refresh);place(mode,2,2);
label('Slice scalar',3);scalar=uidropdown(cg,'Items',scalarItems(),'Value','Speed','ValueChangedFcn',@refresh);place(scalar,3,2);
label('Line color',4);color=uidropdown(cg,'Items',{'Speed','Height','Vorticity','Residual'},'Value','Speed','ValueChangedFcn',@refresh);place(color,4,2);
label('Slice plane',5);plane=uidropdown(cg,'Items',{'XZ (fixed y)','XY (fixed z)','YZ (fixed x)'},'Value','XZ (fixed y)','ValueChangedFcn',@refresh);place(plane,5,2);
sl.slice=control('Slice offset',6,[-1.95 1.95],'slice');
sl.radius=control('Seed radius',7,[.15 1.6],'radius');
sl.seedz=control('Seed z center',8,[-1.3 1.3],'seedz');
sl.count=control('Line count',9,[40 320],'count');
sl.length=control('Arc length cap',10,[.5 12],'length');
sl.iso=control('|omega| level',11,[.02 6],'iso');
sl.alpha=control('Surface opacity',12,[.03 .65],'alpha');
sl.extent=control('Radial view',13,[.3 2],'extent');
sl.cut=control('Cutaway y <=',14,[-1.9 2],'cut');
sl.colorScale=control('Color limit x',15,[.25 2],'colorScale');
quality=uidropdown(cg,'Items',{'Fast (33^3)','Balanced (49^3)','Fine (65^3)'},'Value','Balanced (49^3)','ValueChangedFcn',@refresh);
label('Spatial quality',16);place(quality,16,2);
reset=uibutton(cg,'Text','Reset view','ButtonPushedFcn',@resetView);place(reset,17,1);
imp=uibutton(cg,'Text','Import grid MAT','ButtonPushedFcn',@importGrid);place(imp,17,2);
png=uibutton(cg,'Text','Export PNG','ButtonPushedFcn',@exportPng);place(png,18,1);
mat=uibutton(cg,'Text','Export view MAT','ButtonPushedFcn',@exportMat);place(mat,18,2);
ax3=uiaxes(g);ax3.Layout.Row=2;ax3.Layout.Column=2;
ax2=uiaxes(g);ax2.Layout.Row=2;ax2.Layout.Column=3;
view(ax3,-37,22);grid(ax3,'on');axis(ax3,'equal');colormap(ax3,parula(256));colormap(ax2,parula(256));
bar3=colorbar(ax3);bar2=colorbar(ax2);
foot=uigridlayout(g,[2 4]);foot.Layout.Row=3;foot.Layout.Column=[1 3];
foot.RowHeight={28,38};foot.ColumnWidth={100,'1x',160,140};foot.Padding=[4 0 4 0];
play=uibutton(foot,'Text','Play','ButtonPushedFcn',@togglePlay);play.Layout.Row=1;play.Layout.Column=1;
timeText=uilabel(foot,'Text','','FontWeight','bold');timeText.Layout.Row=1;timeText.Layout.Column=3;
clockMode=uilabel(foot,'Text','Physical time');clockMode.Layout.Row=1;clockMode.Layout.Column=4;
help=uilabel(foot,'Text','Drag: responsive low-detail preview. Release: rebuild full detail. Mouse: rotate / zoom 3-D.');help.Layout.Row=1;help.Layout.Column=2;
time=uislider(foot,'Limits',[m.tmin m.tmax],'Value',state.t,'MajorTicks',linspace(m.tmin,m.tmax,5),...
    'ValueChangingFcn',@(src,e)timeChanged(e.Value,true),'ValueChangedFcn',@(src,e)timeChanged(e.Value,false));
time.Layout.Row=2;time.Layout.Column=[1 4];time.Tag='TimeSlider';
status=uilabel(g,'Text','Preparing real coefficient evaluation...','FontSize',11);
status.Layout.Row=4;status.Layout.Column=[1 3];
tm=timer('ExecutionMode','fixedSpacing','Period',.18,'BusyMode','drop','TimerFcn',@tick);
app=struct('Figure',fig,'TimeSlider',time,'PlayButton',play,'ModelDropdown',candidate,'Axes3D',ax3,'AxesSlice',ax2,'Status',status,...
    'SetTime',@(t)setTime(t),'SetControl',@setControl,'SetMode',@setMode,...
    'SetScalar',@setScalar,'SetCandidate',@setCandidate,'Close',@()closeApp(),...
    'Snapshot',@()snapshot(),'Refresh',@()render(false));
updateImportRestrictions();render(false);

    function label(text,row)
        h=uilabel(cg,'Text',text,'FontSize',11);place(h,row,1);
    end
    function place(h,row,col),h.Layout.Row=row;h.Layout.Column=col;end
    function h=control(text,row,limits,name)
        L=uilabel(cg,'Text',sprintf('%s: %.3g',text,state.(name)),'FontSize',11);place(L,row,1);controlLabels.(name)={L,text};
        h=uislider(cg,'Limits',limits,'Value',state.(name),'MajorTicks',[],...
            'ValueChangingFcn',@(src,e)controlChanged(name,e.Value,true,L,text),...
            'ValueChangedFcn',@(src,e)controlChanged(name,e.Value,false,L,text));place(h,row,2);h.Tag=name;
    end
    function controlChanged(name,val,preview,L,text)
        if strcmp(name,'count'),val=round(val);end
        state.(name)=val;L.Text=sprintf('%s: %.3g',text,val);
        if preview&&toc(lastDraw)<.12,return,end
        render(preview);
    end
    function timeChanged(val,preview)
        if closed,return,end
        state.t=max(m.tmin,min(m.tmax,val));timeText.Text=sprintf('t = %.5f',state.t);
        if preview&&toc(lastDraw)<.12,return,end
        render(preview);
    end
    function setTime(val)
        stop(tm);play.Text='Play';time.Value=val;timeChanged(val,false);
    end
    function setControl(name,val)
        if ~isfield(sl,name),error('ns:Control','Unknown display control.');end
        if val<sl.(name).Limits(1)||val>sl.(name).Limits(2),error('ns:Control','Control outside limits.');end
        sl.(name).Value=val;state.(name)=val;entry=controlLabels.(name);entry{1}.Text=sprintf('%s: %.3g',entry{2},val);render(false);
    end
    function setMode(value),mode.Value=value;render(false);end
    function setScalar(value),scalar.Value=value;render(false);end
    function setCandidate(value),candidate.Value=value;changeCandidate();end
    function s=snapshot()
        s=state;s.candidate=m.id;s.mode=mode.Value;s.scalar=scalar.Value;s.color=color.Value;
        s.plane=plane.Value;s.pde_validated=false;s.data_file=dataFile;
    end
    function refresh(varargin),render(false);end
    function changeCandidate(varargin)
        stop(tm);play.Text='Play';idx=find(strcmp(items,candidate.Value),1);m=models{idx};
        state.t=max(m.tmin,min(m.tmax,state.t));time.Limits=[m.tmin m.tmax];time.Value=state.t;
        cacheTime=NaN;V=[];updateImportRestrictions();render(false);
    end
    function updateImportRestrictions()
        scalar.Items=scalarItems();scalar.Value=scalar.Items{1};
        if strcmp(m.mode,'sampled'),color.Items={'Speed','Height','Vorticity'};else,color.Items={'Speed','Height','Vorticity','Residual'};end
        color.Value='Speed';
    end
    function out=scalarItems()
        out={'Speed','Swirl velocity','Axial velocity','Vorticity'};
        if strcmp(m.mode,'spectral'),out=[out,{'Pressure','Axial pressure force','Full NS residual'}];end
    end
    function togglePlay(varargin)
        if strcmp(tm.Running,'on'),stop(tm);play.Text='Play';render(false);
        else
            if state.t>=m.tmax,state.t=m.tmin;end
            play.Text='Pause';start(tm);
        end
    end
    function tick(varargin)
        if closed||~isvalid(fig),return,end
        state.t=min(m.tmax,state.t+(m.tmax-m.tmin)/100);time.Value=state.t;render(true);
        if state.t>=m.tmax,stop(tm);play.Text='Play';render(false);end
    end
    function resetView(varargin)
        view(ax3,-37,22);state.extent=2;sl.extent.Value=2;state.cut=2;sl.cut.Value=2;render(false);
    end
    function importGrid(varargin)
        [f,pth]=uigetfile('*.mat','Open PR665 [time,x,y,z] grid export');
        if isequal(f,0),return,end
        try,ns_explorer(fullfile(pth,f));catch err,uialert(fig,err.message,'Import error');end
    end
    function exportPng(varargin)
        stop(tm);play.Text='Play';render(false);
        [f,pth]=uiputfile('*.png','Save explorer view',[m.id '_view.png']);if isequal(f,0),return,end
        exportapp(fig,fullfile(pth,f));
    end
    function exportMat(varargin)
        stop(tm);play.Text='Play';render(false);
        [f,pth]=uiputfile('*.mat','Save current sampled volume',[m.id '_view.mat']);if isequal(f,0),return,end
        view_state=snapshot();sampled_volume=V;candidate_metadata=m; %#ok<NASGU>
        save(fullfile(pth,f),'view_state','sampled_volume','candidate_metadata','-v7');
    end
    function render(preview)
        if busy||closed||~isvalid(fig),return,end
        busy=true;guard=onCleanup(@unlock);started=tic;
        try
            if preview,n=25;else
                switch quality.Value,case 'Fast (33^3)',n=33;case 'Fine (65^3)',n=65;otherwise,n=49;end
            end
            if isempty(V)||cacheTime~=state.t||cacheN~=n
                V=ns_volume(m,state.t,n);cacheTime=state.t;cacheN=n;
            end
            timeText.Text=sprintf('t = %.5f',state.t);
            setappdata(fig,'LastRenderError','');
            opts=struct('mode',mode.Value,'color',color.Value,'scalar',scalar.Value,'plane',plane.Value,'preview',preview);
            ns_render_scene(ax3,ax2,bar3,bar2,V,m,state,opts);
            if isempty(V.residual),detail='Imported velocity grid: linear interpolation; NO pressure/NS residual.';
            else,detail=sprintf('Rendered-grid max |R| = %.5g (NOT an independent acceptance test).',max(V.residual(:)));end
            if preview,labelQuality='DRAG/PLAY PREVIEW';else,labelQuality='FULL VIEW';end
            status.Text=sprintf('%s | %d^3 nodes | %.2f s | %s | PDE target NOT met.',labelQuality,n,toc(started),detail);
            drawnow limitrate nocallbacks;lastDraw=tic;
        catch err
            setappdata(fig,'LastRenderError',err.message);status.Text=['Render error: ' err.message];warning('ns:Render','%s',getReport(err,'basic'));
        end
        clear guard
    end
    function unlock(),busy=false;end
    function closeApp(varargin)
        if closed,return,end
        closed=true;
        if isvalid(tm),stop(tm);delete(tm);end
        if isvalid(fig),delete(fig);end
    end
end
