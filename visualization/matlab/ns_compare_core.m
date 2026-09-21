function app=ns_compare_core(varargin)
%NS_COMPARE_CORE Fixed-scale parent/child comparison of real frozen fields.
% Parent is ST061-P. Both children are geometry experiments, NOT accepted PDEs.
% Time is evaluated continuously. Neither axes nor field coordinates are stretched.
% This new UI has NOT been executed in native MATLAB in its creation environment.
p=inputParser;addParameter(p,'Visible','on');parse(p,varargin{:});
D=load(fullfile(fileparts(mfilename('fullpath')),'data','st063_models.mat'));
if ~iscell(D.models),D.models=num2cell(D.models);end
models=D.models;ids=cellfun(@(m)m.id,models,'UniformOutput',false);
parent=models{find(strcmp(ids,'ST061-P'),1)};child=models{find(strcmp(ids,'ST063-G2R'),1)};
state=struct('time',.5,'iso',.25,'radius',.65,'height',.90);busy=false;closed=false;last=tic;
fig=uifigure('Name','ST063 | real axial-core comparison','Position',[30 30 1470 1040],'Visible',p.Results.Visible);
g=uigridlayout(fig,[5 2]);g.RowHeight={40,'1x',200,100,34};g.ColumnWidth={'1x','1x'};
h=uilabel(g,'Text','FROZEN FIELD COMPARISON | identical seeds, physical aspect, camera, isovalue and color scale','FontSize',16,'FontWeight','bold');h.Layout.Row=1;h.Layout.Column=[1 2];
ax1=uiaxes(g);ax1.Layout.Row=2;ax1.Layout.Column=1;ax2=uiaxes(g);ax2.Layout.Row=2;ax2.Layout.Column=2;
profile=uiaxes(g);profile.Layout.Row=3;profile.Layout.Column=[1 2];
c=uigridlayout(g,[3 6]);c.Layout.Row=4;c.Layout.Column=[1 2];c.RowHeight={26,26,26};c.ColumnWidth={100,190,190,'1x',180,120};
h=uilabel(c,'Text','New candidate');at(h,1,1);cd=uidropdown(c,'Items',{'ST063-G1R','ST063-G2R'},'Value','ST063-G2R','ValueChangedFcn',@changeChild);at(cd,1,2);
quantity=uidropdown(c,'Items',{'Axial vorticity','Angular velocity','Full vorticity','Speed','Full NS residual'},'Value','Axial vorticity','ValueChangedFcn',@changeQuantity);at(quantity,1,3);
mode=uidropdown(c,'Items',{'Streamlines + isosurface','Isosurface only','XZ scalar slice'},'Value','Streamlines + isosurface','ValueChangedFcn',@refresh);at(mode,1,4);
viewChoice=uidropdown(c,'Items',{'Central region','Whole support'},'Value','Central region','ValueChangedFcn',@changeView);at(viewChoice,1,5);
play=uibutton(c,'Text','Play','ButtonPushedFcn',@toggle);at(play,1,6);
ih=uilabel(c,'Text','Level = 0.25');at(ih,2,1);is=uislider(c,'Limits',[.01,1.5],'Value',state.iso,'MajorTicks',[],'ValueChangedFcn',@changeIso);is.Layout.Row=2;is.Layout.Column=[2 4];
reset=uibutton(c,'Text','Reset common camera','ButtonPushedFcn',@resetCamera);at(reset,2,5);
export=uibutton(c,'Text','Export PNG','ButtonPushedFcn',@exportPng);at(export,2,6);
h=uilabel(c,'Text','Physical time');at(h,3,1);ts=uislider(c,'Limits',[.25,.75],'Value',.5,'MajorTicks',[],...
 'ValueChangingFcn',@(s,e)setTime(e.Value,true),'ValueChangedFcn',@(s,e)setTime(e.Value,false));ts.Layout.Row=3;ts.Layout.Column=[2 5];
th=uilabel(c,'Text','t = 0.50000');at(th,3,6);
status=uilabel(g,'Text','Initializing fixed-field comparison','FontSize',11);status.Layout.Row=5;status.Layout.Column=[1 2];
for a=[ax1 ax2],view(a,-35,22);grid(a,'on');colormap(a,parula(256));colorbar(a);end
link=[];tm=timer('ExecutionMode','fixedSpacing','Period',.25,'BusyMode','drop','TimerFcn',@tick);fig.CloseRequestFcn=@closeApp;
app=struct('Figure',fig,'SetTime',@(t)setTime(t,false),'SetIso',@setIso,'SetCandidate',@setCandidate,'Close',@closeApp,'TimeSlider',ts,'ParentAxes',ax1,'ChildAxes',ax2,'ProfileAxes',profile);render(false);
 function at(h,r,k),h.Layout.Row=r;h.Layout.Column=k;end
 function setCandidate(s),cd.Value=s;changeChild();end
 function changeChild(varargin),child=models{find(strcmp(ids,cd.Value),1)};render(false);end
 function changeView(varargin)
  if strcmp(viewChoice.Value,'Central region'),state.radius=.65;state.height=.90;else,state.radius=2;state.height=2;end
  render(false);resetCamera();
 end
 function changeQuantity(varargin)
  switch quantity.Value
   case 'Angular velocity',state.iso=.10;
   case 'Speed',state.iso=.15;
   case 'Full NS residual',state.iso=.015;
   case 'Full vorticity',state.iso=.6;
   otherwise,state.iso=.25;
  end
  is.Value=state.iso;render(false);
 end
 function setIso(v),validateattributes(v,{'numeric'},{'scalar','>=',.01,'<=',1.5});state.iso=v;is.Value=v;render(false);end
 function changeIso(src,evt),setIso(src.Value);end
 function refresh(varargin),render(false);end
 function setTime(v,preview)
  validateattributes(v,{'numeric'},{'scalar','>=',.25,'<=',.75});state.time=v;ts.Value=v;th.Text=sprintf('t = %.5f',v);
  if preview&&toc(last)<.15,return,end;render(preview);
 end
 function toggle(varargin)
  if strcmp(tm.Running,'on'),stop(tm);play.Text='Play';render(false);else,if state.time>=.75,state.time=.25;end;play.Text='Pause';start(tm);end
 end
 function tick(varargin)
  setTime(min(.75,state.time+.005),true);if state.time>=.75,stop(tm);play.Text='Play';render(false);end
 end
 function resetCamera(varargin)
  for a=[ax1 ax2],view(a,-35,22);end
 end
 function exportPng(varargin)
  stop(tm);play.Text='Play';render(false);[n,pth]=uiputfile('*.png','Export actual comparison','ST063_comparison.png');if ~isequal(n,0),exportapp(fig,fullfile(pth,n));end
 end
 function render(preview)
  if busy||closed,return,end;busy=true;guard=onCleanup(@unlock);tic;
  try
   % Preserve camera while both views change the PHYSICAL field with time.
   [aa,ee]=view(ax1);if ~isempty(link),delete(link);link=[];end
   if preview,n=27;nz=39;count=48;else,n=49;nz=73;count=200;end
   show(ax1,parent,n,nz,count);show(ax2,child,n,nz,count);
   if strcmp(mode.Value,'XZ scalar slice')
    view(ax1,2);view(ax2,2);
   else
    view(ax1,aa,ee);view(ax2,aa,ee);link=linkprop([ax1 ax2],{'CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle'});
   end
   z=linspace(-state.height,state.height,401);s=zeros(size(z));a=ns_cylindrical(parent,s,z,state.time,false);b=ns_cylindrical(child,s,z,state.time,false);
   % The profile is always real on-axis angular velocity B = u_theta/r limit.
   cla(profile);plot(profile,z,a.B,'LineWidth',1.6);hold(profile,'on');plot(profile,z,b.B,'--','LineWidth',1.6);hold(profile,'off');grid(profile,'on');
   xlim(profile,[-state.height state.height]);ylim(profile,[-.1 .65]);xlabel(profile,'Physical z');ylabel(profile,'On-axis angular velocity B');
   legend(profile,{parent.id,child.id},'Interpreter','none','Location','best');title(profile,'Axial rotation profile: a flatter, longer plateau is not merely a stretched camera');
   ih.Text=sprintf('Level = %.4g',state.iso);th.Text=sprintf('t = %.5f',state.time);
   status.Text=sprintf('%dx%dx%d local samples | %.2f s | Same absolute scales. Midplane disk is retained. Both fields FAIL the original full NS target.',n,n,nz,toc);
   setappdata(fig,'LastRenderError','');drawnow limitrate nocallbacks;last=tic;
  catch e,setappdata(fig,'LastRenderError',e.message);status.Text=['Error: ' e.message];warning('ns:Comparison','%s',e.message);end
  clear guard
 end
 function show(ax,m,n,nz,count)
  cla(ax);hold(ax,'on');rmax=state.radius;hmax=state.height;
  if strcmp(mode.Value,'XZ scalar slice')
   x=linspace(-rmax,rmax,151);z=linspace(-hmax,hmax,181);[X,Z]=meshgrid(x,z);pts=[X(:),zeros(numel(X),1),Z(:)];
   [u,p,w,res]=ns_evaluate(m,pts,state.time); %#ok<ASGLU>
   switch quantity.Value
    case 'Angular velocity',dc=ns_cylindrical(m,X(:).^2,Z(:),state.time,false);v=dc.B;lim=[-.1 .65];
    case 'Axial vorticity',v=w(:,3);lim=[-.2 1.5];
    case 'Full vorticity',v=vecnorm(w,2,2);lim=[0 3];
    case 'Speed',v=vecnorm(u,2,2);lim=[0 1];
    otherwise,v=vecnorm(res,2,2);lim=[0 .06];
   end
   imagesc(ax,x,z,reshape(v,size(X)));ax.YDir='normal';axis(ax,'equal');xlim(ax,[-rmax rmax]);ylim(ax,[-hmax hmax]);view(ax,2);caxis(ax,lim);xlabel(ax,'x');ylabel(ax,'z');title(ax,[m.id ' | physical XZ slice'],'Interpreter','none');hold(ax,'off');return
  end
  x=linspace(-rmax,rmax,n);z=linspace(-hmax,hmax,nz);[xx,yy]=meshgrid(x,x);s=xx(:).^2+yy(:).^2;
  d=ns_cylindrical(m,s,z,state.time,true);[X,Y,Z]=meshgrid(x,x,z);sz=size(X);
  U=reshape(xx(:).*d.A-yy(:).*d.B,sz);V=reshape(yy(:).*d.A+xx(:).*d.B,sz);W=reshape(d.C,sz);
  oz=2*(d.B+s.*d.Bs);speed=sqrt(U.^2+V.^2+W.^2);
  switch quantity.Value
   case 'Angular velocity',value=reshape(d.B,sz);lim=[-.1 .65];
   case 'Axial vorticity',value=reshape(oz,sz);lim=[-.2 1.5];
   case 'Full vorticity',value=reshape(sqrt(s.*d.Bz.^2+s.*(d.Az-2*d.Cs).^2+oz.^2),sz);lim=[0 3];
   case 'Speed',value=speed;lim=[0 1];
   otherwise,value=reshape(sqrt(s.*(d.Ra.^2+d.Rb.^2)+d.Rc.^2),sz);lim=[0 .06];
  end
  if state.iso>min(value(:))&&state.iso<max(value(:))
   sf=isosurface(X,Y,Z,value,state.iso);
   if ~isempty(sf.faces),q=patch(ax,sf,'FaceColor','flat','FaceVertexCData',state.iso*ones(size(sf.vertices,1),1),'EdgeColor','none','FaceAlpha',.22);isonormals(X,Y,Z,value,q);end
  end
  if strcmp(mode.Value,'Streamlines + isosurface')
   [rr,zz,ang]=ndgrid([.07 .14 .28 .45 .58],[-.55 -.275 0 .275 .55],(0:7)*2*pi/8);
   ix=round(linspace(1,numel(rr),count));sx=rr(ix).*cos(ang(ix));sy=rr(ix).*sin(ang(ix));zs=zz(ix);
   F=stream3(X,Y,Z,U,V,W,sx(:),sy(:),zs(:),[.18 1400]);B=stream3(X,Y,Z,-U,-V,-W,sx(:),sy(:),zs(:),[.18 1400]);
   I=griddedInterpolant({x,x,z},value,'linear','none');pts=[];cs=[];
   for k=1:numel(F)
    q=[flipud(trim(B{k},1.5));trim(F{k},1.5)];if size(q,1)<2,continue,end
    pts=[pts;q;nan(1,3)];cs=[cs;I(q(:,2),q(:,1),q(:,3));nan]; %#ok<AGROW>
   end
   if ~isempty(pts),surface(ax,[pts(:,1)';pts(:,1)'],[pts(:,2)';pts(:,2)'],[pts(:,3)';pts(:,3)'],[cs';cs'],'FaceColor','none','EdgeColor','interp','LineWidth',.55);end
  end
  axis(ax,'equal');xlim(ax,[-rmax rmax]);ylim(ax,[-rmax rmax]);zlim(ax,[-hmax hmax]);caxis(ax,lim);grid(ax,'on');xlabel(ax,'x');ylabel(ax,'y');zlabel(ax,'z');title(ax,[m.id ' | instantaneous field'],'Interpreter','none');hold(ax,'off');
 end
 function q=trim(q,cap)
  if size(q,1)<2,return,end;a=[0;cumsum(vecnorm(diff(q),2,2))];k=find(a<=cap,1,'last');q=q(1:k,:);
 end
 function unlock(),busy=false;end
 function closeApp(varargin),if closed,return,end;closed=true;if isvalid(tm),stop(tm);delete(tm);end;if ~isempty(link),delete(link);end;if isvalid(fig),delete(fig);end;end
end
