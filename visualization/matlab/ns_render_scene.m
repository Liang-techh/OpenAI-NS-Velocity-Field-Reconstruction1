function ns_render_scene(ax3,ax2,bar3,bar2,V,m,state,options)
%NS_RENDER_SCENE View-only rendering; frozen coefficients are never modified.
mode=struct('Value',options.mode);color=struct('Value',options.color);
scalar=struct('Value',options.scalar);plane=struct('Value',options.plane);
[az,el]=view(ax3);cla(ax3);hold(ax3,'on');
if contains(mode.Value,'surface'),drawIso();end
if contains(mode.Value,'Streamlines'),drawLines(options.preview);end
if strcmp(mode.Value,'Velocity arrows'),drawVectors();end
xlabel(ax3,'x');ylabel(ax3,'y');zlabel(ax3,'z');axis(ax3,'equal');grid(ax3,'on');
xlim(ax3,[-state.extent state.extent]);ylim(ax3,[-state.extent state.extent]);zlim(ax3,[V.z(1) V.z(end)]);view(ax3,az,el);
title(ax3,sprintf('%s | t=%.4f\n%s',m.id,state.t,mode.Value),'Interpreter','none');
drawSlice(options.preview);hold(ax3,'off');
    function drawIso()
        level=state.iso;
        if level>=max(V.vorticity(:))||level<=min(V.vorticity(:)),return,end
        S=isosurface(V.X,V.Y,V.Z,V.vorticity,level);
        if isempty(S.faces),return,end
        keep=all(reshape(S.vertices(S.faces(:),2),size(S.faces))<=state.cut,2);S.faces=S.faces(keep,:);
        h=patch(ax3,S,'FaceColor',[.3 .62 .72],'EdgeColor','none','FaceAlpha',state.alpha);
        if ~isempty(S.faces),isonormals(V.X,V.Y,V.Z,V.vorticity,h);end
        camlight(ax3,'headlight');lighting(ax3,'gouraud');
        bar3.Label.String='Isosurface |omega| (absolute level)';
    end
    function drawVectors()
        k=1:3:numel(V.x);U=V.U(k,k,k);W=V.W(k,k,k);VV=V.V(k,k,k);
        X=V.X(k,k,k);Y=V.Y(k,k,k);Z=V.Z(k,k,k);U(Y>state.cut)=NaN;
        % Fixed display gain, not per-frame normalization.
        quiver3(ax3,X,Y,Z,.2*U,.2*VV,.2*W,0);bar3.Label.String='Velocity arrows: fixed gain 0.2';
    end
    function drawLines(preview)
        count=round(state.count);L=state.length;
        if preview,count=min(count,48);L=min(L,2.5);end
        nang=ceil(count/25);a=(0:nang-1)*2*pi/nang;
        [rr,zz,aa]=ndgrid(state.radius*[.25 .45 .65 .85 1.05],state.seedz+[-.8 -.4 0 .4 .8],a);
        take=round(linspace(1,numel(rr),count));sx=rr(take).*cos(aa(take));sy=rr(take).*sin(aa(take));sz=zz(take);
        sx=sx(:);sy=sy(:);sz=sz(:);inside=abs(sx)<2&abs(sy)<2&abs(sz)<1.99;sx=sx(inside);sy=sy(inside);sz=sz(inside);
        maxv=min(2400,ceil(L/(.25*min(diff(V.x))))+12);
        F=stream3(V.X,V.Y,V.Z,V.U,V.V,V.W,sx,sy,sz,[.25,maxv]);
        B=stream3(V.X,V.Y,V.Z,-V.U,-V.V,-V.W,sx,sy,sz,[.25,maxv]);
        A=colorVolume();I=griddedInterpolant({V.y,V.x,V.z},A,'linear','none');
        points=[];colors=[];
        for k=1:numel(F)
            q=[flipud(trim(B{k},L/2));trim(F{k},L/2)];
            if size(q,1)<2,continue,end
            if strcmp(color.Value,'Height'),cc=q(:,3);else,cc=I(q(:,2),q(:,1),q(:,3));end
            q(q(:,2)>state.cut,:)=NaN;
            points=[points;q;NaN(1,3)];colors=[colors;cc;NaN]; %#ok<AGROW>
        end
        if ~isempty(points)
            % ONE surface object for all colored segments, not thousands of plot3 calls.
            surface(ax3,[points(:,1)';points(:,1)'],[points(:,2)';points(:,2)'],[points(:,3)';points(:,3)'],...
                [colors';colors'],'FaceColor','none','EdgeColor','interp','LineWidth',.6);
        end
        clim3=colorLimits(color.Value);caxis(ax3,clim3);bar3.Label.String=[color.Value ' (fixed scale)'];
    end
    function A=colorVolume()
        switch color.Value,case 'Speed',A=V.speed;case 'Vorticity',A=V.vorticity;
            case 'Residual',A=V.residual;otherwise,A=V.Z;end
    end
    function lim=colorLimits(name)
        switch name
            case 'Height',lim=[-2 2];
            case {'Speed'},lim=[0 .8];
            case {'Vorticity'},lim=[0 5];
            case 'Pressure',lim=[-.2 .2];
            case {'Swirl velocity','Axial velocity'},lim=[-.7 .7];
            case 'Axial pressure force',lim=[-.5 .5];
            otherwise,lim=[0 .06];
        end
        lim=lim*state.colorScale;
    end
    function drawSlice(preview)
        if preview,nn=61;else,nn=151;end
        coords=linspace(-2,2,nn);[a,b]=meshgrid(coords,coords);pos=state.slice;
        switch plane.Value
            case 'XZ (fixed y)',X=a;Y=pos+zeros(size(a));Z=b;horizontal='x';vertical='z';components=[1 3];
            case 'XY (fixed z)',X=a;Y=b;Z=pos+zeros(size(a));horizontal='x';vertical='y';components=[1 2];
            otherwise,X=pos+zeros(size(a));Y=a;Z=b;horizontal='y';vertical='z';components=[2 3];
        end
        if strcmp(m.mode,'spectral')
            [u,pv,om,rr,pf]=ns_evaluate(m,[X(:),Y(:),Z(:)],state.t);
            switch scalar.Value
                case 'Speed',C=sqrt(sum(u.^2,2));
                case 'Vorticity',C=sqrt(sum(om.^2,2));
                case 'Swirl velocity',C=(-Y(:).*u(:,1)+X(:).*u(:,2))./max(hypot(X(:),Y(:)),eps);
                case 'Axial velocity',C=u(:,3);
                case 'Pressure',C=pv;
                case 'Axial pressure force',C=pf(:,3);
                otherwise,C=sqrt(sum(rr.^2,2));
            end
        else
            u=[interp3(V.X,V.Y,V.Z,V.U,X(:),Y(:),Z(:)),interp3(V.X,V.Y,V.Z,V.V,X(:),Y(:),Z(:)),interp3(V.X,V.Y,V.Z,V.W,X(:),Y(:),Z(:))];
            switch scalar.Value
                case 'Vorticity',C=interp3(V.X,V.Y,V.Z,V.vorticity,X(:),Y(:),Z(:));
                case 'Swirl velocity',C=(-Y(:).*u(:,1)+X(:).*u(:,2))./max(hypot(X(:),Y(:)),eps);
                case 'Axial velocity',C=u(:,3);
                otherwise,C=sqrt(sum(u.^2,2));
            end
        end
        cla(ax2);imagesc(ax2,coords,coords,reshape(C,size(a)));ax2.YDir='normal';hold(ax2,'on');
        vx=reshape(u(:,components(1)),size(a));vy=reshape(u(:,components(2)),size(a));k=1:max(1,round(nn/18)):nn;
        quiver(ax2,a(k,k),b(k,k),.15*vx(k,k),.15*vy(k,k),0,'Color',[.92 .92 .92],'LineWidth',.65);
        hold(ax2,'off');axis(ax2,'equal');xlim(ax2,[-2 2]);ylim(ax2,[-2 2]);
        xlabel(ax2,horizontal);ylabel(ax2,vertical);caxis(ax2,colorLimits(scalar.Value));bar2.Label.String=[scalar.Value ' (fixed scale)'];
        title(ax2,sprintf('%s | offset=%.3f\nExact spectral slice; arrows use fixed gain',scalar.Value,pos));
        if strcmp(m.mode,'sampled'),title(ax2,[scalar.Value ' | interpolated grid slice']);end
    end
    function q=trim(q,cap)
        if isempty(q),return,end
        ds=sqrt(sum(diff(q).^2,2));last=find([0;cumsum(ds)]<=cap,1,'last');q=q(1:last,:);
    end
end
