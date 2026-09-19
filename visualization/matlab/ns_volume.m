function V=ns_volume(m,t,n)
%NS_VOLUME meshgrid-layout volumes. Analytic source evaluated at every node.
% The sampled-grid adapter is linear in time/space and does not supply p or R.
if strcmp(m.mode,'sampled')
    V=sampled(m,t,n);return
end
x=linspace(-2,2,n);z=x;[xx,yy]=meshgrid(x,x);s=xx(:).^2+yy(:).^2;
D=ns_cylindrical(m,s,z,t,true);xx=xx(:);yy=yy(:);sz=[n n n];
V.x=x;V.y=x;V.z=z;[V.X,V.Y,V.Z]=meshgrid(x,x,z);
V.U=reshape(xx.*D.A-yy.*D.B,sz);V.V=reshape(yy.*D.A+xx.*D.B,sz);V.W=reshape(D.C,sz);
V.P=reshape(D.P,sz);V.PforceZ=reshape(-D.Pz,sz);
V.speed=sqrt(V.U.^2+V.V.^2+V.W.^2);
V.swirl=reshape(sqrt(s).*D.B,sz);V.vorticity=reshape(sqrt(s.*D.Bz.^2+s.*(D.Az-2*D.Cs).^2+4*(D.B+s.*D.Bs).^2),sz);
V.residual=reshape(sqrt(s.*(D.Ra.^2+D.Rb.^2)+D.Rc.^2),sz);
V.divergence=reshape(D.divergence,sz);V.time=t;V.approximate=false;
end
function V=sampled(m,t,n)
validateattributes(t,{'numeric'},{'scalar','>=',m.tmin,'<=',m.tmax});
x=linspace(m.x(1),m.x(end),n);y=linspace(m.y(1),m.y(end),n);z=linspace(m.z(1),m.z(end),n);
[V.X,V.Y,V.Z]=meshgrid(x,y,z);V.x=x;V.y=y;V.z=z;
a=find(m.t<=t,1,'last');b=min(a+1,numel(m.t));w=0;
if a~=b,w=(t-m.t(a))/(m.t(b)-m.t(a));end
names={'u','v','w'};out={'U','V','W'};
for k=1:3
    C=m.(names{k});A=reshape(C(a,:,:,:),[numel(m.x),numel(m.y),numel(m.z)]);
    B=reshape(C(b,:,:,:),size(A));F=griddedInterpolant({m.y,m.x,m.z},permute((1-w)*A+w*B,[2 1 3]),'linear','none');
    V.(out{k})=F(V.Y,V.X,V.Z);
end
V.speed=sqrt(V.U.^2+V.V.^2+V.W.^2);r=hypot(V.X,V.Y);
V.swirl=(-V.Y.*V.U+V.X.*V.V)./max(r,eps);
[~,Uy,Uz]=gradient(V.U,x,y,z);[Vx,~,Vz]=gradient(V.V,x,y,z);[Wx,Wy,~]=gradient(V.W,x,y,z);
V.vorticity=sqrt((Wy-Vz).^2+(Uz-Wx).^2+(Vx-Uy).^2);
V.P=[];V.PforceZ=[];V.residual=[];V.divergence=[];V.time=t;V.approximate=true;
end
