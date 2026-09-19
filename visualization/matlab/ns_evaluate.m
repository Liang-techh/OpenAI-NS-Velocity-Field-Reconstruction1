function [u,p,omega,res,pressureForce] = ns_evaluate(m,points,t)
%NS_EVALUATE Exact finite-basis evaluation in Cartesian coordinates (N-by-3).
% No time interpolation. Chunking prevents large temporary allocations.
validateattributes(points,{'numeric'},{'2d','ncols',3,'real','finite'});
validateattributes(t,{'numeric'},{'scalar','real','finite','>=',m.tmin,'<=',m.tmax});
n=size(points,1);u=zeros(n,3);p=zeros(n,1);omega=u;res=u;pressureForce=u;
for first=1:1024:n
    ix=first:min(n,first+1023);x=points(ix,1);y=points(ix,2);z=points(ix,3);s=x.^2+y.^2;
    D=ns_cylindrical(m,s,z,t,false);
    u(ix,:)=[x.*D.A-y.*D.B,y.*D.A+x.*D.B,D.C];p(ix)=D.P;
    omega(ix,:)=[-x.*D.Bz-y.*(D.Az-2*D.Cs),-y.*D.Bz+x.*(D.Az-2*D.Cs),2*(D.B+s.*D.Bs)];
    res(ix,:)=[x.*D.Ra-y.*D.Rb,y.*D.Ra+x.*D.Rb,D.Rc];
    pressureForce(ix,:)=[-2*x.*D.Ps,-2*y.*D.Ps,-D.Pz];
end
end
