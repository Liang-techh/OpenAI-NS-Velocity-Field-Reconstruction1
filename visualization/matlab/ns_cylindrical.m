function D = ns_cylindrical(m,s,z,t,tensor)
%NS_CYLINDRICAL Native continuous-time field and analytic full NS residual.
% tensor=true returns [numel(s),numel(z)], avoiding giant 3-D basis tensors.
% tensor=false returns paired point values; t is a scalar, NOT a frame number.
if nargin<5,tensor=false;end
validateattributes(t,{'numeric'},{'scalar','real','finite','>=',m.tmin,'<=',m.tmax});
s=s(:);z=z(:);
if any(~isfinite(s)|s<0)||any(~isfinite(z)),error('ns:Points','Invalid cylindrical points.');end
if ~tensor&&numel(s)~=numel(z),error('ns:Shape','Paired arrays must have equal lengths.');end
nr=size(m.F,1);nz=size(m.F,2);nt=size(m.F,3);
R=ns_basis(s,'r',0,nr);ZF=ns_basis(z,'z',1,nz);ZG=ns_basis(z,'z',0,nz);
x=4*t-2;T=zeros(nt,1);Td=T;T(1)=1;
if nt>1,T(2)=x;Td(2)=1;end
for k=3:nt,T(k)=2*x*T(k-1)-T(k-2);Td(k)=2*T(k-1)+2*x*Td(k-1)-Td(k-2);end
K=reshape(reshape(m.F,[],nt)*T,nr,nz);Kt=reshape(reshape(m.F,[],nt)*(4*Td),nr,nz);
G=reshape(reshape(m.G,[],nt)*T,nr,nz);Gt=reshape(reshape(m.G,[],nt)*(4*Td),nr,nz);
P=reshape(reshape(m.P,[],nt)*T,nr,nz);
F00=ev(K,ZF,0,0);F10=ev(K,ZF,1,0);F20=ev(K,ZF,2,0);
F01=ev(K,ZF,0,1);F11=ev(K,ZF,1,1);F02=ev(K,ZF,0,2);
D.A=-F01;D.As=-F11;D.Az=-F02;D.At=-ev(Kt,ZF,0,1);
D.AL=-8*F11-4*s.*ev(K,ZF,2,1)-ev(K,ZF,0,3);
D.C=2*F00+2*s.*F10;D.Cs=4*F10+2*s.*F20;D.Cz=2*F01+2*s.*F11;
D.Ct=2*ev(Kt,ZF,0,0)+2*s.*ev(Kt,ZF,1,0);
D.CL=16*F10+32*s.*F20+8*s.^2.*ev(K,ZF,3,0)+2*F02+2*s.*ev(K,ZF,1,2);
D.B=ev(G,ZG,0,0);D.Bs=ev(G,ZG,1,0);D.Bz=ev(G,ZG,0,1);
D.Bt=ev(Gt,ZG,0,0);D.BL=8*D.Bs+4*s.*ev(G,ZG,2,0)+ev(G,ZG,0,2);
D.P=ev(P,ZG,0,0);D.Ps=ev(P,ZG,1,0);D.Pz=ev(P,ZG,0,1);
% Fixed two-parameter forcing: identical to source spacetime.force.
[br,br1]=bump01(s/4);[bz,bz1]=bump01(z.^2/4);
if tensor,bz=bz.';bz1=bz1.';zz=z.';else,zz=z;end
b=br.*bz;bs=br1.*bz/2;bzv=br.*bz1.*zz/2;
g=bump01((2*t-1)^2);a=m.force(1);c=m.force(2);
D.fa=g*(-a*(b+zz.*bzv));D.fb=g*c*(b+s.*bs/2);D.fc=g*a*zz.*(2*b+s.*bs);
D.Ra=D.At+D.A.^2+2*s.*D.A.*D.As+D.C.*D.Az-D.B.^2+2*D.Ps-m.nu*D.AL-D.fa;
D.Rb=D.Bt+2*D.A.*D.B+2*s.*D.A.*D.Bs+D.C.*D.Bz-m.nu*D.BL-D.fb;
D.Rc=D.Ct+2*s.*D.A.*D.Cs+D.C.*D.Cz+D.Pz-m.nu*D.CL-D.fc;
D.divergence=2*D.A+2*s.*D.As+D.Cz;
    function a=ev(co,Z,i,j)
        a=R{i+1}*co;
        if tensor,a=a*Z{j+1}.';else,a=sum(a.*Z{j+1},2);end
    end
end
function [b,b1]=bump01(q)
d=ones(size(q));ix=q<1;d(ix)=1-q(ix);b=zeros(size(q));b(ix)=exp(1-1./d(ix));b1=-b./d.^2;
end
