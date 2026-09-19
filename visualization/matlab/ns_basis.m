function B = ns_basis(v,kind,odd,n)
%NS_BASIS Values and derivatives (orders 0:3) of the pinned ST054 1-D basis.
% No Symbolic Toolbox. Columns 1:7 are Legendre*bump; 8:9 are b/(1-q)^2,4.
% Axial columns 10:12 have the opposite parity. Derivatives are in s or z.
v=v(:); B=cell(1,4);
if strcmp(kind,'r')
    q=v/4; base=bump(q);
    for d=0:3, base(:,d+1)=base(:,d+1)/4^d; end
    L=legendreJets(v/2-1,6);
    for d=0:3
        V=zeros(numel(v),n);
        for j=0:min(n,7)-1
            for k=0:d
                V(:,j+1)=V(:,j+1)+nchoosek(d,k)*base(:,k+1).*L(:,j+1,d-k+1)/2^(d-k);
            end
        end
        for j=7:n-1
            R=rational(q,2*(j-6));V(:,j+1)=R(:,d+1)/4^d;
        end
        B{d+1}=V;
    end
else
    if n~=12, error('ns:Basis','This export expects exactly 12 axial columns.'); end
    base=zchain(bump(v.^2/4),v);L=legendreJets(v/2,13);
    for d=0:3
        V=zeros(numel(v),12);
        for j=0:6
            for k=0:d
                V(:,j+1)=V(:,j+1)+nchoosek(d,k)*base(:,k+1).*L(:,2*j+odd+1,d-k+1)/2^(d-k);
            end
        end
        for j=7:8
            R=zchain(rational(v.^2/4,2*(j-6)),v);
            if odd
                V(:,j+1)=v/2.*R(:,d+1);
                if d>0, V(:,j+1)=V(:,j+1)+d/2*R(:,d); end
            else, V(:,j+1)=R(:,d+1);
            end
        end
        for j=0:2
            for k=0:d
                V(:,j+10)=V(:,j+10)+nchoosek(d,k)*base(:,k+1).*L(:,2*j+(1-odd)+1,d-k+1)/2^(d-k);
            end
        end
        B{d+1}=V;
    end
end
end
function out=bump(q)
d=ones(size(q));inside=q<1;d(inside)=1-q(inside);
y=1./d; b=zeros(size(q));b(inside)=exp(1-y(inside));
out=[b,-b.*y.^2,b.*(y.^4-2*y.^3),b.*(-y.^6+6*y.^5-6*y.^4)];
end
function out=rational(q,k)
d=ones(size(q));inside=q<1;d(inside)=1-q(inside);y=1./d;
b=zeros(size(q));b(inside)=exp(1-y(inside));yc=min(y,1000);
c=zeros(1,k+8);c(k+1)=1;out=zeros(numel(q),4);
for d=0:3
    v=zeros(size(q));
    for pow=0:numel(c)-1, if c(pow+1)~=0, v=v+c(pow+1)*yc.^pow; end, end
    out(:,d+1)=b.*v;out(y>1000,d+1)=0;
    next=zeros(size(c));
    for pow=0:numel(c)-3
        next(pow+2)=next(pow+2)+pow*c(pow+1);
        next(pow+3)=next(pow+3)-c(pow+1);
    end
    c=next;
end
end
function Z=zchain(b,z)
Z=[b(:,1),b(:,2).*z/2,b(:,3).*z.^2/4+b(:,2)/2,...
   b(:,4).*z.^3/8+3*b(:,3).*z/4];
end
function L=legendreJets(x,n)
L=zeros(numel(x),n+1,4);L(:,1,1)=1;
if n==0,return,end
L(:,2,1)=x;L(:,2,2)=1;
for j=2:n
    for d=0:3
        v=x.*L(:,j,d+1);
        if d>0, v=v+d*L(:,j,d);end
        L(:,j+1,d+1)=((2*j-1)*v-(j-1)*L(:,j-1,d+1))/j;
    end
end
end
