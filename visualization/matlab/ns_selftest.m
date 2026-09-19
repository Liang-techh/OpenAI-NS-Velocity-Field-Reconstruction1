function receipt=ns_selftest(dataFile,guiSmoke,outDir)
%NS_SELFTEST Actual native evaluator and optional UI smoke. Not NS acceptance.
if nargin<1||isempty(dataFile),dataFile=fullfile(fileparts(mfilename('fullpath')),'data','st054_models.mat');end
if nargin<2,guiSmoke=false;end
if nargin<3,outDir=fullfile(fileparts(mfilename('fullpath')),'tests','output');end
if ~exist(outDir,'dir'),mkdir(outDir);end
S=load(dataFile);if ~iscell(S.models),S.models=num2cell(S.models);end
receipt=struct('matlab_version',version,'numerical_pass',false,'gui_pass',false,...
    'pde_validated',false,'scope','Native MATLAB software/reference checks, not new PDE acceptance');
errors=[];
for k=1:numel(S.models)
    m=S.models{k};e=zeros(1,3);
    assert(~m.pde_validated);assert(all(size(m.F)==[9 12 8]));
    for j=1:numel(m.ref_times)
        [u,p,w,r]=ns_evaluate(m,m.ref_points,m.ref_times(j)); %#ok<ASGLU>
        ru=squeeze(m.ref_velocity(j,:,:));rp=m.ref_pressure(j,:).';rr=squeeze(m.ref_residual(j,:,:));
        e=max(e,[max(abs(u(:)-ru(:))),max(abs(p(:)-rp(:))),max(abs(r(:)-rr(:)))]);
    end
    assert(e(1)<1e-8&&e(2)<1e-8&&e(3)<1e-8,'ns:Parity','Source parity failed.');
    % Tensor volume and paired evaluation must agree at actual off-axis nodes.
    V=ns_volume(m,.46317,17);ii=[23 141 986 1705 3666];
    [u,p,om,rr]=ns_evaluate(m,[V.X(ii)' V.Y(ii)' V.Z(ii)'],.46317);
    assert(max(abs(u(:,1)-V.U(ii)'))<1e-9);
    assert(max(abs(p-V.P(ii)'))<1e-9);
    assert(max(abs(sqrt(sum(om.^2,2))-V.vorticity(ii)'))<1e-9);
    assert(max(abs(sqrt(sum(rr.^2,2))-V.residual(ii)'))<1e-9);
    assert(max(abs(V.divergence(:)))<1e-9);
    [u,p]=ns_evaluate(m,[2 0 0;0 0 2;0 0 -2;3 0 1],.5);
    assert(all(u(:)==0)&&all(p(:)==0),'Compact support failed.');
    % Independent finite differences of curl and pressure gradients.
    xyz=[.2 .1 .17;.31 -.12 -.5;0 0 .12];h=1e-5;J=zeros(3,3,3);pg=zeros(3,3);
    for d=1:3
        delta=zeros(1,3);delta(d)=h;
        [up,pp]=ns_evaluate(m,xyz+delta,.517);[um,pm]=ns_evaluate(m,xyz-delta,.517);
        J(:,:,d)=(up-um)/(2*h);pg(:,d)=(pp-pm)/(2*h);
    end
    wf=[J(:,3,2)-J(:,2,3),J(:,1,3)-J(:,3,1),J(:,2,1)-J(:,1,2)];
    [~,~,om,~,pf]=ns_evaluate(m,xyz,.517);
    assert(max(abs(wf(:)-om(:)))<1e-6);assert(max(abs(pg(:)+pf(:)))<1e-6);
    a=.39;Q=[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
    u=ns_evaluate(m,xyz,.517);rot=ns_evaluate(m,xyz*Q.',.517);
    assert(max(abs(rot(:)-reshape(u*Q.',[],1)))<1e-9);
    rejected=false;try,ns_evaluate(m,xyz,.751);catch,rejected=true;end;assert(rejected);
    errors=[errors;e]; %#ok<AGROW>
end
receipt.max_velocity_pressure_residual_errors=max(errors,[],1);
% PR665 layout adapter smoke using a labelled unit-test grid (not an NS candidate).
x=linspace(-2,2,7);y=linspace(-2,2,9);z=linspace(-2,2,11);t=[.25,.5,.75];
[X,Y,Z]=ndgrid(x,y,z);u=zeros(3,7,9,11);v=u;w=u;
for i=1:3,u(i,:,:,:)=-Y;v(i,:,:,:)=X;w(i,:,:,:)=t(i)*Z;end
sampleFile=fullfile(outDir,'adapter_unit_test.mat');save(sampleFile,'x','y','z','t','u','v','w');
D=ns_import_grid(sampleFile);G=ns_volume(D.models{1},.413,13);
assert(max(abs(G.U(:)+G.Y(:)))<1e-12);assert(max(abs(G.W(:)-.413*G.Z(:)))<1e-12);
assert(isempty(G.residual)&&isempty(G.P));delete(sampleFile);
receipt.numerical_pass=true;
if guiSmoke
    a=ns_explorer(dataFile,'Visible','off');cleanup=onCleanup(@()a.Close());
    assert(isempty(getappdata(a.Figure,'LastRenderError')));
    a.SetTime(.63127);a.SetControl('count',80);a.SetControl('slice',.3);
    a.SetScalar('Full NS residual');a.SetMode('Streamlines + surface');
    assert(isempty(getappdata(a.Figure,'LastRenderError')));
    a.SetCandidate('ST054-M3');a.SetMode('Velocity arrows');
    assert(isempty(getappdata(a.Figure,'LastRenderError')));
    a.SetCandidate('ST054-Q2');a.SetMode('Streamlines');a.SetControl('count',200);
    a.SetScalar('Vorticity');a.SetTime(.5);
    assert(isempty(getappdata(a.Figure,'LastRenderError')));
    a.Figure.Visible='on';drawnow;pause(1);
    exportapp(a.Figure,fullfile(outDir,'matlab_explorer.png'));
    receipt.gui_pass=true;receipt.snapshot=a.Snapshot();
    clear cleanup
end
fid=fopen(fullfile(outDir,'matlab_test_receipt.json'),'w');clean=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(receipt),'char');fprintf('%s\n',jsonencode(receipt));
end
