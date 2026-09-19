function receipt=ns_cross_language_parity(dataFile,fixtureDir,outDir)
%NS_CROSS_LANGUAGE_PARITY Compare native MATLAB ST054 velocity to frozen Python output.
% Engineering/cross-language parity only; this is not PDE or visual acceptance.
if nargin<1||isempty(dataFile),dataFile=fullfile(fileparts(mfilename('fullpath')),'data','st054_models.mat');end
if nargin<2||isempty(fixtureDir),error('ns:Fixture','fixtureDir is required');end
if nargin<3||isempty(outDir),outDir=fullfile(tempdir,'ns-st054-cross-language-parity');end
if ~exist(outDir,'dir'),mkdir(outDir);end

manifest=jsondecode(fileread(fullfile(fixtureDir,'fixture_manifest.json')));
expectedTimes=[.271 .389 .503 .617 .739];
tol=1e-8;
assert(strcmp(string(manifest.schema),"st054_python_matlab_offgrid_parity_fixture_v1"));
assert(manifest.seed==9174091&&manifest.points_per_time==96);
assert(max(abs(double(manifest.times(:)).'-expectedTimes))<1e-15);
assert(abs(double(manifest.engineering_atol)-tol)<eps(tol));
assert(~logical(manifest.pde_validated));
assert(~logical(manifest.visual_correspondence_verified));
assert(~logical(manifest.openai_field_identified));

S=load(dataFile);if ~iscell(S.models),S.models=num2cell(S.models);end
assert(strcmp(string(S.schema),"ns_matlab_spectral_v1"));
assert(strcmp(string(S.source_commit),"c77492a48e9c0f13d4d51244987c57c28519409b"));
assert(~logical(S.pde_validated));
modelIDs=["ST054-Q2","ST054-M3"];
maxErrors=zeros(1,numel(modelIDs));
rowsPerModel=zeros(1,numel(modelIDs));

for k=1:numel(modelIDs)
    modelID=modelIDs(k);m=[];
    for j=1:numel(S.models)
        if strcmp(string(S.models{j}.id),modelID),m=S.models{j};break;end
    end
    assert(~isempty(m),'ns:Model','Published model missing from MAT payload.');
    assert(~m.pde_validated);
    filename=lower(strrep(char(modelID),'-','_'));
    filename=[filename '.csv'];
    table=readmatrix(fullfile(fixtureDir,filename));
    assert(isequal(size(table),[96*numel(expectedTimes),7]),'ns:Fixture','Unexpected fixture shape.');
    assert(all(isfinite(table),'all'),'ns:Fixture','Fixture contains non-finite values.');
    assert(all(hypot(table(:,1),table(:,2))<1.70+1e-14));
    assert(all(abs(table(:,3))<1.70+1e-14));
    rowsPerModel(k)=size(table,1);
    localMax=0;
    for ti=1:numel(expectedTimes)
        time=expectedTimes(ti);
        ix=abs(table(:,4)-time)<1e-15;
        assert(nnz(ix)==96,'ns:Fixture','Unexpected points-per-time count.');
        actual=ns_evaluate(m,table(ix,1:3),time);
        expected=table(ix,5:7);
        localMax=max(localMax,max(abs(actual-expected),[],'all'));
    end
    assert(localMax<=tol,'ns:Parity','Python/native-MATLAB ST054 velocity parity failed.');
    maxErrors(k)=localMax;
end

receipt=struct(...
    'schema','st054_python_native_matlab_offgrid_parity_v1',...
    'matlab_version',version,...
    'model_ids',{cellstr(modelIDs)},...
    'rows_per_model',rowsPerModel,...
    'times',expectedTimes,...
    'engineering_atol',tol,...
    'max_component_error',maxErrors,...
    'parity_pass',all(maxErrors<=tol),...
    'scope','cross_language_velocity_parity_only_not_pde_or_visual_acceptance',...
    'pde_validated',false,...
    'visual_correspondence_verified',false,...
    'source_correspondence_verified',false,...
    'paper_exact',false,...
    'openai_field_identified',false,...
    'blowup_proved',false);
fid=fopen(fullfile(outDir,'cross_language_parity_receipt.json'),'w');
if fid<0,error('ns:Receipt','Could not open receipt output.');end
clean=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(receipt),'char');
fprintf('%s\n',jsonencode(receipt));
end
