function D=ns_import_grid(filename)
%NS_IMPORT_GRID Import PR665-compatible velocity grids WITHOUT claiming ST054 identity.
% Required x,y,z,t,u,v,w; component layout [time,x,y,z]. Never guess permutations.
S=load(filename);keys={'x','y','z','t','u','v','w'};
for k=1:numel(keys),if ~isfield(S,keys{k}),error('ns:GridSchema','Missing %s.',keys{k});end,end
for k=1:4
    v=double(S.(keys{k}));if ~isvector(v)||numel(v)<2||any(~isfinite(v))||any(diff(v(:))<=0)
        error('ns:GridSchema','Coordinates/time must be finite strictly increasing vectors.');end
    S.(keys{k})=v(:).';
end
shape=[numel(S.t),numel(S.x),numel(S.y),numel(S.z)];
for k=5:7
    v=double(S.(keys{k}));if ~isequal(size(v),shape)||any(~isfinite(v(:)))
        error('ns:GridSchema','Expected component layout [time,x,y,z].');end
    S.(keys{k})=v;
end
[~,name]=fileparts(filename);m=S;m.mode='sampled';m.id=['Imported: ' name];m.tmin=S.t(1);m.tmax=S.t(end);m.pde_validated=0;
m.original_raw_sha256='Unverified grid import; consult the original export receipt';
m.note='Linear space/time interpolation; curl from grid differences. Pressure/full NS residual unavailable.';
D=struct('schema','ns_matlab_imported_grid_v1','models',{{m}},'pde_validated',0);
end
