function ns_callback_test(dataFile,outDir)
%NS_CALLBACK_TEST Exercise actual time-slider and playback callback functions.
a=ns_explorer(dataFile,'Visible','off');clean=onCleanup(@()a.Close());
cb=a.TimeSlider.ValueChangingFcn;cb(a.TimeSlider,struct('Value',.41937));
a.TimeSlider.Value=.41937;
cb=a.TimeSlider.ValueChangedFcn;cb(a.TimeSlider,struct('Value',.41937));
st=a.Snapshot();assert(abs(st.t-.41937)<1e-12);
assert(isempty(getappdata(a.Figure,'LastRenderError')));
cb=a.PlayButton.ButtonPushedFcn;cb(a.PlayButton,[]);pause(.7);cb(a.PlayButton,[]);
st=a.Snapshot();assert(st.t>.41937);
assert(isempty(getappdata(a.Figure,'LastRenderError')));
receipt=struct('slider_callbacks_pass',true,'playback_timer_pass',true,...
    'last_time',st.t,'pde_validated',false);
if ~exist(outDir,'dir'),mkdir(outDir);end
fid=fopen(fullfile(outDir,'callback_receipt.json'),'w');guard=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(receipt),'char');fprintf('%s\n',jsonencode(receipt));
end
