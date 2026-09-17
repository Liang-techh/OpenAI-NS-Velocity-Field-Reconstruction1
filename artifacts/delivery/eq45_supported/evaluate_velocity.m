% Interpolated sampled candidate; not exact off-grid evaluation.
S = load('velocity.mat');
Fu = griddedInterpolant({S.times,S.x,S.y,S.z},S.u,'linear','none');
Fv = griddedInterpolant({S.times,S.x,S.y,S.z},S.v,'linear','none');
Fw = griddedInterpolant({S.times,S.x,S.y,S.z},S.w,'linear','none');
velocity = @(x,y,z,t) [Fu(t,x,y,z),Fv(t,x,y,z),Fw(t,x,y,z)];
disp(velocity(0.1,0,0.1,0.5));
