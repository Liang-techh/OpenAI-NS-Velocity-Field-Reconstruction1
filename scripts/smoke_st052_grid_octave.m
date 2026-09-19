function smoke_st052_grid_octave(mat_path, png_path, receipt_path)
%SMOKE_ST052_GRID_OCTAVE Load the frozen ST052 MAT grid in a real Octave runtime.
%
% This is a delivery/visualization software smoke only.  It verifies that the
% MATLAB-v5 payload produced by CR-A9-071 can be loaded by GNU Octave, checks
% the documented array layout, computes a deterministic mid-plane omega_z
% diagnostic, and writes a fixed PNG.  It does not define a visual acceptance
% threshold, validate Navier-Stokes, or identify the OpenAI field.

  S = load(mat_path);
  required = {'x', 'y', 'z', 't', 'u', 'v', 'w'};
  for k = 1:numel(required)
    assert(isfield(S, required{k}), ['missing MAT variable: ', required{k}]);
  end

  x = S.x(:);
  y = S.y(:);
  z = S.z(:);
  t = S.t(:);

  assert(isequal(size(x), [33, 1]));
  assert(isequal(size(y), [33, 1]));
  assert(isequal(size(z), [33, 1]));
  assert(isequal(size(t), [5, 1]));
  assert(isequal(size(S.u), [5, 33, 33, 33]));
  assert(isequal(size(S.v), [5, 33, 33, 33]));
  assert(isequal(size(S.w), [5, 33, 33, 33]));

  expected_t = [0.25; 0.375; 0.5; 0.625; 0.75];
  assert(isequal(t, expected_t));
  assert(x(1) == -2.0 && x(end) == 2.0);
  assert(y(1) == -2.0 && y(end) == 2.0);
  assert(z(1) == -2.0 && z(end) == 2.0);

  assert(all(isfinite(S.u(:))));
  assert(all(isfinite(S.v(:))));
  assert(all(isfinite(S.w(:))));

  [~, ti] = min(abs(t - 0.5));
  [~, kz] = min(abs(z));
  U = squeeze(S.u(ti, :, :, kz));
  V = squeeze(S.v(ti, :, :, kz));
  W = squeeze(S.w(ti, :, :, kz));

  assert(isequal(size(U), [33, 33]));
  assert(isequal(size(V), [33, 33]));
  assert(isequal(size(W), [33, 33]));

  dx = x(2) - x(1);
  dy = y(2) - y(1);
  dVdx = (V(3:end, 2:end-1) - V(1:end-2, 2:end-1)) / (2.0 * dx);
  dUdy = (U(2:end-1, 3:end) - U(2:end-1, 1:end-2)) / (2.0 * dy);
  omega_z = dVdx - dUdy;

  speed = sqrt(U.^2 + V.^2 + W.^2);
  max_speed = max(speed(:));
  omega_z_rms = sqrt(mean(omega_z(:).^2));
  assert(isfinite(max_speed) && max_speed > 0.0);
  assert(isfinite(omega_z_rms));

  % Fixed, target-free visualization smoke.  The image is only proof that a
  % MATLAB-compatible consumer can load and visualize the exported numerical
  % field; no pixel/camera/source-correspondence objective is evaluated here.
  try
    graphics_toolkit('gnuplot');
  catch
    % Keep the runtime-selected toolkit if gnuplot is unavailable.
  end
  set(0, 'defaultfigurevisible', 'off');
  fig = figure('visible', 'off');
  imagesc(y(2:end-1), x(2:end-1), omega_z);
  axis xy equal tight;
  colorbar;
  xlabel('y');
  ylabel('x');
  title('ST052-M temporal child: mid-plane \omega_z at t=0.5');
  print(fig, png_path, '-dpng', '-r120');
  close(fig);

  info = dir(png_path);
  assert(~isempty(info));
  assert(info.bytes > 0);

  fid = fopen(receipt_path, 'w');
  assert(fid >= 0);
  fprintf(fid, 'task_id=CR-A9-072\n');
  fprintf(fid, 'octave_version=%s\n', version());
  fprintf(fid, 'mat_file_loaded=true\n');
  fprintf(fid, 'grid_shape=33,33,33\n');
  fprintf(fid, 'time_count=5\n');
  fprintf(fid, 'mid_time=%.17g\n', t(ti));
  fprintf(fid, 'mid_z=%.17g\n', z(kz));
  fprintf(fid, 'max_speed=%.17g\n', max_speed);
  fprintf(fid, 'omega_z_rms=%.17g\n', omega_z_rms);
  fprintf(fid, 'render_bytes=%d\n', info.bytes);
  fprintf(fid, 'actual_octave_runtime_executed=true\n');
  fprintf(fid, 'actual_matlab_runtime_executed=false\n');
  fprintf(fid, 'fixed_visualization_smoke_passed=true\n');
  fprintf(fid, 'visualization_ready=false\n');
  fprintf(fid, 'visual_correspondence_verified=false\n');
  fprintf(fid, 'pde_validated=false\n');
  fclose(fid);

  fprintf('CR-A9-072 Octave MAT load/render smoke passed.\n');
  fprintf('max_speed=%.17g omega_z_rms=%.17g render_bytes=%d\n', ...
          max_speed, omega_z_rms, info.bytes);
end
