"""Refine the entrance hotspot after curvature optimization."""
from coupled_moment_fd_refinement import run


if __name__ == '__main__':
    run(slice_filename='wide_taper_curvature_optimize.json',
        output_name='wide_taper_curvature_fd_refinement.json',
        X=1.01, eta=.3)
