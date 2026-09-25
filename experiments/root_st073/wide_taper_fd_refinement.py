"""Refine the widened physical-lift hotspot."""
from coupled_moment_fd_refinement import run


if __name__ == '__main__':
    run(slice_filename='wide_taper_five_moment_slice.json',
        output_name='wide_taper_fd_refinement.json', X=2.975)
