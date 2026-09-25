"""Physical full-momentum comparison for the width-0.05 moment lift."""
from coupled_moment_physical_screen import run


if __name__ == '__main__':
    run(slice_filename='wide_taper_five_moment_slice.json',
        output_name='wide_taper_physical_screen.json',
        sample_X=(1.01, 1.025, 1.1, 2., 2.975, 2.99))
