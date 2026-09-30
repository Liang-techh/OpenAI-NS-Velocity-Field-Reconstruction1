# Shared continuous axial bump basis

The provider defines one smooth MP bump for point values, first/second derivatives, weighted full integrals, partial primitives, matrix entries and the energy Gram integral. Support endpoints follow the function definition.

160/200-digit matrix refinement differs by about 1.87e-161. The continuous matrix differs from the previous binary64-node canonical matrix by about 2.21e-19. Thus exact cancellation in the stored quadrature graph does not prove cancellation for this continuous field. Primitive derivative replay refines from 1.01e-17 to 6.30e-19.

Run `python experiments/root_st073/lei_ren_part1_paper_continuous_axial_basis.py`.

Next: use this provider for the actual coefficient solve, values and accumulated means together, including Z derivatives. Uniform quadrature enclosures, incoming primitive provenance, pulse energy, and global finite energy remain open. This provider is not yet installed in the global profile.
