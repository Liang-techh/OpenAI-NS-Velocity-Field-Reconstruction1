# Uniform C1 inverse for the fresh five-bump defect family

The directed full-support bump weights and fresh R110 five-row defect family now admit a certified implicit C1 coefficient function h(Z) on [0.49,0.51], for fixed accepted construction parameters. The proof gives existence and uniqueness inside a uniform coefficient box, and encloses h and h_Z. It does not yet install or independently validate the repaired physical reference-annulus field.

## Map and interval inputs

The coefficient order is (c1,c2,xi1,xi2,xi3). The paper's map is Lh+Q(h,Am), with the two axial and three angular disjoint bumps and the mixed/quadratic terms in equation10.8. All linear and quadratic weights use the previously directed integrals of the actual normalized bump family. Exact full-support normalized masses are1. No midpoint of a source defect or map coefficient replaces its interval.

The amplitude Am is exp(13.4)/(1+Z^2), from the same fixed logPstar14 datum. The defect values and first axial derivatives retain the full controlled source, flat-kernel and restoration errors.

## Uniform contraction

Choose a fixed preconditioner R from the inverse of the midpoint linear matrix. This is only a numerical choice of preconditioner; R need not equal the exact inverse. Directed bounds for I-RL prove that R is nonsingular. The true fixed-point map is

    Phi(h) = (I-RL)h - R(Q(h,Am)+d).

On the box |c_i|<=2e-14 and |xi_i|<=1e-32, its interval image lies strictly inside the box. The weighted infinity norm of I-R(L+DQ) is at most0.00271269494639, strictly below1. Banach contraction therefore gives one true solution for each accepted axial center. The bounds are uniform in the center, so they are not isolated point solves.

Twelve successive interval intersections tighten the solution enclosure without projecting inputs. The resulting axial controls are approximately c1=-1.00e-14 and c2=7.76e-15. Angular controls are of order1e-36, with their exact intervals stored in the receipt.

## First axial derivative

Implicit differentiation gives

    (L+DQ) h_Z = -d_Z - partial_Z Q.

Only the amplitude-dependent row4 contributes partial_Z Q, through (Am^-2)_Z sum(GG_i c_i^2). A directed Neumann contraction bounds this linear solve, then interval intersections tighten h_Z. Smoothness and invertibility provide an implicit C1 family. Raw residual intervals containing zero are recorded only as diagnostics; they are not used as the existence proof.

An independent prescribed nonlinear solution fixture checks all ten value/derivative coefficients against the certified boxes. A data set outside the initial box is correctly rejected. A read-only mathematical review covers the contraction and derivative argument.

## Scope and next construction step

This proves a local uniform five-bump map inverse and supplies coefficient enclosures. The actual repaired field still needs a callable interval evaluator for normalized bumps and partial cumulative integrals, preserving the same P0, F0 and source moments. Verify its terminal moment identities, divergence recovery and annular stress before extending to flatten/heat matching. Whole-axis coverage, higher axial derivatives, source parameter errors, temporal recursion, flat remainder and oscillatory correction remain incomplete.
