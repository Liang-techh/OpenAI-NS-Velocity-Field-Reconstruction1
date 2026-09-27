"""Eliminate linear momentum controls while preserving moment equalities.

All residual rows supplied here must already carry sqrt(physical weights).
The constraints act on unscaled physical control coefficients. This small
operator is independent of how the physical moment rows are assembled.
"""
import numpy as np


class ConstrainedTangent:
    def __init__(self, weighted_design, moment_rows, rcond=1e-10):
        design = np.asarray(weighted_design,float)
        rows = np.asarray(moment_rows,float)
        if design.ndim != 2 or rows.ndim != 2 or design.shape[1] != rows.shape[1]:
            raise ValueError('Design and moment rows must share control columns')
        self.rows = rows
        self.scales = np.maximum(np.linalg.norm(design,axis=0),1e-30)
        self.design = design/self.scales
        scaled_rows = rows/self.scales
        row_scales = np.maximum(np.linalg.norm(scaled_rows,axis=1),1e-30)
        normalized_rows = scaled_rows/row_scales[:,None]
        u,s,vh = np.linalg.svd(normalized_rows,full_matrices=True)
        rank = int(np.sum(s > rcond*s[0])) if len(s) else 0
        if rank != len(rows):
            raise ValueError('Dependent moment rows require explicit target compatibility handling')
        self.particular = ((vh[:rank].T/s[:rank])@u[:,:rank].T)/row_scales
        self.nullspace = vh[rank:].T
        reduced = self.design@self.nullspace
        ru,rs,rvh = np.linalg.svd(reduced,full_matrices=False)
        keep = rs > rcond*rs[0] if len(rs) else np.zeros(0,dtype=bool)
        self.Q = ru[:,keep]
        self.inverse_null = (rvh[keep].T/rs[keep])@self.Q.T
        self.target_response = self.design@self.particular
        self.metadata = dict(moment_rank=rank,null_design_rank=int(keep.sum()),
                             moment_singular_values=s.tolist(),
                             null_design_singular_values=rs.tolist(),rcond=rcond)

    def solve(self, weighted_residual, moment_target, residual_jacobian=None,
              target_jacobian=None):
        residual = np.asarray(weighted_residual,float)
        target = np.asarray(moment_target,float)
        particular = self.particular@target
        shifted = residual+self.design@particular
        controls_scaled = particular-self.nullspace@(self.inverse_null@shifted)
        controls = controls_scaled/self.scales
        corrected = shifted-self.Q@(self.Q.T@shifted)
        result = dict(controls=controls,residual=corrected,
                      moment_error=self.rows@controls-target,
                      objective=float(corrected@corrected))
        if residual_jacobian is not None:
            jac = np.asarray(residual_jacobian,float)
            if target_jacobian is not None:
                jac = jac+self.target_response@np.asarray(target_jacobian,float)
            result['objective_gradient'] = 2*jac.T@corrected
        return result


def quadratic_moment_target(baseline, forms, x):
    """Target E control = -(baseline + x^T H_k x), and its Jacobian."""
    forms = np.asarray(forms,float)
    forms = .5*(forms+forms.swapaxes(-1,-2))
    x = np.asarray(x,float)
    return (-np.asarray(baseline,float)-np.einsum('i,kij,j->k',x,forms,x),
            -2*np.einsum('kij,j->ki',forms,x))
