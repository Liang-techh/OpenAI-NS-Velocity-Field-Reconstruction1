"""Second-Z propagation through the finite Section 9.23 comparison."""
from lei_ren_part1_paper_pressure_width_axial_comparison import AxialPressureWidthComparison
from lei_ren_part1_paper_axial_second_jet import AxialSecondJet
from lei_ren_part1_paper_axial_dual import AxialDual


class SecondAxialPressureWidthComparison(AxialPressureWidthComparison):
    automatic_Z_second=True
    required_component_Z_depth=3
    required_component_Z_row_length=4

    def _dual(self,value=0,tangent=0,second=0):
        if isinstance(value,AxialSecondJet):
            if tangent!=0 or second!=0:raise ValueError('Do not overwrite second-jet derivatives')
            return value
        if isinstance(value,AxialDual) or isinstance(tangent,AxialDual):
            raise TypeError('First-Z data cannot supply an unknown second derivative')
        return AxialSecondJet(self._base_jet(value),self._base_jet(tangent),self._base_jet(second))

    def jet(self,value=0):
        if isinstance(value,AxialSecondJet):
            if value.orders!=(self.pressure_order,self.width_order):
                raise ValueError('Second-jet orders do not match comparison')
            return value
        return self._dual(value)

    def _differentiate_rows(self,rows):
        result=[]
        for n,row in enumerate(rows):
            if len(row)<4:
                raise ValueError(f'Second-Z comparison requires four Taylor coefficients in radial row {n}')
            values=[c.value if isinstance(c,AxialSecondJet) else c for c in row]
            output=[]
            for k,value in enumerate(values):
                first=(k+1)*values[k+1] if k+1<len(values) else 0
                second=(k+1)*(k+2)*values[k+2] if k+2<len(values) else 0
                output.append(self._dual(value,first,second))
            result.append(output)
        return result

    def metadata(self):
        result=super().metadata()
        result.update(automatic_Z_second=True,axial_dual_ring=False,
                      axial_second_jet_ring=True,second_Z_remainder_enclosed=False,
                      second_of_first_Z_driver_uses_third_core_Z_row=True)
        return result
