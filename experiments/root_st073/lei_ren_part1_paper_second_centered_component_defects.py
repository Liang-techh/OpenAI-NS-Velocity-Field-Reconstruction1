"""Second-Z assembly of the same labeled centered five-moment defects."""
import mpmath as mp
from lei_ren_part1_paper_axial_second_jet import AxialSecondJet
from lei_ren_part1_paper_centered_component_defects import CenteredComponentDefects,KERNEL_SPECS
from lei_ren_part1_paper_flat_shape_component import compose_flat_shape_defect


class SecondCenteredComponentDefects(CenteredComponentDefects):
    branch='second_centered_component_defects'

    def _dual(self,value=0,tangent=0,second=0):
        if isinstance(value,AxialSecondJet):return value
        return AxialSecondJet(self._jet(value),self._jet(tangent),self._jet(second))

    def _source(self,Z):
        with mp.workdps(self.work_precision):
            z=mp.mpf(Z);key=mp.nstr(z,self.work_precision)
            if key in self._z_cache:return self._z_cache[key]
            source=self.switches.evaluate_R(110,z)
            fields=source['second_jet_fields'];moments=source['second_jet_moments']
            for name in ('F','Uz','P'):
                if not isinstance(fields[name],AxialSecondJet):
                    raise TypeError('Second centered defects require all three source slots')
            if any(not isinstance(v,AxialSecondJet) for v in moments.values()):
                raise TypeError('Second centered moments cannot come from first-only data')
            for name in ('P0','P0_Z','P0_ZZ'):
                if name not in source:raise KeyError(f'Explicit common source {name} required')
            result=dict(Z=z,switch=source,F110=fields['F'],Uz110=fields['Uz'],
                        moments=moments,P110=fields['P'],
                        P0=self._dual(source['P0'],source['P0_Z'],source['P0_ZZ']))
            self._z_cache[key]=result
            return result

    def _kernel_records(self,B):
        return {name:compose_flat_shape_defect(k,m,B.value,B.tangent,self.T,
                    B_ZZ=B.second,precision=self.precision,order=self.order,window=self.window)
                for name,(k,m) in KERNEL_SPECS.items()}

    @staticmethod
    def _flat_term(record,n):
        return AxialSecondJet(record['value_terms'][n],record['tangent_terms'][n],record['second_terms'][n])

    def evaluate(self,Z):
        with mp.workdps(self.work_precision):
            result=super().evaluate(Z)
            result['d_ZZ']={row:v.second for row,v in result['d_dual'].items()}
            result['parts_ZZ']={label:{row:v.second for row,v in rows.items()}
                                for label,rows in result['parts_dual'].items()}
            source=self._source(Z);z=self._dual(source['Z'],1)
            Am=mp.exp(self.logPstar-mp.mpf('.6'))/(1+z*z)
            uref=mp.exp(-self.logC)/(1+z*z)
            B=(mp.sqrt(220)*source['F110']).log()+self.logC+(1+z*z).log()
            result.update(P0_ZZ=source['P0'].second,P110_ZZ=source['P110'].second,
                          F110_ZZ=source['F110'].second,Uz110_ZZ=source['Uz110'].second,
                          Am_ZZ=Am.second,uref110_ZZ=uref.second,B_ZZ=B.second,
                          g_ZZ=source['Uz110'].second,
                          d_second_jet=result['d_dual'],parts_second_jet=result['parts_dual'])
            result['P0_receipt']['P0_ZZ']=source['P0'].second
            return result

    __call__=evaluate

    def metadata(self):
        result=super().metadata()
        result.update(automatic_second_Z_derivative=True,second_Z_remainder_enclosed=False,
                      source_second_Z_required=True,finite_Z_difference_used=False)
        return result
