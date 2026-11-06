"""Independent LTE population check and conservative source-only screening."""
import numpy as np


def independent_ion_populations(partition,atmosphere):
    """Log-space Saha recursion, separate from native multiplicative recursion."""
    import lightweaver.constants as c
    from scipy.special import logsumexp
    T,ne,nH=atmosphere.temperature,atmosphere.ne,atmosphere.nHTot
    if T.min()<partition.Tpf[0] or T.max()>partition.Tpf[-1]:raise ValueError('Partition domain does not cover column')
    logU=np.array([np.interp(T,partition.Tpf,row) for row in partition.pf])
    stages=len(partition.ionPot);relative=np.zeros((stages,len(T)))
    logpref=np.log(2.)+1.5*np.log(2*np.pi*c.MElectron*c.KBoltzmann*T/c.HPlanck**2)-np.log(ne)
    for stage in range(1,stages):
        relative[stage]=relative[stage-1]+logpref+logU[stage]-logU[stage-1]-partition.ionPot[stage-1]/(c.KBoltzmann*T)
    return np.exp(relative-logsumexp(relative,axis=0))*(partition.abundance*nH)[None,:]


def source_eligibility(item,*,changed_on_main=False):
    if changed_on_main:return 'MAIN_CATALOG_CHANGED'
    if item['treatment']not in ('not_in_trace_subset','trace_FeI'):return item['treatment'].upper()
    if item['exact_identity_matches']!=1:return 'SOURCE_IDENTITY_NOT_UNIQUE'
    if int(item['canonical_row']['hfs_n_components'])!=1:return 'COMBINED_COMPONENTS'
    if not item['canonical_gf_matches_any_source']:return 'CANONICAL_SOURCE_GF_DIFFER'
    row=item['source_rows'][0]
    gamma,loggamma=float(row['turbospectrum_rad']),float(row['rad'])
    if not np.isfinite([gamma,loggamma]).all() or gamma<=0 or abs(np.log10(gamma)-loggamma)>.005:
        return 'NATURAL_WIDTH_UNQUALIFIED'
    return 'SOURCE_ELIGIBLE_FOR_NATURAL_ONLY_SCREEN'
