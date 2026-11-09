"""Finite reference/support accounting; diagnostic only, never adoption."""
import numpy as np
from pipeline.bride_continuum_audit import trace_fit


def finite_accounting(corners,key):
    """Path-averaged 2x2 changes; row is support, column is denominator.

    Dynamic fits can change clipping in each corner. These finite effects are
    correlated and order dependent before averaging, not physical error terms.
    A refused corner is retained and makes the decomposition unavailable.
    """
    values=[corners[k].get(key) for k in ('before_before','before_after','after_before','after_after')]
    if any(v is None for v in values):return dict(available=False,reason='REFUSED_CORNER')
    if not np.isfinite(values).all():raise ValueError('Nonfinite crossed statistic')
    a,b,c,d=values
    support=.5*((c-a)+(d-b));denominator=.5*((b-a)+(d-c))
    return dict(available=True,total=d-a,support_path_average=support,
        denominator_path_average=denominator,interaction=d-c-b+a,
        closure=support+denominator-(d-a),
        support_first=c-a,denominator_after_support=d-c,
        denominator_first=b-a,support_after_denominator=d-b)


def crossed_fits(wave,flux,centre,references,supports,*,clip):
    """Cross two declared masks with two references without moving pixels.

    All masks refer to the same delivered observation grid. ``clip=False``
    keeps the coefficient support exactly fixed. ``clip=True`` freezes the
    initial selection only; subsequent clipping belongs to that corner.
    """
    if set(references)!={'before','after'} or set(supports)!={'before','after'}:
        raise ValueError('Require both explicit crossed identities')
    result={}
    for support in ('before','after'):
        mask=np.asarray(supports[support])
        if mask.dtype!=bool:raise ValueError('Crossed support must be a Boolean mask')
        for denominator in ('before','after'):
            trace,continuum=trace_fit(wave,flux,centre,references[denominator],mask,clip=clip)
            if not clip and trace['applied']:
                np.testing.assert_array_equal(trace['fit_indices'],np.flatnonzero(mask))
            result[support+'_'+denominator]=(trace,continuum)
    return result
