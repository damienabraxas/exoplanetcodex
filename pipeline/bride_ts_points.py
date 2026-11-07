"""Input-only adapter for isolated TS opacity diagnostics, never atmosphere adoption."""
import numpy as np


def point_input(temperature_K, electron_pressure_Pa):
    """Serialize positive finite point arrays for default-REAL TS in K/dyn cm^-2."""
    t = np.asarray(temperature_K, dtype=float)
    pe = np.asarray(electron_pressure_Pa, dtype=float)
    if t.ndim != 1 or len(t) < 2 or pe.shape != t.shape or not np.isfinite(t).all() or not np.isfinite(pe).all() or np.any(t <= 0) or np.any(pe <= 0):
        raise ValueError('Require matched finite positive T/Pe point arrays')
    # Explicitly require representability at the pinned producer's REAL precision.
    with np.errstate(over='ignore', under='ignore'):
        producer = np.column_stack((t, pe/.1)).astype(np.float32)
    if not np.isfinite(producer).all() or np.any(producer <= 0):
        raise ValueError('Input not representable in native TS REAL')
    return str(len(t))+'\n'+''.join(f'{a:.17e} {b/.1:.17e}\n' for a, b in zip(t, pe))


def point_driver_source(source):
    """Patch only input loading in a copy of pinned BABSMA, before opacity calls.

    Original MARCS file initializes composition and native driver infrastructure.
    Its geometry/tau are placeholders for opacity evaluation, not a newly solved
    atmosphere. Every point state is required, with no truncation/fallback.
    """
    declaration = '      PROGRAM BABSMA\n'
    anchor = '      NEWT=2\n'
    if source.count(declaration) != 1 or source.count(anchor) != 1 or 'BRIDE_' in source or 'point_state.txt' in source:
        raise ValueError('Unexpected pinned BABSMA input insertion sites')
    source = source.replace(declaration, declaration+'      INTEGER BRIDE_N,BRIDE_I,BRIDE_IOS\n')
    block = """C     Isolated diagnostic T/Pe override; not an equilibrium atmosphere.
      OPEN(97,FILE='point_state.txt',STATUS='OLD',IOSTAT=BRIDE_IOS)
      IF(BRIDE_IOS.NE.0) STOP 71
      BRIDE_N=0
      READ(97,*,IOSTAT=BRIDE_IOS) BRIDE_N
      IF(BRIDE_IOS.NE.0.OR.BRIDE_N.NE.NTAU) STOP 72
      DO BRIDE_I=1,NTAU
        READ(97,*,IOSTAT=BRIDE_IOS) T(BRIDE_I),PE(BRIDE_I)
        IF(BRIDE_IOS.NE.0) STOP 73
        IF(T(BRIDE_I).LE.0.OR.PE(BRIDE_I).LE.0) STOP 74
        IF(T(BRIDE_I).NE.T(BRIDE_I)) STOP 75
        IF(PE(BRIDE_I).NE.PE(BRIDE_I)) STOP 75
      ENDDO
      READ(97,*,IOSTAT=BRIDE_IOS) BRIDE_N
      IF(BRIDE_IOS.GE.0) STOP 76
      CLOSE(97)
"""
    return source.replace(anchor, block+anchor)
