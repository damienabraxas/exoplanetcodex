"""Read native absdat wavelength knots exactly as composition_mod does.

KVADL controls interpolation order, not wavelength frame. Native tables and HOP
queries use Angstrom; the experimental bridge converts these coordinates to
vacuum nm exactly once. No table values or opacity coefficients are modified.
"""
import re
import numpy as np


def read_background_knots(path):
    lines=path.read_text().splitlines();tables=[]
    for index,line in enumerate(lines):
        match=re.search(r'\bNLATB\s*=\s*([+-]?\d+)',line,re.I)
        if match is None:
            if re.search(r'\bNLATB\b',line,re.I):raise ValueError('Malformed native table knot count')
            continue
        count=int(match[1])
        if not 0<count<=10000:raise ValueError('Unsupported native table knot count')
        values=[]
        for row in lines[index+1:]:
            for token in row.split():
                try:value=float(token.lower().replace('d','e'))
                except ValueError as error:raise ValueError('Malformed native wavelength record')from error
                values.append(value)
                if len(values)==count:break
            if len(values)==count:break
        wave=np.asarray(values,float)
        if len(wave)!=count or not np.isfinite(wave).all()or np.any(wave<0)or np.any(np.diff(wave)<=0):raise ValueError('Invalid native table wavelengths')
        # Match composition_mod::get_table_sp's explicit max(1 AA, lambda).
        wave=np.maximum(wave,1.)*.1
        tables.append(dict(label=lines[index-1].strip(),header_line_1based=index+1,wave_vac_nm=wave))
    if not tables:raise ValueError('No native background wavelength tables')
    return tables
