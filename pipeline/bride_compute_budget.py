"""Explicit single-case compute allowance; never an accuracy-tolerance override."""
import copy


def extended_controls(original,declaration,case,ng_control):
    keys={'schema','case','ng_control','physical_plan','reference_run',
          'baseline_seconds_per_case','seconds_per_case','reason'}
    if (set(declaration)!=keys or declaration['schema']!='codex.bride_compute_extension/1'
            or declaration['case']!=case or declaration['ng_control']!=ng_control
            or ng_control!='record' or original['ng'] is not True):
        raise ValueError('Outside declared enabled-Ng compute extension')
    seconds=declaration['seconds_per_case']
    if (type(seconds) is not int or type(declaration['baseline_seconds_per_case']) is not int
            or original['seconds_per_case']!=declaration['baseline_seconds_per_case']
            or seconds<=original['seconds_per_case']):
        raise ValueError('Invalid declared compute allowance')
    result=copy.deepcopy(original);result['seconds_per_case']=seconds
    return result
