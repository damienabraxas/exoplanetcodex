"""Require the inherited rate-preserving operator; never silently use backend default."""
def validate_operator_contract(plan,source_report):
    expected=source_report['iteration_plugin']['receipt']
    item=plan.get('iteration_plugin')
    if item is None or item.get('resource')!='experiments.bride_full_cap_plugin' or item.get('build_receipt')!=expected:
        raise ValueError('Coupled solve requires the explicitly pinned inherited rate-preserving operator')
    return item
