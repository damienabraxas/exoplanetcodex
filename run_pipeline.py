"""
run_pipeline.py
===============
Single entry point for the full Exoplanet Codex analysis pipeline.

Usage:
    python run_pipeline.py --star 55cancri
    python run_pipeline.py --star 55cancri --validate-only
    python run_pipeline.py --star sol
    python run_pipeline.py --list-stars

Naming convention — pipeline scripts use subject_action pattern:
    spectra_*    → operates on raw/normalized spectra
    lines_*      → operates on spectral line data
    params_*     → derives stellar parameters
    abundances_* → derives elemental abundances
    uncertainty_*→ propagates measurement uncertainties
    ratios_*     → computes/interprets abundance ratios
"""

import argparse
import sys

from pipeline.spectra_acquire   import run as acquire
from pipeline.spectra_normalize import run as normalize
from pipeline.params_stellar    import run as stellar_params
from pipeline.lines_load        import run as load_lines
from pipeline.lines_fit         import run as fit_lines
from pipeline.abundances_derive import run as derive
from pipeline.uncertainty_stack import run as uncertainty
from pipeline.ratios_interpret  import run as interpret

STARS = {
    '55cancri': 'HD 75732',
    'sol':      'solar_asteroid',
    'hd89307':  'HD 89307',
    'gl581':    'Gliese 581',
}


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Exoplanet Codex — stellar abundance pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--star', choices=STARS.keys(),
                        help='Target star identifier')
    parser.add_argument('--validate-only', action='store_true',
                        help='Run only through stellar parameter determination')
    parser.add_argument('--list-stars', action='store_true',
                        help='List available star identifiers and exit')
    args = parser.parse_args()

    if args.list_stars:
        print("\nAvailable stars:")
        for key, hd in STARS.items():
            print(f"  {key:<12} {hd}")
        print()
        sys.exit(0)

    if not args.star:
        parser.print_help()
        sys.exit(1)

    star_id = STARS[args.star]
    print(f"\n{'='*60}")
    print(f"  Exoplanet Codex Pipeline")
    print(f"  Target: {star_id} ({args.star})")
    print(f"{'='*60}\n")

    acquire(star_id)
    normalize(star_id)
    stellar_params(star_id)
    load_lines(star_id)
    fit_lines(star_id)

    if not args.validate_only:
        derive(star_id)
        uncertainty(star_id)
        interpret(star_id)

    print(f"\n✓  Pipeline complete for {star_id}.")


if __name__ == '__main__':
    main()
