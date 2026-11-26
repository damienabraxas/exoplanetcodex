#!/bin/sh
set -eu
root=/mnt/codex-data/codex/rya1354
src="$root/runs/mcdhf_layers_v78/n6_E1"
out="$root/runs/mcdhf_e1_v79"
refs="$root/data/reference/li2023_oi/raw"
mkdir "$out"
cd "$out"
exec >driver.stdout 2>driver.stderr
status=FAILED_RETAINED
stage=preflight
receipt() {
    printf '{"status":"%s","stage":"%s","scientific_admission":"HOLD","finished_utc":"%s","files":["driver.sh","launch_sha.txt","driver.stdout","driver.stderr","rtransition.stdin","rtransition.stdout","rtransition.stderr","lower.upper.t","loggf_vs_li2023.csv","comparison.md"]}\n'         "$status" "$stage" "$(date -u +%FT%TZ)" > receipt.tmp
    mv receipt.tmp receipt.json
}
trap receipt EXIT
cp "$0" driver.sh
grep -Fq 'BIOTRANSFORM: Execution complete.' "$src/rbiotransform.stdout"
for f in isodata rtransition.stdin lower.c lower.w lower.m lower.bm lower.bw lower.TB lower.lsj.lbl upper.c upper.w upper.m upper.bm upper.bw upper.TB upper.lsj.lbl; do
    test -s "$src/$f"
    cp "$src/$f" "$f"
    cmp -s "$src/$f" "$f"
done
cp "$root/runs/mcdhf_layers_v78/rya1395_s0_e1_report.py" report.py
cp "$root/evidence/mcdhf_layers_v78_review/compare.py" compare.py
bin="$root/toolchains/grasp/build/bin/rtransition"
sha256sum driver.sh report.py compare.py isodata rtransition.stdin lower.* upper.* "$bin" "$refs/tablea1.tsv" "$refs/tablea2.tsv" > launch_sha.txt
stage=rtransition
timeout 24h "$bin" < rtransition.stdin > rtransition.stdout 2> rtransition.stderr
grep -Fq 'RTRANSITION: Execution complete.' rtransition.stdout
test -s lower.upper.t
stage=transition_report
python3 report.py . "$refs/tablea1.tsv" "$refs/tablea2.tsv" loggf_vs_li2023.csv
test "$(wc -l < loggf_vs_li2023.csv)" -gt 1
stage=comparison
python3 compare.py "$root/runs/mcdhf_layers_v77/n5_E1/loggf_vs_li2023.csv" loggf_vs_li2023.csv > comparison.tmp
mv comparison.tmp comparison.md
status=E1_COMPLETED_REVIEW_REQUIRED
