# RYA-419 — Sirius 500 GB data-drive runbook (RUN ON SIRIUS)

> **EXECUTED 2026-06-23** (over SSH, once RYA-80 key-auth + NOPASSWD sudo landed).
> Actual hardware (supersedes the stale RYA-113 "256 + 500" assumption):
> - `sdb` 476.9 GiB (~512 GB) = **OS drive** — `sdb1` 1G vfat `/boot/efi`, `sdb2` ext4 `/` (395G free).
> - `sda` 465.8 GiB (**500 GB**) = **data drive** — was BLANK; formatted GPT + `sda1` ext4 (label `codex-data`)
>   per Ryan's sign-off. **UUID `75f9235a-0d9f-4b77-9eb8-3eb9474a2c88`**, mounted at `/mnt/codex-data`
>   (458G usable). fstab line (UUID, `nofail`) verified via `mount -a`; all 5 smoke checks pass incl. the
>   negative test (unmounted → CRITICAL/exit 1). Skeleton created. The steps below are the as-run procedure.



**This is a Linux/Sirius (HP ProBook 450) ops runbook.** The drive is physically on
Sirius; the commands below use `lsblk`/`/etc/fstab`/`mountpoint`/`systemctl` and `sudo`,
which do not exist / do not apply on the MacBook. They were **not** run from the Mac
session — run them on Sirius (or via SSH once RYA-80 key-auth is set up). The repo half
of this issue (`scripts/check_data_mount.sh`, `scripts/init_codex_data_skeleton.sh`) is
already committed and pulled onto Sirius with the repo.

Layout reality (supersedes the stale RYA-113 description): **256 GB main SSD** (OS + repo
+ backend) · **500 GB data drive** (spectra, VALD, grids, outputs).

## Step 0 — discover (report only, no changes)
```bash
lsblk -o NAME,SIZE,FSTYPE,MOUNTPOINT,UUID,LABEL
sudo blkid
df -h
```
Record the ~500 GB partition's **UUID** and **FSTYPE**. Every step below uses the UUID,
never `/dev/sdX` (device names renumber across reboots).

## Step 1 — filesystem decision (DESTRUCTIVE GUARD)
- Partition **has** a filesystem (ext4/exfat/ntfs) or any data → **use as-is. Do NOT format.**
- Partition is **blank/unformatted** → **STOP. Do not `mkfs`.** Report device + size to Ryan
  with a recommendation (ext4 for a Linux-only data drive) and wait for explicit go.
  Formatting is irreversible — never auto-format.

## Step 2 — mount point + first mount (by UUID)
```bash
sudo mkdir -p /mnt/codex-data
sudo mount UUID=<UUID> /mnt/codex-data
df -h /mnt/codex-data    # MUST show the 500 GB drive, NOT the 256 GB root
```

## Step 3 — persistence (fstab, UUID, nofail)
```bash
sudo cp /etc/fstab /etc/fstab.bak.$(date +%Y%m%d)
echo "UUID=<UUID>  /mnt/codex-data  <FSTYPE>  defaults,nofail,x-systemd.device-timeout=10  0  2" | sudo tee -a /etc/fstab
sudo systemctl daemon-reload
sudo mount -a
```
`nofail` so a missing drive never blocks boot. For **exfat/ntfs** append
`,uid=$(id -u),gid=$(id -g)` to the options (user owns files); for **ext4** use Step 4.

## Step 4 — ownership (ext4 only)
```bash
sudo chown -R "$(id -un):$(id -gn)" /mnt/codex-data
```

## Step 5 — mount-guard sentinel (the silent-fallback kill)
```bash
touch /mnt/codex-data/.codex_mounted     # exists ONLY when the drive is mounted
bash scripts/check_data_mount.sh         # -> OK: /mnt/codex-data mounted and verified.
```
(Repo script already present. Wiring it into the Python preflight is deferred to the
transfer / `CODEX_DATA_ROOT` work — RYA-264/332.)

## Step 6 — folder skeleton
```bash
bash scripts/init_codex_data_skeleton.sh   # guard-gated; idempotent; canonical tokens
```

## Smoke test (run all 5, incl. the negative test)
```bash
# 1
df -h /mnt/codex-data
# 2
bash scripts/check_data_mount.sh                        # OK ...
# 3 — fstab line mounts cleanly without reboot
sudo umount /mnt/codex-data && sudo mount -a && bash scripts/check_data_mount.sh
# 4 — NEGATIVE: the guard MUST catch an unmounted path
sudo umount /mnt/codex-data && bash scripts/check_data_mount.sh; echo "exit=$?"
#    Expect: "CRITICAL: ... not a mount point ... Refusing"  and  exit=1
sudo mount -a
# 5
find /mnt/codex-data -maxdepth 2 -type d | sort
```

**CRITICAL stop conditions:** `df` resolves to the 256 GB root → STOP (fstab/mount wrong).
Negative test (4) does not exit 1 → the guard is broken → STOP before any data is written.
About to `mkfs` a blank drive → STOP, get sign-off.

Report back (Linear): device, UUID, FSTYPE, total/free size; whether a format was needed
(and that none was done without sign-off); verbatim output of all 5 smoke checks.
