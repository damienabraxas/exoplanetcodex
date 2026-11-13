# Pinned QLG profile diagnostic

The driver calls the exact upstream QLG-bot and fp64 reference routines from
OpenMULTI3D application commit32085d26c46e9919a986397b5514a9728abe2074.
Original source and reference binaries remain unchanged. This is a kernel probe,
not a transfer/SE calculation or a model/default choice. The native active non-H
line path directly calls QLG-bot; standalone HOP line-list options do not select it.

GPL-3.0-only, like the upstream application. License text is retained in
[../hop_trace_provider/COPYING](../hop_trace_provider/COPYING); upstream attribution
is in the pinned application's AUTHORS and QLG module. The module is compiled
from that verified external checkout, not copied into this repository.
