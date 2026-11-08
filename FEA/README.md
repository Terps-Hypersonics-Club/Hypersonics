# FEA

**Start in [`THERMAL_PIPELINE/`](THERMAL_PIPELINE/README.md).** It has everything for the
trajectory heating runs, step by step, set up as the first full run.

| Folder | What |
|---|---|
| `THERMAL_PIPELINE/` | Current workflow: CFD runs → ANSYS loads → transient thermal → exports |
| `WAV_RID_9B/` | Current vehicle. Use `WAV_RID_9B_ANSUS_INTERFACED.dsco` (overlaps fixed, bodies named) |
| `ANSYS_MATS/` | Material files for Engineering Data (see below) |
| `WAV_RID_9_FINS/` | Earlier fins version of the vehicle, kept for reference |
| `Old/` | Earlier work. `Old/radiative_pipeline/` is the previous pipeline |

## Materials (`ANSYS_MATS/`)

Engineering Data → File → Import Engineering Data, one file at a time:

| File | Material | Bodies |
|---|---|---|
| `WAV_RID_9_materials.xml` | Inconel 718 (solution treated & aged) | INCONEL_* |
| | 6061-T6 | PAYLOAD_AL, AVIONICS_AL |
| | Tungsten, solid | TUNGSTEN |
| | Ti-6Al-4V (Grade 5) | (no titanium bodies in WAV_RID_9B) |
| `CC.xml` | Carbon/carbon, NASA Langley 9 | CCC_* |
| `Min_K, Conductivity.xml` | Min-K TE1400 board, in air | MINK*, MIN-K* (insulation) |

Min-K: the file has density and conductivity only. Add **Specific Heat = 1089 J/kg·K**
(team EduPack value) in Engineering Data before solving.
`granta_records/` holds the original single-material exports.

ANSYS project files (`*.wbpj` and their `_files` folders) are not on GitHub because of size.
