# FEA

**Start in [`THERMAL_PIPELINE/`](THERMAL_PIPELINE/README.md).** It has everything for the
trajectory heating runs, step by step.

| Folder | What |
|---|---|
| `THERMAL_PIPELINE/` | Current workflow: trajectory → CFD cases → ANSYS loads → transient thermal |
| `WAV_RID_9_FINS/` | Current vehicle: SolidWorks part, Discovery file, trajectory |
| `ANSYS_MATS/` | Granta Selector material exports for Engineering Data (see below) |
| `Old/` | Earlier work, kept for reference. `Old/radiative_pipeline/` is the previous version of the pipeline |

## Materials (`ANSYS_MATS/`)

Import **`WAV_RID_9_materials.xml`** in Engineering Data: File → Import Engineering Data.
It has all four, ready to assign (density, conductivity, specific heat vs. temperature):

| Bodies | Material |
|---|---|
| INCONEL_* | Nickel-chromium alloy, INCONEL 718, solution treated & aged |
| Titanium spar | Ti - 6 Al - 4 V (Grade 5), annealed |
| PAYLOAD_AL, AVIONICS_AL | 6061, T6 temper |
| TUNGSTEN | Tungsten, solid |

Still needed: the CCC (carbon-carbon) and INSULATION bodies (Granta Selector, same way).
`granta_records/` holds the original single-material exports.
