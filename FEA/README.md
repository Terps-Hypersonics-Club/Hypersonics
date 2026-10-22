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

Import in Engineering Data: Engineering Data Sources → "Click here to add a new library" → pick the XML.

| Body | File | Note |
|---|---|---|
| Inconel parts | `Inconel_718_solution_treated.xml` | |
| Titanium spar | `Ti6Al4V_Grade5_Annealed.xml` | No specific heat: copy it in from `Ti6Al4V_GeneralProperties.xml` |
| Payload, avionics | `Al6061_T6.xml` | |
| Nose | `Tunsten_solid_Binkele (Conductivity).xml` | No specific heat: copy it in from `Tungsten_solid_NO_conductivity.xml` |
| ZrB2 | `ZrB2, general properties.xml` | |

Still needed: the CCC (carbon-carbon) and INSULATION bodies.

ANSYS project files (`*.wbpj` and their `_files` folders) are not on GitHub because of size.
Share them through Drive.
