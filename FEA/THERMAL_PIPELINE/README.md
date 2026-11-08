# THERMAL_PIPELINE

CFD heating → ANSYS transient thermal along the trajectory, for WAV_RID_9B.
This folder is set up exactly as Kevin's first full run (Nov 2026). Do the steps in order.
Needs MATLAB with the Statistics and Machine Learning Toolbox, and ANSYS 2026 R1.

## What is already here
| Path | What |
|---|---|
| `cfd/` | The 5 CHAMPS runs used (radiative wall, WAV_RID_9B geometry) |
| `geometry/WAV_RID_9B.stl` | Half-model outer surface, mm (lines the CFD up with ANSYS) |
| `inputs/Trajectory_WAV_RID_COMP_aero.csv` | Trajectory |
| `../WAV_RID_9B/WAV_RID_9B_ANSUS_INTERFACED.dsco` | Vehicle geometry: overlaps already fixed, bodies named |

## 1. Assign CFD runs to the trajectory (MATLAB)
Run `step1_plan_cfd_cases.m`. The `fixedCases` table at the top lists the 5 runs in `cfd/`
(file, Mach, AoA, altitude). Each trajectory time gets the nearest run.
New CFD runs: add a row per file. Empty `fixedCases = {}` plans brand-new cases instead.

### Using your own CFD runs
1. Put the surface VTKs in `cfd/` (remove the ones you are not using).
2. `step1_plan_cfd_cases.m`, `fixedCases`: one row per file, with its Mach, AoA [deg] and
   altitude [km]. The scripts assume each run's freestream is the standard atmosphere at that
   altitude; if you ran other Pref/Tref values, use the altitude that matches them.
3. `step2_build_ansys_loads.m`, `wallBC`: `'radiative'` or `'isothermal'` (then `T_wall_cfd` =
   your wall temperature). All runs must use the same wall BC; step 2 stops if one does not match.
4. All runs must be on the same CHAMPS geometry file, and that geometry must match the model in
   ANSYS. If you changed the vehicle, export a new STL of it into `geometry/` and set `stlFile`.
   Step 2 stops if the runs' surface meshes differ, and warns if CFD and STL are more than a few mm apart.

## 2. Build the loads (MATLAB)
Run `step2_build_ansys_loads.m`. Settings at the top are already right for this data:
`stlFile` = `geometry/WAV_RID_9B.stl`, `stlUnits = 'mm'`, `wallBC = 'radiative'`.
It writes `outputs/WAV_RID_9_FINS_faces.csv` and `..._time_table.csv`.
A warning about a few-mm gap between CFD and STL is normal.

## 3. ANSYS model (once)
1. Workbench: Transient Thermal → Geometry → Import `../WAV_RID_9B/WAV_RID_9B_ANSUS_INTERFACED.dsco`.
   Geometry Properties: Surface Bodies off.
2. Edit in Discovery. If the bodies span both sides of y = 0, cut in half: select all bodies →
   Split with the centerline plane (the plane the fins mirror across) → delete the +y half.
   Then look for leftover small fin fragments and delete them.
3. Materials: import `../ANSYS_MATS/WAV_RID_9_materials.xml`, `CC.xml` and `Min_K, Conductivity.xml`
   (see `../README.md` for which body gets which). Min-K: add Specific Heat = 1089 J/kg·K.
4. Mesh. Thin fins and leading edges may need Patch Independent + a smaller size.
5. Named selection **`Skin`**: every outer face the flow touches, **top and bottom** (not the cut face).
6. Optional, for full-vehicle pictures: Named selection of the cut faces (Worksheet: Face,
   Location Y = 0) → Model → Insert Symmetry → Symmetry Region (Symmetric, normal Y).
   Results then show both halves with Display → Visual Expansion.

## 4. Apply loads and solve (ANSYS)
1. Delete any old Imported Load and "Radiation (trajectory script)". Initial Temperature 315 K.
2. In `step3_apply_loads_mechanical.py`, set `EXPORTS` to your `outputs` folder.
3. Automation > Scripting, in the **Shell** (bottom prompt, not the Run button):
   `exec(open(r"<path>\step3_apply_loads_mechanical.py").read())`
4. On the new Imported Convection: Tabular Loading = **Ramped**,
   Weighting = **Distance Based Average**, Number of Points = **1**. Right-click → **Import Load**.
   Every outer face must be colored; gray means it is missing from `Skin`.
5. Analysis Settings: step 1 max time step 2 s; steps 2-5 initial 0.5 s, max 10 s.
   Output Controls: Store Results At = **All Time Points**, Calculate Thermal Flux = Yes.
   Expect a few hours of solving and ~40-60 GB of results.
6. Solve. If step 3 stopped, send `outputs/apply_trajectory_loads.log` to Kevin.

## 5. Export results (ANSYS)
File → Options → Export → Include Node Location = Yes. Then in the Shell:
`exec(open(r"<path>\step4_export_temperatures.py").read())`
Writes temperature (`T_*.txt`) and total heat flux (`q_*.txt`) every 10 s to `outputs/snapshots/`.
The 0 s file shows the last time point (Mechanical quirk); ignore it.

---
Method: each CFD run gives h = q_w / (T_aw − T_w,CFD), scaled to the trajectory's density and
speed between runs. ANSYS applies q = h (T_aw − T_structure) plus surface radiation (ε 0.85).
**Preliminary:** the CFD wall was radiative (hot), so body heating is likely low early in flight,
and 514-766 s (Mach 3.2-5, AoA 2-8°) has no close CFD run. Isothermal 315 K CFD runs fix the
first issue: set `wallBC = 'isothermal'` in step 2.
