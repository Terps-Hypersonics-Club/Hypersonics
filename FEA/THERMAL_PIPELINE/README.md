# THERMAL_PIPELINE

Trajectory → CFD → ANSYS transient thermal. Do the steps in order.
Needs MATLAB with the Statistics and Machine Learning Toolbox, CHAMPS, and ANSYS 2026 R1.

## 1. Plan the CFD runs (MATLAB)
Run `step1_plan_cfd_cases.m`. It prints a table of CHAMPS runs:
```
VTK_file     mach    aoa    Pref[Pa]   Tref[K]
case01.vtk   7.958   0.182  39715.7    241.12
...
```
Fewer or more runs: change `nCases` at the top (now 10).

## 2. Run CHAMPS (one run per row)
Copy `inputs/CHAMPS_input_TEMPLATE.sdf` and set:
- `mach`, `aoa`, `Pref`, `Tref` from the row
- `filename` to the WAV_RID_9_FINS CHAMPS geometry (same file for every run)
- wall: **isothermal 315 K** (the template is radiative, change `radbc`)

Rename each surface output to its `VTK_file` name and put it in `cfd/`.

## 3. Build the loads (MATLAB)
1. Export the vehicle's outer surface as STL from the **same file you import into ANSYS**.
   Save it as `geometry/WAV_RID_9_FINS.stl`. If it is in mm, set `stlUnits = 'mm'` in step 2.
2. Run `step2_build_ansys_loads.m`. It writes two CSVs to `outputs/`.

## 4. ANSYS Mechanical
Geometry (once):
1. Workbench: Transient Thermal → right-click Geometry → Import →
   `../WAV_RID_9_FINS/WAV_RID_9_FINS_FINAL.SLDPRT`. In Geometry Properties, untick Surface Bodies.
2. Right-click Geometry → Edit in Discovery. There should be 23 solid bodies. Delete any loose surfaces.
3. Cut in half: select all bodies → Split with the vertical plane through the centerline
   (the plane the fins mirror across) → delete one half.
4. Close Discovery. Leave bodies separate: Mechanical bonds touching bodies with contacts.
   The cut face needs no load (no heat flow across it = symmetry).

Set up first:
- Materials assigned (see `../README.md`), mesh done
- Transient Thermal: Initial Temperature 315 K
- Named selection **`Skin`** = outer heated faces (not the cut face)
- Analysis Settings > Output Controls: store results at ~100 points, not all (a full run was 47 GB)

Then:
1. In `step3_apply_loads_mechanical.py`, set `EXPORTS` to your `outputs` folder.
2. Automation > Scripting, in the Shell (bottom prompt):
   `exec(open(r"<path>\step3_apply_loads_mechanical.py").read())`
3. On the new Imported Convection: Tabular Loading = Ramped, Mapping = Nearest,
   Outside Option = Ignore. Right-click > **Import Load**. Check the contours look like heating.
4. **Solve.**

If step 3 stops, open `outputs/apply_trajectory_loads.log` and send it to Kevin.

## 5. Export temperatures (optional)
Set `OUT` and `TIMES` in `step4_export_temperatures.py`, run it like step 3.

---
Method in one line: each CFD run gives h = q_w / (T_aw − 315 K), scaled to the trajectory
between runs. ANSYS applies q = h (T_aw − T_structure), plus radiation (ε 0.85).
Conservative: h is the cold-wall value.
