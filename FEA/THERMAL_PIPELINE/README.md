# THERMAL_PIPELINE: trajectory heating for WAV_RID_9 (fins body)

Everything needed to go from a trajectory to a transient thermal solve in ANSYS.
Run the steps in order. Nothing outside this folder is needed.

```
trajectory ──step1──▶ CFD case list ──CHAMPS (you)──▶ cfd/caseNN.vtk ──step2──▶ loads CSVs ──step3──▶ ANSYS solve ──step4──▶ temperature exports
```

## Folder

| Path | What |
|---|---|
| `step1_plan_cfd_cases.m` | MATLAB. Splits the trajectory into CFD cases, prints what to run in CHAMPS |
| `step2_build_ansys_loads.m` | MATLAB. Turns the CHAMPS outputs into ANSYS convection loads |
| `step3_apply_loads_mechanical.py` | ANSYS Mechanical script. Imports the loads, radiation and time steps |
| `step4_export_temperatures.py` | ANSYS Mechanical script. Exports temperatures at chosen times after the solve |
| `inputs/Trajectory_WAV_RID_COMP_aero.csv` | Trajectory (time, altitude, Mach, AoA, ...) |
| `inputs/CHAMPS_input_TEMPLATE.sdf` | CHAMPS input to copy per case (see step A) |
| `reference/` | The Mach 8 sea-level runs at 315 K and 1500 K wall (background, not used by the scripts) |
| `lib/` | Helper functions the MATLAB steps call. Do not edit |
| `cfd/` | Put the CHAMPS surface outputs here as `case01.vtk`, `case02.vtk`, ... |
| `geometry/` | Put the body's outer-surface STL here (step 0) |
| `outputs/` | Everything the scripts write |

## Step 0: geometry STL (once per geometry)

Step 2 lines the CFD points up with the ANSYS geometry using an STL of the outer
surface. Export it from SolidWorks or Discovery from the **same file you import into
ANSYS**, in metres, as `geometry/WAV_RID_9_FINS.stl`. Full body is fine even if
ANSYS uses a half model.

## Step 1: plan the CFD cases (MATLAB)

Run `step1_plan_cfd_cases.m`. Settings at the top:
- `trajFile`: the trajectory CSV
- `nCases`: how many CHAMPS runs you can afford (10 now)

It writes `outputs/Trajectory_..._cfd_cases.csv` and prints a table like:

```
VTK_file        mach      aoa   Pref[Pa]   Tref[K]
case01.vtk     7.958    0.182    39715.7    241.12
...
```

## Step A: run CHAMPS (one run per row)

For each row, copy `inputs/CHAMPS_input_TEMPLATE.sdf` and set
`mach`, `aoa`, `Pref`, `Tref` from that row. **Wall BC: isothermal at 315 K**,
not radiative. The template is a radiative run (`radiative = true` in `radbc`);
switch it the same way as the Mach 8 sea-level 315 K runs in `reference/` were set up
(ask whoever ran those for the exact `radbc` lines).

Use the same geometry file in CHAMPS for every case. Save each surface output as
`cfd/case01.vtk`, `cfd/case02.vtk`, ... matching the `VTK_file` column.

## Step 2: build the ANSYS loads (MATLAB)

Run `step2_build_ansys_loads.m`. Check `runName`, `stlFile`, `T_wall_cfd` (315) at the top.
It stops with a message if a VTK is missing, is not a 315 K isothermal run, or is on a
different CHAMPS geometry. Writes:
- `outputs/WAV_RID_9_FINS_faces.csv`: h and T_aw for every CFD point and case
- `outputs/WAV_RID_9_FINS_time_table.csv`: which case applies when, and its scale factors

## Step 3: apply the loads (ANSYS Mechanical)

Before running, the model needs: materials assigned, mesh done, a Transient Thermal
analysis with an initial temperature (315 K), and a named selection **`Skin`** on the
heated outer faces (not the symmetry cut face).

1. Open `step3_apply_loads_mechanical.py`, set `EXPORTS` to your `outputs` folder.
2. Mechanical: Automation > Scripting, in the Shell (bottom prompt) run
   `exec(open(r"<path>\step3_apply_loads_mechanical.py").read())`
3. On the new Imported Convection: Tabular Loading = Ramped, Mapping = Nearest,
   Outside Option = Ignore. Right-click > Import Load.
4. Check the imported contours at a few times (one per case window).
5. Solve.

If it stops, read `outputs/apply_trajectory_loads.log` (it prints what it found).

## Step 4: export temperatures (ANSYS Mechanical, after the solve)

Set `OUT` and `TIMES` in `step4_export_temperatures.py`, run it the same way.
Needs File > Options > Export > Text File Export > Include Node Location = Yes.

## What the method assumes

- **Heating = convection with h and T_aw**, so it drops as the structure heats.
  `h = q_w / (T_aw - 315 K)` from each isothermal run, T_aw from the local CFD edge state.
- **h is the cold-wall value.** Exact at 315 K, too high by up to ~1.5x on the body once the
  wall nears 1500 K (measured with the two runs in `reference/`). So results are conservative.
- **Between cases** each case is scaled to the trajectory's density and speed:
  h x (rho V)^0.8 ratio, T_aw x recovery-temperature ratio. Turbulent form everywhere
  (nose and leading edges would scale more weakly). Mach and AoA differences between a case
  and the trajectory are not corrected; `nCases` controls how large they get.
- Ideal gas (gamma 1.4), CHAMPS Euler + Zoby boundary layer, surface radiation
  epsilon 0.85 to 250 K, no ablation.
- In the Mach 8 runs about 8% of the surface aft of x = 1.04 m has q_w = 0 in the CFD,
  so it gets no heating. Worth a look in ParaView.
