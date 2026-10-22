# apply_heatflux_intervals.py
# Run inside ANSYS Mechanical: Automation > Scripting, from the Shell (the prompt
# under the editor), which runs without the editor's debugger:
#   exec(open(r"C:\Users\Kevin Kyle\Documents\Hypersonics\FEA\ANSYS_SCRIPTS\apply_heatflux_intervals.py").read())
#
# Applies raw CFD wall heat flux, one case per time interval, from
# export_cfd_heatflux_intervals.m outputs:
#   1. External Data from <run>_heatflux_faces.csv (X/Y/Z [m], q_<case> [W/m^2])
#   2. Imported Heat Flux on the SKIN named selection. Data rows hold each case
#      for its whole interval and switch to the next case JUMP_S after the
#      boundary (a near-step, since Mechanical ramps between rows)
#   3. (optional) Radiation on SKIN, so the skin can shed heat
#   4. (optional) Analysis Settings steps that end at every interval boundary
#
# Same structure as apply_trajectory_loads.py (the h/T_aw convection version);
# only the load type and the row layout differ. The script prints the data
# table's columns before writing and stops if they are not the expected ones:
# paste that printout back.
#
# Safe to re-run: it first deletes the Imported Load group it created earlier
# (matched by the faces file path) and its Radiation. It STOPS if the analysis
# already has another Imported Load group (for example the convection loads
# from apply_trajectory_loads.py), because both together would double the heating.
import os

# ---------------- user settings ----------------
EXPORTS     = r"C:\Users\Kevin Kyle\Documents\Hypersonics\FEA\Exports"
RUN_NAME    = "WAV_RID_traj9"
SKIN_NS     = "Skin"          # named selection of the heated outer faces
ANALYSIS    = None            # analysis name, or None = first Transient Thermal

JUMP_S      = 0.001           # [s] ramp from one case to the next, after each boundary
ALLOW_OTHER_IMPORTED = False  # True = keep other Imported Load groups (risk of double heating)

ADD_RADIATION = True
EMISSIVITY    = 0.85          # CHAMPS used 0.85; set your surface's real value
T_RADIATE_K   = 250.0         # radiation sink temperature [K]

SETUP_STEPS   = True          # one step per interval
T_INITIAL_S   = 0.05          # initial time step after every boundary [s]
T_MIN_S       = 0.001
T_MAX_S       = 2.0
# ------------------------------------------------

FACES = os.path.join(EXPORTS, RUN_NAME + "_heatflux_faces.csv")
INTERVALS = os.path.join(EXPORTS, RUN_NAME + "_heatflux_intervals.csv")
RAD_NAME = "Radiation (heat flux script)"
VT = MechanicalEnums.ExternalData.VariableType

LOG = os.path.join(EXPORTS, "apply_heatflux_intervals.log")   # progress log (the Shell hides print output)
open(LOG, "w").close()


def say(msg):
    print(msg)
    with open(LOG, "a") as fh:
        fh.write(msg + "\n")


def fail(msg):
    say("STOPPED: " + msg)
    raise Exception(msg)


# ---------------- read the CSVs ----------------
for f in (FACES, INTERVALS):
    if not os.path.isfile(f):
        fail("missing " + f + " (run export_cfd_heatflux_intervals.m first)")

with open(FACES) as fh:
    header = fh.readline().strip().split(",")
qCols = [(i, c) for i, c in enumerate(header) if c.startswith("q_")]
if header[1:4] != ["X", "Y", "Z"] or not qCols:
    fail("unexpected header in %s: %s" % (FACES, header))

# plain split, not csv.DictReader: Mechanical's IronPython reports the reader's
# normal end-of-file StopIteration as an error (the file has no quoted fields)
with open(INTERVALS) as fh:
    lines = [ln.strip() for ln in fh.read().splitlines() if ln.strip()]
cols = lines[0].split(",")
ivals = [dict(zip(cols, ln.split(","))) for ln in lines[1:]]
for iv in ivals:
    if iv["qColumn"] not in header:
        fail("intervals file references column %s, which %s does not have" % (iv["qColumn"], FACES))
    iv["t0"] = float(iv["t_start_s"])
    iv["t1"] = float(iv["t_end_s"])
for k in range(1, len(ivals)):
    if abs(ivals[k]["t0"] - ivals[k - 1]["t1"]) > 1e-9:
        fail("intervals do not abut: %s ends at %s, %s starts at %s"
             % (ivals[k - 1]["Case"], ivals[k - 1]["t1"], ivals[k]["Case"], ivals[k]["t0"]))

# rows: hold each case over its interval, switch JUMP_S after the boundary
#   t0_1: q1 | t1_1: q1 | t1_1+J: q2 | t1_2: q2 | t1_2+J: q3 | ... | t1_N: qN
rows = [(ivals[0]["t0"], ivals[0]["qColumn"])]
for k in range(len(ivals)):
    rows.append((ivals[k]["t1"], ivals[k]["qColumn"]))
    if k + 1 < len(ivals):
        rows.append((ivals[k]["t1"] + JUMP_S, ivals[k + 1]["qColumn"]))
for k in range(1, len(rows)):
    if rows[k][0] <= rows[k - 1][0]:
        fail("row times must increase; is an interval shorter than JUMP_S? (%s)" % (rows,))
say("Faces file: %d cases; %d intervals, t = %s to %s s; %d table rows"
    % (len(qCols), len(ivals), ivals[0]["t0"], ivals[-1]["t1"], len(rows)))
for iv in ivals:
    say("  %-12s %8.2f - %8.2f s  -> %s" % (iv["Case"], iv["t0"], iv["t1"], iv["qColumn"]))

# ---------------- analysis and named selection ----------------
analysis = None
for a in Model.Analyses:
    if ANALYSIS is not None and a.Name == ANALYSIS:
        analysis = a
    elif ANALYSIS is None and a.AnalysisType == AnalysisType.Transient and a.PhysicsType == PhysicsType.Thermal:
        analysis = analysis or a
if analysis is None:
    fail("no transient thermal analysis found")
skin = [ns for ns in Model.NamedSelections.Children if ns.Name == SKIN_NS] if Model.NamedSelections else []
if not skin:
    fail("named selection '%s' not found" % SKIN_NS)
skin = skin[0]
say("Analysis: %s, scope: %s" % (analysis.Name, skin.Name))

# ---------------- remove what an earlier run of this script left ----------------
# Keep a live object selected in the tree after every change: Ansys's own
# context-menu code assumes one is (see apply_trajectory_loads.py).
tree = ExtAPI.DataModel.Tree


def keep(obj=None):
    tree.Activate(obj if obj is not None else analysis)


keep()
kids = list(analysis.Children)
others = []
for g in [c for c in kids if c.DataModelObjectCategory == DataModelObjectCategory.ImportedLoadGroup]:
    fs = g.GetExternalDataFiles()
    paths = [str(fs[i].FilePath).lower() for i in range(fs.Count)]
    if FACES.lower() in paths:
        say("Deleting earlier Imported Load group '%s'" % g.Name)
        g.Delete()
        keep()
    else:
        others.append("%s (%s)" % (g.Name, ", ".join(paths)))
if others and not ALLOW_OTHER_IMPORTED:
    fail("the analysis already has other Imported Load group(s): %s. Delete or suppress them "
         "(applying convection and heat flux together doubles the heating), or set "
         "ALLOW_OTHER_IMPORTED = True" % "; ".join(others))
for b in [c for c in kids if c.DataModelObjectCategory == DataModelObjectCategory.Radiation]:
    if b.Name == RAD_NAME:
        say("Deleting earlier '%s'" % b.Name)
        b.Delete()
        keep()

# ---------------- external data ----------------
if not hasattr(VT, "HeatFlux"):
    fail("this Mechanical version has no ExternalData VariableType.HeatFlux; available: %s"
         % ", ".join([n for n in dir(VT) if not n.startswith("_")]))
files = Ansys.Mechanical.ExternalData.ExternalDataFileCollection()
files.SaveFilesWithProject = False
edf = Ansys.Mechanical.ExternalData.ExternalDataFile()
edf.Identifier = "File1"
edf.FilePath = FACES
edf.IsMainFile = True
edf.ImportSettings = Ansys.Mechanical.ExternalData.ImportSettingsFactory.GetSettingsForFormat(
    MechanicalEnums.ExternalData.ImportFormat.Delimited)
st = edf.ImportSettings
st.SkipRows = 1
st.SkipFooter = 0
st.Delimiter = ","
st.AverageCornerNodesToMidsideNodes = False
st.UseColumn(1, VT.XCoordinate, "m", "X")
st.UseColumn(2, VT.YCoordinate, "m", "Y")
st.UseColumn(3, VT.ZCoordinate, "m", "Z")
for i, c in qCols:
    st.UseColumn(i, VT.HeatFlux, "W m^-2", c)
files.Add(edf)

group = analysis.AddImportedLoadExternalData()
keep(group)
group.ImportExternalDataFiles(files)
keep(group)
say("External Data: %s (columns X, Y, Z [m], %d heat flux [W/m^2])" % (FACES, len(qCols)))

hf = group.AddImportedHeatFlux()
keep(hf)
hf.Location = skin
keep(hf)

# ---------------- data table ----------------
# GetTableByName returns an Ansys.ACT.Automation.Mechanical.WorksheetTable:
#   .Columns (column names), .Count, tbl[i] -> WorksheetRow, row[colIndex] get/set,
#   .Add(WorksheetRow()). No try/except on purpose (see apply_trajectory_loads.py).
tbl = hf.GetTableByName("")
names = [str(n) for n in tbl.Columns]
row0 = tbl[0]
say("Table: %d column(s), %d row(s)" % (len(names), tbl.Count))
for k in range(len(names)):
    v = row0[k]
    say("  col %d %-28s row0 = %-30s (%s)" % (k, "'" + names[k] + "'", str(v), type(v).__name__))
low = [n.lower() for n in names]


def find(*must, **kw):
    bad = kw.get("bad", ())
    hits = [k for k in range(len(low))
            if all([m in low[k] for m in must]) and not any([b in low[k] for b in bad])]
    if len(hits) == 1:
        return hits[0]
    return None


iQ    = find("heat flux", bad=("scale", "offset"))
iTime = find("analysis time")
iQS   = find("scale", bad=("offset",))      # optional: left at its default of 1
if None in (iQ, iTime):
    fail("could not identify the heat flux / analysis-time columns from the list above; "
         "paste the 'col' lines back")

items = [str(x) for x in tbl.DropDownItems(iQ)]
say("  dropdown col %d: %s" % (iQ, items))
qItems = {}
for it in items:
    qItems[it.split(":")[-1].split("@")[0].strip()] = it
missing = [c for _, c in qCols if c not in qItems]
if missing:
    fail("dropdown has no entry for %s; see the dropdown line in %s" % (missing, LOG))
timeIsQty = "Quantity" in type(row0[iTime]).__name__


while tbl.Count < len(rows):
    tbl.Add(Ansys.ACT.Automation.Mechanical.WorksheetRow())
for k in range(len(rows)):
    t, col = rows[k]
    row = tbl[k]
    row[iQ] = qItems[col]
    row[iTime] = Quantity(t, "s") if timeIsQty else t

# read back every row to prove the writes landed
for k in range(len(rows)):
    got = tbl[k]
    say("row %d: time %s | heat flux %s" % (k, got[iTime], got[iQ]))
if tbl.Count != len(rows):
    fail("table has %d rows, expected %d" % (tbl.Count, len(rows)))

say("Imported Heat Flux: %d rows set up on %s (not imported yet: set the mapping options first)" % (len(rows), skin.Name))

# ---------------- radiation ----------------
if ADD_RADIATION:
    rad = analysis.AddRadiation()
    keep(rad)
    rad.Name = RAD_NAME
    rad.Location = skin
    rad.Emissivity.Output.DiscreteValues = [Quantity(EMISSIVITY, "")]
    rad.AmbientTemperature.Output.DiscreteValues = [Quantity(T_RADIATE_K, "K")]
    say("Radiation: emissivity %.2f to %.0f K on %s" % (EMISSIVITY, T_RADIATE_K, skin.Name))

# ---------------- analysis settings: a step per interval ----------------
if SETUP_STEPS:
    ends = [iv["t1"] for iv in ivals]
    s = analysis.AnalysisSettings
    keep(s)
    s.NumberOfSteps = len(ends)
    # end times must increase step to step; new steps default to 1, 2, 3 ... s, so set
    # them last-first (otherwise Mechanical silently rejects an early large end time)
    for n in range(len(ends), 0, -1):
        s.SetStepEndTime(n, Quantity(ends[n - 1], "s"))
    for n in range(1, len(ends) + 1):
        s.SetAutomaticTimeStepping(n, AutomaticTimeStepping.On)
        s.SetDefineBy(n, TimeStepDefineByType.Time)
        s.SetInitialTimeStep(n, Quantity(T_INITIAL_S, "s"))
        s.SetMinimumTimeStep(n, Quantity(T_MIN_S, "s"))
        s.SetMaximumTimeStep(n, Quantity(T_MAX_S, "s"))
    say("Analysis Settings: %d steps ending at %s s" % (len(ends), ", ".join(["%.1f" % e for e in ends])))

say("Done. By hand on Imported Heat Flux: Tabular Loading = Ramped, mapping = Nearest, "
    "Outside Option = Ignore; then right-click > Import Load and check contours at a time inside each interval.")
