# step3_apply_loads_mechanical.py   (THERMAL_PIPELINE, see README.md)
# Run inside ANSYS Mechanical: Automation > Scripting. Best from the Shell
# (the prompt under the editor), which runs without the editor's debugger.
# Set EXPORTS below to YOUR THERMAL_PIPELINE\outputs folder, then paste into the Shell
# (with this file's real path on your machine):
#   exec(open(r"C:\...\Hypersonics\FEA\THERMAL_PIPELINE\step3_apply_loads_mechanical.py").read())
# The editor's Run button attaches a debugger that also stops on errors inside
# Ansys's own UI code, which Ansys normally ignores.
#
# Builds the time-dependent CFD heating for the transient thermal analysis from
# step2_build_ansys_loads.m outputs:
#   1. External Data from <run>_faces.csv (X/Y/Z [m], h_<case> [W/m^2K], Taw_<case> [K])
#   2. Imported Convection on the SKIN named selection, one data row per
#      <run>_time_table.csv row:  film coeff = h_<case> * hScale,
#                                 ambient    = Taw_<case> * TawScale
#   3. (optional) Radiation on SKIN
#   4. (optional) Analysis Settings steps that end at every case switch
#
# Written against the 2026 R1 Mechanical API (names checked in the installed
# API reference and, for the undocumented data table, by reflection on
# Ansys.ACT.WB1.dll). The script prints the table's columns before writing and
# stops if they are not the expected ones: paste that printout back.
#
# Safe to re-run: it first deletes the Imported Load group (matched by the faces
# file path) and the Radiation it created on an earlier run.
# Do not press Save in Mechanical's script editor unless you edited the script
# there: the editor writes its open copy back over this file.
import os

# ---------------- user settings ----------------
EXPORTS     = r"C:\Users\Kevin Kyle\Documents\Hypersonics\FEA\THERMAL_PIPELINE\outputs"   # CHANGE to your path
RUN_NAME    = "WAV_RID_9_FINS"    # same as runName in step2_build_ansys_loads.m
SKIN_NS     = "Skin"          # named selection of the heated outer faces
ANALYSIS    = None            # analysis name, or None = first Transient Thermal

ADD_RADIATION = True
EMISSIVITY    = 0.85          # CHAMPS used 0.85; set your surface's real value
T_RADIATE_K   = 250.0         # radiation sink temperature [K] (sky/freestream)

SETUP_STEPS   = True          # one step per case window, ending at each switch
T_INITIAL_S   = 0.05          # initial time step after every switch [s]
T_MIN_S       = 0.001
T_MAX_S       = 2.0
# ------------------------------------------------

FACES = os.path.join(EXPORTS, RUN_NAME + "_faces.csv")
TABLE = os.path.join(EXPORTS, RUN_NAME + "_time_table.csv")
RAD_NAME = "Radiation (trajectory script)"
VT = MechanicalEnums.ExternalData.VariableType


LOG = os.path.join(EXPORTS, "apply_trajectory_loads.log")   # progress log (the Shell hides print output)
open(LOG, "w").close()


def say(msg):
    print(msg)
    with open(LOG, "a") as fh:
        fh.write(msg + "\n")


def fail(msg):
    say("STOPPED: " + msg)
    raise Exception(msg)


# ---------------- read the CSVs ----------------
for f in (FACES, TABLE):
    if not os.path.isfile(f):
        fail("missing " + f + " (run build_ansys_loads.m first)")

with open(FACES) as fh:
    header = fh.readline().strip().split(",")
hCols   = [(i, c) for i, c in enumerate(header) if c.startswith("h_")]
tawCols = [(i, c) for i, c in enumerate(header) if c.startswith("Taw_")]
if header[1:4] != ["X", "Y", "Z"] or not hCols or len(hCols) != len(tawCols):
    fail("unexpected header in %s: %s" % (FACES, header))

# plain split, not csv.DictReader: Mechanical's IronPython reports the reader's
# normal end-of-file StopIteration as an error (the file has no quoted fields)
with open(TABLE) as fh:
    lines = [ln.strip() for ln in fh.read().splitlines() if ln.strip()]
cols = lines[0].split(",")
rows = [dict(zip(cols, ln.split(","))) for ln in lines[1:]]
for r in rows:
    for c in ("hColumn", "TawColumn"):
        if r[c] not in header:
            fail("time table references column %s, which %s does not have" % (r[c], FACES))
say("Faces file: %d cases; time table: %d rows, t = %s to %s s"
    % (len(hCols), len(rows), rows[0]["Time_s"], rows[-1]["Time_s"]))

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
# context-menu code (DSPages/Python/context_menu_visibility.py) assumes one is,
# and the script editor's debugger pauses on the error it throws when none is.
tree = ExtAPI.DataModel.Tree


def keep(obj=None):
    tree.Activate(obj if obj is not None else analysis)


keep()
kids = list(analysis.Children)
for g in [c for c in kids if c.DataModelObjectCategory == DataModelObjectCategory.ImportedLoadGroup]:
    fs = g.GetExternalDataFiles()
    paths = [str(fs[i].FilePath).lower() for i in range(fs.Count)]
    if FACES.lower() in paths:
        say("Deleting earlier Imported Load group '%s'" % g.Name)
        g.Delete()
        keep()
for b in [c for c in kids if c.DataModelObjectCategory == DataModelObjectCategory.Radiation]:
    if b.Name == RAD_NAME:
        say("Deleting earlier '%s'" % b.Name)
        b.Delete()
        keep()

# ---------------- external data ----------------
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
for i, c in hCols:
    st.UseColumn(i, VT.HeatTransferCoefficient, "W m^-2 K^-1", c)
for i, c in tawCols:
    st.UseColumn(i, VT.Temperature, "K", c)
files.Add(edf)

group = analysis.AddImportedLoadExternalData()
keep(group)
group.ImportExternalDataFiles(files)
keep(group)
say("External Data: %s (columns X, Y, Z [m], %d h [W/m^2K], %d T_aw [K])" % (FACES, len(hCols), len(tawCols)))

conv = group.AddImportedConvection()
keep(conv)
conv.Location = skin
keep(conv)

# ---------------- data table ----------------
# GetTableByName returns an Ansys.ACT.Automation.Mechanical.WorksheetTable:
#   .Columns (column names), .Count, tbl[i] -> WorksheetRow, row[colIndex] get/set,
#   .Add(WorksheetRow()). No try/except here on purpose: Mechanical's script editor
#   pauses on every exception, even handled ones, and that freezes the main window.
tbl = conv.GetTableByName("")
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


iFilm = find("film", bad=("scale", "offset"))
iAmb  = find("ambient", bad=("scale", "offset"))
iTime = find("analysis time")
iFilmS = find("film", "scale")
iAmbS  = find("ambient", "scale")
if None in (iFilm, iAmb, iTime, iFilmS, iAmbS):
    fail("could not identify film / ambient / analysis-time / scale columns from the list above; "
         "paste the 'col' lines back")


def choices(col):
    # the table's own dropdown entries for a column, keyed by the bare CSV column name
    # (an entry looks like 'File1:h_M8_a0', possibly with extra decoration)
    items = [str(x) for x in tbl.DropDownItems(col)]
    say("  dropdown col %d: %s" % (col, items))
    out = {}
    for it in items:
        out[it.split(":")[-1].split("@")[0].strip()] = it
    return out


filmItems = choices(iFilm)
ambItems = choices(iAmb)
missing = [c for _, c in hCols if c not in filmItems] + [c for _, c in tawCols if c not in ambItems]
if missing:
    fail("dropdowns have no entry for %s; see the dropdown lines in %s" % (missing, LOG))
timeIsQty = "Quantity" in type(row0[iTime]).__name__
scaleIsQty = "Quantity" in type(row0[iFilmS]).__name__


def num(x, unit, isQty):
    if isQty:
        return Quantity(x, unit)
    return x


while tbl.Count < len(rows):
    tbl.Add(Ansys.ACT.Automation.Mechanical.WorksheetRow())
for k in range(len(rows)):
    r = rows[k]
    row = tbl[k]
    row[iFilm] = filmItems[r["hColumn"]]
    row[iAmb] = ambItems[r["TawColumn"]]
    row[iTime] = num(float(r["Time_s"]), "s", timeIsQty)
    row[iFilmS] = num(float(r["hScale"]), "", scaleIsQty)
    row[iAmbS] = num(float(r["TawScale"]), "", scaleIsQty)

# read back the first and last rows to prove the writes landed
for k in (0, len(rows) - 1):
    got = tbl[k]
    say("row %d: time %s | film %s x %s | ambient %s x %s" % (k, got[iTime], got[iFilm], got[iFilmS], got[iAmb], got[iAmbS]))
if tbl.Count != len(rows):
    fail("table has %d rows, expected %d" % (tbl.Count, len(rows)))

say("Imported Convection: %d rows set up on %s (not imported yet: set the mapping options first)" % (len(rows), skin.Name))

# ---------------- radiation ----------------
if ADD_RADIATION:
    rad = analysis.AddRadiation()
    keep(rad)
    rad.Name = RAD_NAME
    rad.Location = skin
    rad.Emissivity.Output.DiscreteValues = [Quantity(EMISSIVITY, "")]
    rad.AmbientTemperature.Output.DiscreteValues = [Quantity(T_RADIATE_K, "K")]
    say("Radiation: emissivity %.2f to %.0f K on %s" % (EMISSIVITY, T_RADIATE_K, skin.Name))

# ---------------- analysis settings: a step per case window ----------------
if SETUP_STEPS:
    ends = []
    for k in range(1, len(rows)):
        if rows[k]["CaseIdx"] != rows[k - 1]["CaseIdx"]:
            ends.append(float(rows[k]["Time_s"]))
    ends.append(float(rows[-1]["Time_s"]))
    s = analysis.AnalysisSettings
    keep(s)
    s.NumberOfSteps = len(ends)
    # end times must increase step to step; new steps default to 1, 2, 3 ... s, so set
    # them last-first (otherwise Mechanical silently rejects step 1 = 71.5 s while step 2 = 2 s)
    for n in range(len(ends), 0, -1):
        s.SetStepEndTime(n, Quantity(ends[n - 1], "s"))
    for n in range(1, len(ends) + 1):
        s.SetAutomaticTimeStepping(n, AutomaticTimeStepping.On)
        s.SetDefineBy(n, TimeStepDefineByType.Time)
        s.SetInitialTimeStep(n, Quantity(T_INITIAL_S, "s"))
        s.SetMinimumTimeStep(n, Quantity(T_MIN_S, "s"))
        s.SetMaximumTimeStep(n, Quantity(T_MAX_S, "s"))
    say("Analysis Settings: %d steps ending at %s s" % (len(ends), ", ".join(["%.1f" % e for e in ends])))

say("Done. By hand on Imported Convection: Tabular Loading = Ramped, mapping = Nearest, "
    "Outside Option = Ignore; then right-click > Import Load and check contours at a few Analysis Times.")
