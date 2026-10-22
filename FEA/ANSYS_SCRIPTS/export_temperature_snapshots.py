# export_temperature_snapshots.py
# Run inside ANSYS Mechanical from the Scripting Shell (no editor debugger):
#   exec(open(r"C:\Users\Kevin Kyle\Documents\Hypersonics\FEA\ANSYS_SCRIPTS\export_temperature_snapshots.py").read())
#
# Exports the temperature field at each time in TIMES to its own text file
# (Node Number, X, Y, Z, Temperature). Uses a separate "Snapshot export"
# Temperature result so your own results are left alone.
# Needs File > Options > Export > Text File Export > Include Node Location = Yes.
import os

OUT   = r"C:\Users\Kevin Kyle\Documents\Hypersonics\FEA\Exports\snapshots"
TIMES = [3.2, 31.252, 197.6, 526.0, 851.16]      # [s] stored result times
NAME  = "Snapshot export"

if not os.path.isdir(OUT):
    os.makedirs(OUT)
LOG = os.path.join(OUT, "export_log.txt")
open(LOG, "w").close()


def say(msg):
    print(msg)
    with open(LOG, "a") as fh:
        fh.write(msg + "\n")


analysis = None
for a in Model.Analyses:
    if analysis is None and a.AnalysisType == AnalysisType.Transient and a.PhysicsType == PhysicsType.Thermal:
        analysis = a
sol = analysis.Solution
tree = ExtAPI.DataModel.Tree

res = None
for c in sol.Children:
    if c.Name == NAME:
        res = c
if res is None:
    res = sol.AddTemperature()
    res.Name = NAME
tree.Activate(res)
res.By = SetDriverStyle.Time

for t in TIMES:
    res.DisplayTime = Quantity(t, "s")
    res.EvaluateAllResults()
    path = os.path.join(OUT, "T_%07.2fs.txt" % t)
    res.ExportToTextFile(path)
    say("t = %8.2f s  ->  %s" % (t, path))

say("Done: %d snapshots in %s" % (len(TIMES), OUT))
