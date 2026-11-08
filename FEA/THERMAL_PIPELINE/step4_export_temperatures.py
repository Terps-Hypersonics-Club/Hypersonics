# step4_export_temperatures.py   (THERMAL_PIPELINE, see README.md)
# Run inside ANSYS Mechanical from the Scripting Shell (no editor debugger), after
# the solve. Set OUT and TIMES below first, then:
#   exec(open(r"C:\...\Hypersonics\FEA\THERMAL_PIPELINE\step4_export_temperatures.py").read())
#
# Exports the temperature field (T_*.txt) and the total heat flux field (q_*.txt) at
# each time in TIMES (Node Number, X, Y, Z, value). Uses separate "Snapshot export"
# result objects so your own results are left alone.
# Needs File > Options > Export > Text File Export > Include Node Location = Yes.
import os

OUT   = r"C:\Users\Kevin Kyle\Documents\Hypersonics\FEA\THERMAL_PIPELINE\outputs\snapshots"   # CHANGE to your path
END_S = 824.0                                    # [s] last time in the solve
TIMES = [float(t) for t in range(0, int(END_S), 10)] + [END_S]   # every 10 s of flight, plus the end
NAME  = "Snapshot export"
# Mechanical only has results at the times it stored (Output Controls). If a requested
# time is not stored it shows the nearest stored set; the log records which time it used.

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

def result(name, add):
    # reuse our own result object if an earlier run made it, so user results stay untouched
    for c in sol.Children:
        if c.Name == name:
            return c
    r = add()
    r.Name = name
    return r


res = result(NAME, sol.AddTemperature)
flux = result(NAME + " (heat flux)", sol.AddTotalHeatFlux)   # needs Output Controls > Calculate Thermal Flux = Yes
for r in (res, flux):
    tree.Activate(r)
    r.By = SetDriverStyle.Time

for t in TIMES:
    for r, tag in ((res, "T"), (flux, "q")):
        r.DisplayTime = Quantity(t, "s")
        r.EvaluateAllResults()
        used = r.DisplayTime.Value        # the time Mechanical actually shows
        path = os.path.join(OUT, "%s_%06.1fs.txt" % (tag, t))
        r.ExportToTextFile(path)
        say("%s t = %6.1f s (result at %8.3f s), max %.4g  ->  %s" % (tag, t, used, r.Maximum.Value, path))

say("Done: %d snapshots in %s" % (len(TIMES), OUT))
