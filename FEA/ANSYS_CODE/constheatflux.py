from ansys.mechanical.core import Mechanical
import os

model_path = r"Odd_seed_Point_101_round.STEP"
csv_path   = r"heatflux_mach5dot5.csv"   # columns: FaceID, HeatFlux or Flux_W_per_m2
log_path = os.path.join(os.getcwd(), "PyMechanical_log.txt")
if os.path.exists(log_path):
    os.remove(log_path)

mech = Mechanical()
mech.wait_till_mechanical_is_ready()

script = r'''
import csv, traceback, os
from System.IO import File

LOG_PATH   = r"__LOG__"
MODEL_PATH = r"__MODEL__"
CSV_PATH   = r"__CSV__"
MAC_PATH   = os.path.join(ExtAPI.DataModel.Project.ProjectDirectory, "apply_flux.mac")

def write_log(msg):
    try:
        with open(LOG_PATH, "a") as f:
            f.write(msg + "\n")
    except:
        pass
    ExtAPI.Log.WriteMessage(msg)

write_log("=== Transient heat-flux simulation started ===")

try:
    write_log("Project dir: " + ExtAPI.DataModel.Project.ProjectDirectory)

    # 1) Geometry import
    if not File.Exists(MODEL_PATH):
        write_log("Geometry file not found: " + MODEL_PATH)
        raise RuntimeError("Geometry missing")
    gi = Model.GeometryImportGroup.AddGeometryImport()
    gi.Import(MODEL_PATH)
    write_log("Geometry imported successfully.")

    # 2) Analysis
    analysis = Model.AddTransientThermalAnalysis()
    write_log("Transient thermal analysis created.")

    # 3) Mesh
    write_log("Generating mesh...")
    Model.Mesh.GenerateMesh()
    write_log("Mesh generated successfully.")

    # 3.1) MeshData handle and capability dump
    try:
        mesh = ExtAPI.DataModel.MeshDataByName("Global")
    except:
        mesh = getattr(Model, "MeshData", None)
    if mesh is None:
        raise RuntimeError("Mesh topology API is unavailable on this build")

    try:
        mesh.BuildGeoMeshCache()
        write_log("Built geo-mesh cache")
    except Exception as e:
        write_log("BuildGeoMeshCache failed: " + str(e))

    try:
        mesh_attrs = sorted([x for x in dir(mesh) if not x.startswith("_")])
        write_log("MeshData capabilities: " + ", ".join(mesh_attrs))
    except Exception as e:
        write_log("Could not list MeshData attrs: " + str(e))

    # 4) Read CSV (FaceID, HeatFlux or Flux_W_per_m2)
    if not File.Exists(CSV_PATH):
        write_log("CSV not found: " + CSV_PATH)
        raise RuntimeError("CSV missing")

    face_flux = []
    with open(CSV_PATH, "r") as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            rid = row.get("FaceID")
            q   = row.get("HeatFlux") or row.get("Flux_W_per_m2")
            if rid and q:
                try:
                    face_flux.append((int(float(rid)), float(q)))
                except:
                    pass

    if not face_flux:
        raise RuntimeError("No valid rows found. Expect headers: FaceID, HeatFlux or Flux_W_per_m2")
    write_log("Read {} face rows from CSV".format(len(face_flux)))

    # 5) helpers to fetch node ids for a geometric face id
    def nodes_for_face_refid(ref_id):
        # Preferred: bulk getter
        get_nodes_bulk = getattr(mesh, "GetNodeIdsFromRegionIds", None)
        if get_nodes_bulk:
            try:
                arr = get_nodes_bulk([int(ref_id)])
                if arr and len(arr) > 0 and arr[0] is not None:
                    return [int(n) for n in arr[0]]
            except Exception as e:
                write_log("GetNodeIdsFromRegionIds failed for FaceID {}: {}".format(ref_id, e))
        # Fallback: region object
        get_region = getattr(mesh, "MeshRegionById", None)
        if get_region:
            try:
                reg = get_region(int(ref_id))
                node_ids = getattr(reg, "NodeIds", None)
                if node_ids is not None:
                    return [int(n) for n in node_ids]
                nodes_coll = getattr(reg, "Nodes", None)
                if nodes_coll is not None:
                    return [int(n.Id) for n in nodes_coll]
            except Exception as e:
                write_log("MeshRegionById failed for FaceID {}: {}".format(ref_id, e))
        # Nothing found
        return []

    # APDL writer utilities
    def write_node_select(mf, node_ids, chunk=64):
        if not node_ids:
            mf.write("nsel,none\n")
            return
        it = iter(node_ids)
        first = next(it, None)
        if first is None:
            mf.write("nsel,none\n")
            return
        mf.write("nsel,s,node,,{}\n".format(first))
        acc = []
        for n in it:
            acc.append(n)
            if len(acc) >= chunk:
                mf.write("nsel,a,node,,{}\n".format(",".join(str(x) for x in acc)))
                acc = []
        if acc:
            mf.write("nsel,a,node,,{}\n".format(",".join(str(x) for x in acc)))

    # 6) Build the macro that uses ESURF to stamp per-element-face flux
    mapped = 0
    missed = 0
    with open(MAC_PATH, "w") as mf:
        mf.write("/prep7\n")
        for i, (face_id, q) in enumerate(face_flux, 1):
            try:
                node_ids = nodes_for_face_refid(face_id)
                if not node_ids:
                    missed += 1
                    continue
                mf.write("nsel,all\n")
                mf.write("esel,all\n")
                write_node_select(mf, node_ids)
                mf.write("esln,s\n")
                mf.write("sfe,all,1,HFLUX,,{:.12e}\n".format(q))
                mf.write("esurf\n")
                mapped += 1
                if i % 100 == 0:
                    mf.flush()
                    write_log("processed {} of {} faces".format(i, len(face_flux)))
                    ExtAPI.DoEvents()
            except Exception as e:
                missed += 1
                write_log("Face {} macro build error: {}".format(face_id, e))
        mf.write("! end macro\n")

    write_log("Macro written to: " + MAC_PATH.replace("\\", "/"))
    write_log("Faces processed: {}, faces missed: {}".format(mapped, missed))

    # 7) Execute macro once before solve
    cmd = analysis.AddCommandSnippet()
    cmd.Name  = "ApplyPerFaceFlux_ESURF"
    cmd.Input = r"""
/prep7
*use,'{}'
""".format(MAC_PATH.replace("\\","/"))

    # 8) Time stepping
    try:
        from Ansys.Mechanical.DataModel.Enums import AutomaticTimeStepping
        s = analysis.AnalysisSettings
        s.AutomaticTimeStepping = AutomaticTimeStepping.On
        s.NumberOfSteps   = 100
        s.StepEndTime     = Quantity(1.0, "s")
        s.InitialTimeStep = Quantity(0.01, "s")
        s.MinimumTimeStep = Quantity(0.001, "s")
        s.MaximumTimeStep = Quantity(0.1, "s")
        write_log("Adaptive time stepping set")
    except:
        s = analysis.AnalysisSettings
        s.NumberOfSteps = 100
        s.StepEndTime   = Quantity(1.0, "s")
        s.TimeStep      = Quantity(0.01, "s")
        write_log("Fixed time stepping set")

    # 9) Solve
    write_log("Solving transient thermal case...")
    analysis.Solve(True)
    write_log("=== Solve complete ===")

    # 10) Post
    temp = analysis.Solution.AddTemperature()
    analysis.Solution.EvaluateAllResults()
    write_log("Temperature evaluated successfully.")

except Exception as e:
    write_log("CRITICAL ERROR: " + str(e))
    write_log(traceback.format_exc())

write_log("=== Script finished ===")
'''

script = script.replace("__LOG__", log_path)\
               .replace("__MODEL__", os.path.abspath(model_path))\
               .replace("__CSV__", os.path.abspath(csv_path))

mech.run_python_script(script)
print(f"Log saved to: {log_path}")
input("Program finished. Press Enter to exit...")
