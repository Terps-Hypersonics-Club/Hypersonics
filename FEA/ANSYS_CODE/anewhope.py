import numpy as np
import pandas as pd
import pyvista as pv
from ansys.mapdl.core import launch_mapdl
from ansys.mapdl.core import LOG
from ansys.mapdl.core.errors import MapdlRuntimeError

LOG.setLevel("DEBUG")
LOG.log_to_file("mylog.log")

# start mapdl
mapdl = launch_mapdl(run_location=r"C:\Users\kdk\Documents\PyMechanical\MAPDL_runs\run10", jobname="thermal_run",loglevel="DEBUG")
print(mapdl)

geometry_path = "Odd_seed_Point_101_round.x_t"
hfluxdata_path = "heatflux.csv"
mesh_path = "mymesh.cdb"

mapdl.clear()
print("adding mesh")
mapdl.prep7() # Ensure you are in the preprocessor
#mapdl.run('TOLLER, 0.0001') # Set tolerance for cleaner import (optional but recommended)
#mapdl.parain(geometry_path) # CORRECT COMMAND for IGES files
# Load the Mechanical mesh you exported
mapdl.input("mymesh.cdb")
print("Mesh loading complete")
#print(mapdl.geometry)

# Define element attributes
# Second-order tetrahedral elements (SOLID187)
#mapdl.prep7()
#mapdl.et(1, "SOLID87") # Orignally SOLID187 but switched to SOLID70 for thermal analysis but SOLID70 is brick elements, instead using tetra thermal elements in the form of SOLID87
mapdl.prep7() # Need another prep because mesh loading kills prep
print("Assigning Material Values")
# Define material properties of structural steel
E = 2e11  # Youngs modulus
NU = 0.3  # Poisson's ratio
CTE = 1.2e-5  # Coeff. of thermal expansion
CP = 500
K = 80
RHO = 7850
mapdl.mp("EX", 1, E)
mapdl.mp("PRXY", 1, NU)
mapdl.mp("ALPX", 1, CTE)
mapdl.mp("KXX", 1, K)      # W/mK
mapdl.mp("C",   1, CP)     # J/kgK
mapdl.mp("DENS",1, RHO)    # kg/m^3

print("Material Values Assigned")

#print("Meshing")
#mapdl.run('SHPP, WARN') # Tells the solver to accept poor element shapes as WARNINGS, not ERRORS.
#mapdl.vmesh("all")
#print("Meshing complete")

# Attempting plotting now!
mapdl.eplot(show_edges=True, smooth_shading=True,
                show_node_numbering=True)


# Save mesh as VTK object
print(mapdl.mesh)
grid = mapdl.mesh.grid  # save mesh as a VTK object

# Import csv file and save data to a NumPy array
hflux_file = pd.read_csv(hfluxdata_path, sep=",", header=None, low_memory=False)
print(hflux_file)
hflux_data = hflux_file.values  # Save data to a NumPy array
nd_hflux_data = hflux_data[1:, 1:].astype(float)  # Change data type to Float

# Map temperature data to FE mesh
# Convert imported data into PolyData format
print("Polyfitting function: ")
wrapped = pv.PolyData(nd_hflux_data[:, :3])  # Convert NumPy array to PolyData format
wrapped["heat flux"] = nd_hflux_data[
    :, 3
]  # Add a scalar variable 'heat flux' to PolyData
print("Polyfitting complete. Starting data mapping: ")


# Perform data mapping
inter_grid = grid.interpolate(
    wrapped,
    sharpness=5,
    radius=0.0001,
    strategy="closest_point",
    progress_bar=True,
)  # Map the imported data to MAPDL grid
inter_grid.plot(show_edges=False)  # Plot the interpolated data on MAPDL grid
hflux_load_val = pv.convert_array(
    pv.convert_array(inter_grid.active_scalars)
)  # Save temperatures interpolated to each node as NumPy array
node_num = inter_grid.point_data["ansys_node_num"]  # Save node numbers as NumPy array


print("Data mapping complete, preparing to enter solution mode: ")
# Enter /SOLU processor to apply loads and BCs
mapdl.finish() # Exit from processor (PREP7)
mapdl.slashsolu() # Enter into solution mode
mapdl.tunif(293.15)
mapdl.antype("trans")        # transient
mapdl.trnopt("full")         # full transient (or "on" depending on version)
mapdl.kbc(0)                 # ramped loads (1 = stepped)
mapdl.time(10.0)             # total time of analysis
mapdl.nsubst(50, 50, 1)      # 50 time steps (min,max,step change factor)

print("Solution mode entered, assigning heat flux at each node: ")
# Enter non-interactive mode to assign heatflux at each node using imported data
with mapdl.non_interactive:
    for node, hflux in zip(node_num,hflux_load_val):
        mapdl.bf(node, "HGEN", hflux) # HGEN is defined as heat generation by nodes... other option is HEAT which is just heat flux without consideration of nodal area 

print("Heat Flux assigned, solving the model now")

try:
    output = mapdl.solve()
    print(output)
except MapdlRuntimeError as e:
    # MAPDL printed an ERROR about element shape, but often still proceeds.
    print("MAPDL reported element shape errors during SOLVE.")
    print("Continuing anyway – check .out/.rst if you’re worried about quality.")
    print(e)

print("Model solved. Entering post-processing")
# Enter post-processor
mapdl.post1()
# number of substeps
# ===== SETTINGS YOU SHOULD TWEAK =====
print("Detecting available result sets...")
available_steps = []
i = 1
while True:
    try:
        mapdl.set(1, i)
        available_steps.append(i)
        i += 1
    except MapdlRuntimeError:
        break

nsteps = len(available_steps)
print(f"Found {nsteps} result substeps")   # <-- set this to the number of substeps you actually solved
gif_name = "temp_pyvista.gif"
# =====================================

print("Building PyVista temperature animation...")

# Ensure full model is selected
mapdl.allsel()

# First: scan all steps to get global min/max temperature for consistent colors
tmin, tmax = float("inf"), float("-inf")
for i in range(1, nsteps + 1):
    mapdl.set(1, i)
    temp = mapdl.post_processing.nodal_temperature()
    tmin = min(tmin, float(temp.min()))
    tmax = max(tmax, float(temp.max()))

print(f"Global temperature range: {tmin:.2f} to {tmax:.2f}")

# Get the PyVista grid from MAPDL mesh
grid = mapdl.mesh.grid  # unstructured grid with nodes/elements

# Optional: pick a nice camera orientation (you can change this later)
cpos = "xy"  # or "xz", "yz", or a custom camera tuple

# Create PyVista plotter for off-screen rendering
plotter = pv.Plotter(off_screen=True)
plotter.open_gif(gif_name)

for i in range(1, nsteps + 1):
    print(f"Rendering frame {i}/{nsteps}...")
    mapdl.set(1, i)
    temp = mapdl.post_processing.nodal_temperature()

    # Attach temperature to the grid as point data
    grid.point_data["TEMP"] = temp

    # Draw mesh with temperature field
    plotter.add_mesh(
        grid,
        scalars="TEMP",
        clim=[tmin, tmax],
        scalar_bar_args={"title": "Temperature"},
        show_edges=False,
        smooth_shading=True,
    )

    plotter.add_text(f"Time step {i}/{nsteps}", font_size=12)

    # Fix camera position
    plotter.camera_position = cpos

    # Write this frame to the GIF
    plotter.write_frame()

    # Clear actors for the next frame
    plotter.clear()

plotter.close()
print(f"PyVista GIF saved as: {gif_name}")

mapdl.set(1, nsteps)

# 1) Get MAPDL mesh as PyVista grid
grid = mapdl.mesh.grid

# 2) Get nodal temperatures from MAPDL and attach to grid
temp = mapdl.post_processing.nodal_temperature()
grid.point_data["TEMP"] = temp

# 3) Sample the temperature field at your centroid locations
#    `wrapped` should contain the centroid points you used for the original flux
sampled = wrapped.sample(grid)   # same geometry as wrapped, with TEMP interpolated

centroid_temps = sampled.point_data["TEMP"]  # numpy array, one value per centroid

print("Centroid temperatures shape:", centroid_temps.shape)

# Optional: save centroids + temps to CSV for post processing
centroids = wrapped.points  # xyz of each centroid
df = pd.DataFrame({
    "x": centroids[:, 0],
    "y": centroids[:, 1],
    "z": centroids[:, 2],
    "T": centroid_temps,
})
df.to_csv("centroid_temperatures.csv", index=False)
print("Saved centroid temperatures to centroid_temperatures.csv")

print("Launching interactive 3D viewer with time slider...")

mapdl.post1()
mapdl.allsel()

# Detect available result sets
from ansys.mapdl.core.errors import MapdlRuntimeError

available_steps = []
i = 1
while True:
    try:
        mapdl.set(1, i)
        available_steps.append(i)
        i += 1
    except MapdlRuntimeError:
        break

nsteps = len(available_steps)
print(f"Found {nsteps} result substeps")

# Global min and max temperature
tmin, tmax = float("inf"), float("-inf")
for step in available_steps:
    mapdl.set(1, step)
    temp = mapdl.post_processing.nodal_temperature()
    tmin = min(tmin, float(temp.min()))
    tmax = max(tmax, float(temp.max()))
print(f"Global temperature range: {tmin:.2f} to {tmax:.2f}")

grid = mapdl.mesh.grid

# Start at first step
mapdl.set(1, available_steps[0])
temp0 = mapdl.post_processing.nodal_temperature()
grid.point_data["TEMP"] = temp0

plotter = pv.Plotter()
actor = plotter.add_mesh(
    grid,
    scalars="TEMP",
    clim=[tmin, tmax],
    scalar_bar_args={"title": "Temperature"},
    show_edges=False,
    smooth_shading=True,
)
plotter.camera_position = "iso"

def update_time(value):
    step = int(round(value))
    print(f"Slider moved to step {step}")
    try:
        mapdl.set(1, step)
        temp = mapdl.post_processing.nodal_temperature()
        grid.point_data["TEMP"] = temp
        actor.mapper.scalar_range = (tmin, tmax)
        plotter.update_scalar_bar_range([tmin, tmax])
        plotter.render()
    except MapdlRuntimeError:
        print(f"No result for step {step}")

plotter.add_slider_widget(
    callback=update_time,
    rng=[available_steps[0], available_steps[-1]],
    value=available_steps[0],
    title="Time step",
    style="modern",
)

plotter.show()


print("Code Finished")
input("Press Enter to Exit APDL")

mapdl.exit()
input("APDL and Program finished")