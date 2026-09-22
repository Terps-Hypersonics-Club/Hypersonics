from ansys.mechanical.core import launch_mechanical

# Launch Mechanical with GUI
mech = launch_mechanical(batch=False)
mech.wait_till_mechanical_is_ready()

proj_directory = mech.run_python_script("ExtAPI.DataModel.Project.ProjectDirectory")
print(proj_directory)
proj_directory = mech.project_directory
print(proj_directory)

input("Press Enter to close Mechanical...")