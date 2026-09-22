from ansys.mechanical.core import launch_mechanical

mech = launch_mechanical(batch=False)
mech.wait_till_mechanical_is_ready()

script = r"""
print("Hello from Mechanical UI!")
print("Running inside Mechanical's Python environment.")
"""
mech.run_python_script(script)

input("Press Enter to close Mechanical...")