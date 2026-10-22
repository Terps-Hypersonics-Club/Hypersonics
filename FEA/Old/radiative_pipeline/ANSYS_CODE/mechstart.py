from ansys.mechanical.core import launch_mechanical

# Launch the Mechanical instance (with GUI so you can see it)
mech = launch_mechanical(batch=False)
mech.wait_till_mechanical_is_ready()

print("Mechanical is open and ready!")
input("Press Enter to close Mechanical...")