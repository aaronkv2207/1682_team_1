### Wing Sizing Workflow

Run runner.py to create geometry case files and run at your desired specified operating conditions. (Note, this will need to be modified when running JVL--will need to impose trim, Cmalpha = 0, constraint and implicitly define deflection as a free-variable.)

Determine whether the constraints below are satisfied: 
- Takeoff: How much lift can the wing produce under, reasonable, maximum augmentation? --> determine takeoff velocity (1.1 x Vstall) & corresponding x_TO for each S_ref tested
- Landing: Can the aircraft be trimmed and controlled at high lift without exceeding control deflection limits?

Finally, post-process saved-outputs for each of the cases at takeoff and at landing.