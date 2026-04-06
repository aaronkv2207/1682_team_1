### Wing Sizing Workflow

Run `runner.py` to create geometry case files and run at your specified operating conditions. (Note, this will need to be modified when running JVL--will need to impose trim, Cmalpha = 0, constraint and implicitly define deflection as a free-variable.)

Determine whether the constraints below are satisfied: 
- Takeoff (V1+V2 sizing): How much lift can the wing produce under, reasonable, maximum augmentation? --> determine corresponding x_TO for each S_ref tested.
- Cruise (V1+V2 sizing): Characterize trade with cruise performance.
- Landing (V2 sizing): Can the aircraft be trimmed and controlled at high lift without exceeding control deflection limits? (Will primarily be assessed with aft-end control surfaces in later stages.)

Finally, post-process saved outputs in `Python/aero_workspace/design_trades` for each of the cases.