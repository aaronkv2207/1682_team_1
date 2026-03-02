import numpy as np
g=9.81
#static surfaces
S_w=39
S_htail=9.1
S_vtail=4.46

#control surfaces both wings included)
S_flaps=10.42
S_ailerons=3.08
S_elevators=3.25
S_rudder=3.06

def W_batt(E_batt, E_m_batt):
    return (E_batt * g) / E_m_batt

def W_motor(Q_max, Q_m_motor):
    return (Q_max * g) / Q_m_motor

def W_wiring(W_wire):
    return W_wire

def W_elec(W_electronics):
    return W_electronics

#W_total = W_power_system(battery=W_batt(E, Em),motor=W_motor(Q, Qm),wiring=W_wiring(Ww),electronics=W_elec(We))



def W_power_system(**pow_comps):
    return sum(pow_comps.values())

def W_air(*args_air):
    #W_wing, W_fuse, W_tail
    def W_wing(*args_wing):
       # W_skin, W_spar, W_act
            #W_skin sclaes with wing area
            #W_act scales with max control surface hinge torques
        return sum(args_wing)
    return sum(args_air)

def W_fix(*args_fix):
    #W_pax, W_cargo
    return sum(args_fix)

def MTOW(W_pow, W_air, W_fix)
    return sum(W_pow, W_air, W_fix)
