import numpy as np
from ambiance import Atmosphere
from conceptual_design import V_STALL, W_S, S, ureg


class CruiseModel:
    def __init__(self, s_ref, weight, v_cruise, h_cruise, AR, e, Cd0) -> None:
        self.s_ref = s_ref * ureg("m^2")
        self.weight = weight * ureg("N")  # newtons
        self.v_cruise = v_cruise.magnitude * ureg("m/s")
        self.h_cruise = h_cruise
        self.density = Atmosphere(h=h_cruise.magnitude).density[0] * ureg("kg/m^3")
        self.AR = AR
        self.e = e
        self.q = 0.5 * self.density * (self.v_cruise**2)
        # self.thrust = thrust
        self.Cd0 = Cd0

    def cl(self):
        L = self.weight  # assumes I have weight as a function of time
        return L / (self.q * self.s_ref)

    def cd_induced(self):
        # see eqn: https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/induced.html
        return (self.cl() ** 2) / (np.pi * self.AR * self.e)

    def cd_parasitic(self):
        return self.Cd0

    def cd_total(self):
        # Current model will not account for viscous drag (not significant in cd calculations)
        return self.cd_induced().magnitude + self.cd_parasitic()

    def drag_total(self):
        return self.cd_total() * self.q * self.s_ref


# Runner script
if __name__ == "__main__":
    # Design Variables
    e = 0.7  # TODO: determine if this value is a reasonable guess
    Cd0 = parastic_drag()
    AR = 10  # placeholder
    # TODO: Update with actual values by calling relevant class
    cruise_cls = CruiseModel(
        s_ref=10,  # TODO: update to varied model
        weight=1000,  # TODO: update to varied model
        v_cruise=V_CRUISE,  # # TODO: update to varied model; # should sweep over an array of cruise values
        h_cruise=AircraftConfig.alt_cruise_m,
        AR=AR,
        e=e,
        Cd0=Cd0,
    )

    CD_total = cruise_cls.cd_total()
    L_over_D = cruise_cls.cl() / CD_total
