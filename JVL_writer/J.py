import os
import subprocess
import tempfile
import warnings
from pathlib import Path
from typing import Dict, List, Union
import aerosandbox as asb
from aerosandbox import AVL
from typing import List
from pathlib import Path
import aerosandbox.numpy as np


class JWing(asb.Wing):
    def __init__(
        self, name, xsecs, JetParam=None, symmetric=True, JetSpacing=None, **kwargs
    ):
        super().__init__(name=name, xsecs=xsecs, symmetric=symmetric, **kwargs)
        self.JetParam = JetParam
        if JetSpacing is None:
            self.JetSpacing = {
                "Nujet": 0.2,
                "Cspu": 0.0,
                "Nwjet": 12,
                "Cewsp": -2.0,
            }
        else:
            self.JetSpacing = JetSpacing


class JVL(AVL):
    def __init__(
        self,
        airplane,
        op_point,
        xyz_ref=[0, 0, 0],
        ground_effect=False,
        ground_effect_height=0.0,
        AVL_spacing_parameters=None,
        avl_command=".\\jvl2.20",
    ):
        super().__init__(
            airplane=airplane,
            op_point=op_point,
            # xyz_ref=airplane.xyz_ref,
            ground_effect=ground_effect,
            ground_effect_height=ground_effect_height,
            avl_command=avl_command,
        )

    def write_jvl(
        self,
        filepath,
        CLAF=True,
        j=True,
    ) -> None:
        """
        Writes a .avl file corresponding to this airplane to a filepath.

        For use with the AVL vortex-lattice-method aerodynamics analysis tool by Mark Drela at MIT.
        AVL is available here: https://web.mit.edu/drela/Public/web/avl/

        Args:
            filepath: filepath (including the filename and .avl extension) [string]
                If None, this function returns the .avl file as a string.

        Returns: None

        """

        def clean(s):
            """
            Removes leading and trailing whitespace from each line of a multi-line string.
            """
            return "\n".join([line.strip() for line in s.split("\n")])

        airplane = self.airplane

        jvl_file = ""

        airplane_options = self.get_options(airplane)

        jvl_file += clean(
            f"""\
            {airplane.name}
            #Mach
            0        ! AeroSandbox note: This is overwritten later to match the current OperatingPoint Mach during the AVL run.
            #IYsym   IZsym   Zsym
            0       {1 if self.ground_effect else 0}   {self.ground_effect_height}
            #Sref    Cref    Bref
            {airplane.s_ref} {airplane.c_ref} {airplane.b_ref}
            #Xref    Yref    Zref
            {self.xyz_ref[0]} {self.xyz_ref[1]} {self.xyz_ref[2]}
            # CDp
            {airplane_options["profile_drag_coefficient"]}
            """
        )

        control_surface_counter = 0
        airfoil_counter = 0

        for wing in airplane.wings:
            wing_options = self.get_options(wing)

            spacing_line = f"{wing_options['chordwise_resolution']}   {self.AVL_spacing_parameters[wing_options['chordwise_spacing']]}"
            if wing_options["wing_level_spanwise_spacing"]:
                spacing_line += f"   {wing_options['spanwise_resolution']}   {self.AVL_spacing_parameters[wing_options['spanwise_spacing']]}"
                if isinstance(wing, JWing):
                    spacing_line += f"   {wing.JetSpacing['Nujet']}   {wing.JetSpacing['Cspu']}   {wing.JetSpacing['Nwjet']}   {wing.JetSpacing['Cewsp']}"

            jvl_file += clean(
                f"""\
                #{"=" * 79}
                SURFACE
                {wing.name}
                #Nchord  Cspace  [ Nspan Sspace  Nujet Cusp  Nwjet Cwsp ]
                {spacing_line}
                
                """
            )

            if wing_options["component"] is not None:
                jvl_file += clean(
                    f"""\
                    COMPONENT
                    {wing_options["component"]}
                        
                    """
                )

            if wing.symmetric:
                jvl_file += clean(
                    """\
                    YDUPLICATE
                    0
                        
                    """
                )

            if wing_options["no_wake"]:
                jvl_file += clean(
                    """\
                    NOWAKE
                    
                    """
                )

            if wing_options["no_alpha_beta"]:
                jvl_file += clean(
                    """\
                    NOALBE
                    
                    """
                )

            if wing_options["no_load"]:
                jvl_file += clean(
                    """\
                    NOLOAD
                    
                    """
                )

            if j:
                if isinstance(wing, JWing) and wing.JetParam is not None:
                    JetParam = wing.JetParam
                    jvl_file += clean(
                        f"""\
                            JETPARAM
                            #hdisk   fh   djet0   djet1   djet3 dxdsk dndsk
                            {JetParam.hdisk:.3f} {JetParam.fh:.3f} {JetParam.djet0:.3f} {JetParam.djet1:.3f} {JetParam.djet3:.6f} {JetParam.dxdsk:.6f} {JetParam.dndsk:.6f}
                            """
                    )

            ### Build up a buffer of the control surface strings to write to each section
            control_surface_commands: List[List[str]] = [[] for _ in wing.xsecs]
            for i, xsec in enumerate(wing.xsecs):
                for surf in xsec.control_surfaces:
                    xhinge = (
                        surf.hinge_point if surf.trailing_edge else -surf.hinge_point
                    )
                    sign_dup = 1 if surf.symmetric else -1

                    command = (
                        clean(
                            f"""\
                            CONTROL
                            #name, gain, Xhinge, XYZhvec, SgnDup
                            {surf.name} 1 {xhinge:.8g} 0 0 0 {sign_dup}
                            """
                        )
                        + "\n"
                    )

                    control_surface_commands[i].append(command)
                    # control_surface_commands[i + 1].append(command)

            ### Write the commands for each wing section
            for i, xsec in enumerate(wing.xsecs):
                xsec_options = self.get_options(xsec)

                xsec_def_line = f"{xsec.xyz_le[0]:.8g} {xsec.xyz_le[1]:.8g} {xsec.xyz_le[2]:.8g} {xsec.chord:.8g} {xsec.twist:.8g}"
                if not wing_options["wing_level_spanwise_spacing"]:
                    xsec_def_line += f"   {xsec_options['spanwise_resolution']}   {self.AVL_spacing_parameters[xsec_options['spanwise_spacing']]}"

                if xsec_options["cl_alpha_factor"] is None:
                    claf_line = f"{1 + 0.77 * xsec.airfoil.max_thickness()}  # Computed using rule from avl_doc.txt"
                else:
                    claf_line = f"{xsec_options['cl_alpha_factor']}"

                af_filepath = Path(str(filepath) + f".af{airfoil_counter}")
                airfoil_counter += 1
                xsec.airfoil.repanel(50).write_dat(
                    filepath=af_filepath, include_name=True
                )

                jvl_file += clean(
                    f"""\
                    #{"-" * 50}
                    SECTION
                    #Xle    Yle    Zle     Chord   Ainc  [Nspanwise   Sspace]
                    {xsec_def_line}
                    
                    AFIL
                    {af_filepath.name.split("/")[-1]}
                    
                    """
                )
                if CLAF:
                    jvl_file += clean(
                        f"""\
                        CLAF
                        {claf_line}
                        
                        """
                    )
                for control_surface_command in control_surface_commands[i]:
                    jvl_file += control_surface_command

                # **JETCONTROL and JETPARAM Handling**
                if isinstance(xsec, WingJSec):
                    if j:
                        if JetParam is None and len(xsec.JetControls) > 0:
                            raise ValueError(
                                "JetControl defined without JetParam in JWing section."
                            )
                        for param in xsec.JetControls:
                            jvl_file += (
                                clean(
                                    f"""\
                                    JETCONTROL
                                    #Djet  1.0   1.0  ! name, gain, SgnDup
                                    {param.jet_name} {param.gain:.3f} {param.sgn_dup:.3f}
                                    """
                                )
                                + "\n"
                            )

        filepath = Path(filepath)
        for i, fuse in enumerate(airplane.fuselages):
            fuse_filepath = Path(str(filepath) + f".fuse{i}")
            self.write_avl_bfile(fuselage=fuse, filepath=fuse_filepath)
            fuse_options = self.get_options(fuse)

            jvl_file += clean(
                f"""\
                #{"=" * 50}
                BODY
                {fuse.name}
                {fuse_options["panel_resolution"]} {self.AVL_spacing_parameters[fuse_options["panel_spacing"]]}
                
                BFIL
                {fuse_filepath.name}
                
                TRANSLATE
                0 {np.mean([x.xyz_c[1] for x in fuse.xsecs]):.8g} 0
                
                """
            )

        if filepath is not None:
            with open(filepath, "w+") as f:
                f.write(jvl_file)

    def run(
        self,
        run_command: str = None,
        trim_Cm_to_zero: bool = False,
        trim_variable=None,
        flap_deflections=None,
        blowing=None,
        xyz_cg=[None],
    ) -> Dict[str, float]:
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)

            if self.working_directory is not None:
                directory = Path(self.working_directory)

            output_filename = "output.txt"

            # Write JVL geometry
            airplane_file = "airplane.jvl"
            self.write_jvl(directory / airplane_file)

            # Build keystroke script
            keystroke_file_contents = self._default_keystroke_file_contents(
                trim_Cm_to_zero, trim_variable, flap_deflections, blowing, xyz_cg
            )
            if run_command is not None:
                keystroke_file_contents += [run_command]
            keystroke_file_contents += [
                "x",
                "st",
                f"{output_filename}",
                "o",
                "",
                "",
                "quit",
            ]

            keystrokes = "\n".join(keystroke_file_contents)

            command = [self.avl_command, airplane_file]

            # Execute solver
            try:
                proc = subprocess.Popen(
                    command,
                    cwd=directory,
                    stdin=subprocess.PIPE,
                    # stdout=None if self.verbose else subprocess.DEVNULL,
                    # stderr=None if self.verbose else subprocess.DEVNULL,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )

                outs, errs = proc.communicate(input=keystrokes, timeout=self.timeout)
                # print("STDOUT:", outs)
                # print("STDERR:", errs)
                return_code = proc.poll()

            except subprocess.TimeoutExpired:
                proc.kill()
                outs, errs = proc.communicate()

                warnings.warn(
                    "JVL run timed out. Increase timeout if needed.",
                    stacklevel=2,
                )

            # Read output file
            try:
                with open(directory / output_filename, "r") as f:
                    output_data = f.read()
            except FileNotFoundError:
                raise FileNotFoundError(
                    "JVL did not produce an output file. Check executable path or geometry."
                )

            # Parse results
            res = self.parse_unformatted_data_output(
                output_data, data_identifier=" =", overwrite=False
            )

            # Normalize keys
            for key_to_lowerize in ["Mach"]:
                res[key_to_lowerize.lower()] = res.pop(key_to_lowerize)

            for key in list(res.keys()):
                if "tot" in key:
                    res[key.replace("tot", "")] = res.pop(key)

            # Compute dimensional forces
            q = self.op_point.dynamic_pressure()
            S = self.airplane.s_ref
            b = self.airplane.b_ref
            c = self.airplane.c_ref

            res["velocity"] = self.op_point.velocity
            res["p"] = res["pb/2V"] * (2 * self.op_point.velocity / b)
            res["q"] = res["qc/2V"] * (2 * self.op_point.velocity / c)
            res["r"] = res["rb/2V"] * (2 * self.op_point.velocity / b)

            res["L"] = q * S * res["CL"]
            res["Y"] = q * S * res["CY"]
            res["D"] = q * S * res["CD"]

            res["l_b"] = q * S * b * res["Cl"]
            res["m_b"] = q * S * c * res["Cm"]
            res["n_b"] = q * S * b * res["Cn"]

            try:
                res["Clb Cnr / Clr Cnb"] = (
                    res["Clb"] * res["Cnr"] / (res["Clr"] * res["Cnb"])
                )
            except ZeroDivisionError:
                res["Clb Cnr / Clr Cnb"] = np.nan

            # Force vectors
            res["F_w"] = [-res["D"], res["Y"], -res["L"]]

            res["F_b"] = self.op_point.convert_axes(
                *res["F_w"], from_axes="wind", to_axes="body"
            )

            res["F_g"] = self.op_point.convert_axes(
                *res["F_b"], from_axes="body", to_axes="geometry"
            )

            # Moment vectors
            res["M_b"] = [res["l_b"], res["m_b"], res["n_b"]]

            res["M_g"] = self.op_point.convert_axes(
                *res["M_b"], from_axes="body", to_axes="geometry"
            )

            res["M_w"] = self.op_point.convert_axes(
                *res["M_b"], from_axes="body", to_axes="wind"
            )

            return res

    def _default_keystroke_file_contents(
        self,
        trim_Cm_to_zero: bool = False,
        trim_variable=None,
        flap_deflections: dict | None = None,
        blowing: dict | None = None,
        xyz_cg: tuple[float, float, float] | None = None,
    ) -> List[str]:
        """_summary_

        Args:
            trim_Cm_to_zero (bool, optional): Defaults to False.
            trim_variable (str, optional): Defaults to "d6".
            blowing (dict, optional): Dictionary mapping jet variables to their respective magnitudes. Defaults to None.
            xyz_cg ([float, float, float]): xyz_cg
        Returns:
            List[str]: _description_
        """
        run_file_contents = []

        # Disable graphics
        run_file_contents += [
            "plop",
            "g",
            "",
        ]

        run_file_contents += [
            "mass",
            "m.mass",
            "mset 0",
        ]

        # Enter oper mode
        run_file_contents += [
            "oper",
        ]

        # Direct p, q, r to be in body axes, to match ASB convention
        run_file_contents += [
            "o",
            "r",
            "",
        ]

        # Sets Xcg, Ycg, Zcg
        if xyz_cg.any():
            run_file_contents += [
                "c1",
                f"x {xyz_cg[0]}",
                f"y {xyz_cg[1]}",
                f"z {xyz_cg[2]}",
                "",
            ]

        # Set parameters
        run_file_contents += [
            "m",
            f"mn {float(self.op_point.mach())}",
            f"v {float(self.op_point.velocity)}",
            f"d {float(self.op_point.atmosphere.density())}",
            "g 9.81",
            "",
        ]

        # Set analysis state
        p_bar = self.op_point.p * self.airplane.b_ref / (2 * self.op_point.velocity)
        q_bar = self.op_point.q * self.airplane.c_ref / (2 * self.op_point.velocity)
        r_bar = self.op_point.r * self.airplane.b_ref / (2 * self.op_point.velocity)

        run_file_contents += [
            f"a a {float(self.op_point.alpha)}",
            f"b b {float(self.op_point.beta)}",
            f"r r {float(p_bar)}",
            f"p p {float(q_bar)}",
            f"y y {float(r_bar)}",
        ]

        # ADDED: Blowing Modulation

        if blowing:
            for jet_name, magnitude in blowing.items():
                if jet_name == "Tcp":
                    run_file_contents += [
                        f"{'J1'} {'JT'} {float(magnitude)}",  # TODO: temporary general blowing parameter
                    ]
                else:
                    run_file_contents += [
                        f"{jet_name} {jet_name} {float(magnitude)}",  # Set index and value
                    ]

        if flap_deflections:
            for flap_name, magnitude in flap_deflections.items():
                run_file_contents += [
                    f"{flap_name} {flap_name} {float(magnitude)}",  # Set index and value
                ]
        else:
            # Set control surface deflections
            run_file_contents += ["d1 d1 1"]

        # ADDED: Trim functionality
        if trim_Cm_to_zero:
            run_file_contents += [
                f"{trim_variable} pm 0",
            ]

        return run_file_contents


class JetParam:
    def __init__(
        self,
        name="JET",
        hdisk=0.45,
        fh=1.0,
        djet0=-2.0,
        djet1=-0.2,
        djet3=-0.0003,
        dxdsk=0.0,
        dndsk=-0.0,
    ):
        """
        Jet Param class

        :param hdisk: Disk height scaling factor
        :param fh: Type of propulsor (0 for long-duct, 1 for no-duct)
        :param djet0: Jet deviation angle due to TE wedge angle
        :param djet1: Incomplete jet turning due to finite jet height/flap ratio
        :param djet3: BL separation effect on jet turning
        """
        self.name = name
        self.hdisk = hdisk
        self.fh = fh
        self.djet0 = djet0
        self.djet1 = djet1
        self.djet3 = djet3
        self.dxdsk = dxdsk
        self.dndsk = dndsk


class JetControl:
    def __init__(self, jet_name, gain, sgn_dup):
        """
        Jet Control class

        :param jet_name: Name of the jet control variable
        :param gain: Control gain (Delta_Vjet/Vinf per jet variable)
        :param sgn_dup: Symmetric (1) or differential (-1) blowing
        :param hdisk: Disk height scaling factor
        :param fh: Type of propulsor (0 for long-duct, 1 for no-duct)
        :param djet0: Jet deviation angle due to TE wedge angle
        :param djet1: Incomplete jet turning due to finite jet height/flap ratio
        :param djet3: BL separation effect on jet turning
        """
        self.jet_name = jet_name
        self.gain = gain
        self.sgn_dup = sgn_dup


class WingJSec(asb.WingXSec):
    def __init__(
        self, xyz_le, chord, twist, airfoil, control_surfaces=None, JetControls=[]
    ):
        """
        Extended WingXSec class with Jet Control (JETCONTROL & JETPARAM)

        :param xyz_le: Leading edge position [x, y, z]
        :param chord: Chord length at this section
        :param twist: Twist angle at this section
        :param airfoil: AeroSandbox Airfoil object
        :param control_surfaces: List of control surfaces
        :param JetControls: List of JetControl objects
        """
        super().__init__(
            xyz_le=xyz_le,
            chord=chord,
            twist=twist,
            airfoil=airfoil,
            control_surfaces=control_surfaces,
        )
        self.JetControls = JetControls
