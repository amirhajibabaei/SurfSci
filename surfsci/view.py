# +
import nglview as ngl
from ase import Atoms


def view(traj, radiusScale=1.0):
    if type(traj) == Atoms:
        traj = [traj]
    view = ngl.show_asetraj(traj)
    view.add_unitcell()
    view.add_spacefill()
    view.remove_ball_and_stick()
    view.camera = "orthographic"
    view.parameters = {"clipDist": 0}
    view.update_spacefill(
        radiusType="covalent", radiusScale=radiusScale, color_scale="rainbow"
    )
    return view
