from tempfile import NamedTemporaryFile

import MDAnalysis as mda
import xarray as xr


def test_store_rogfp2(gcs_fs, uuid_rogfp2):
    uri_topo = uuid_rogfp2 + "/topo/mol.prmtop"
    uri_coord = uuid_rogfp2 + "/outputs/06_relax_npt_r.nc"
    file_topo = NamedTemporaryFile(suffix=".prmtop")
    file_coord = NamedTemporaryFile(suffix=".nc")
    with gcs_fs.open_input_stream(uri_topo) as fs_topo:
        file_topo.write(fs_topo.readall())
    with gcs_fs.open_input_file(uri_coord) as fs_coord:
        dataset = xr.open_dataset(fs_coord, engine="scipy")
        dataset.to_netcdf(file_coord.name)
    u = mda.Universe(file_topo.name, file_coord.name)
    assert u.atoms.positions.shape == (37381, 3)
