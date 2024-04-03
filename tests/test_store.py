import os

import MDAnalysis as mda
import xarray as xr


def test_store_rogfp2(gcs_fs, gcs_uuid, gcs_cache_dir):

    test_dir = os.path.join(gcs_cache_dir, "rogfp2-relax")
    os.makedirs(test_dir, exist_ok=True)
    file_topo = os.path.join(test_dir, "mol.prmtop")
    if not os.path.exists(file_topo):
        uri_topo = os.path.join(gcs_uuid, "rogfp2-relax/mol.prmtop")
        with gcs_fs.open_input_stream(uri_topo) as fs_topo:
            with open(file_topo, "w", encoding="utf-8") as f:
                f.write(fs_topo.readall())

    file_coord = os.path.join(test_dir, "05_relax_nvt_r.nc")
    if not os.path.exists(file_coord):
        uri_coord = os.path.join(gcs_uuid, "rogfp2-relax/05_relax_nvt_r.nc")
        with gcs_fs.open_input_stream(uri_coord) as fs_coord:
            dataset = xr.open_dataset(fs_coord, engine="scipy")
            dataset.to_netcdf(file_coord)

    u = mda.Universe(file_topo, file_coord)
    assert u.atoms.positions.shape == (37381, 3)
