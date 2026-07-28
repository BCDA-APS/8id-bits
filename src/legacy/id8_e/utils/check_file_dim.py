import h5py

def check_h5_shape(path):

    f = h5py.File(path)

    data_dimensions = f['entry']['data']['data'].shape

    return data_dimensions