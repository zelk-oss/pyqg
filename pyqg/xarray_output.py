import numpy as np
try:
    import xarray as xr
except ImportError:
    raise ImportError(
        "Xarray output in Pyqg requires the Xarray package, which is not installed on your system. " 
        "Please install Xarray in order to activate this feature. "
        "Instructions at http://xarray.pydata.org/en/stable/getting-started-guide/installing.html#instructions"
    )
    
from pyqg.errors import DiagnosticNotFilledError

# Define dict for variable dimensions
spatial_dims = ('time','lev','y','x')
time_dims = ('time')
dim_database = {
    'q': spatial_dims,
    'u': spatial_dims,
    'v': spatial_dims,
    'times': time_dims
}

# dict for variable dimensions
var_attr_database = {
    'q':     { 'units': 's^-1',      'long_name': 'potential vorticity in real space',},
    'u':     { 'units': 'm s^-1',    'long_name': 'zonal velocity anomaly',},
    'v':     { 'units': 'm s^-1',    'long_name': 'meridional velocity anomaly',},
    'times':     { 'units': 's',    'long_name': 'times',}
}


# dict for coordinate dimensions
coord_database = {
    'time': ('time'),
    'lev': ('lev'),
    'lev_mid': ('lev_mid'),
    'x': ('x'),
    'y': ('y')
}

# dict for coordinate attributes 
coord_attr_database = {
    'time': {'long_name': 'model time', 'units': 's',},
    'lev': {'long_name': 'vertical levels',},
    'lev_mid': {'long_name': 'vertical level interface',},
    'x': {'long_name': 'real space grid points in the x direction', 'units': 'grid point',},
    'y': {'long_name': 'real space grid points in the y direction', 'units': 'grid point',}
}

# list for dataset attributes
attribute_database = [
    'beta',
    'delta',
    'del2',
    'dt',
    'filterfac',
    'L',
    'M',
    'ntd',
    'nx',
    'ny',
    'nz',
    'rd',
    'rho',
    'rek',
    'nek',
    'cphi',
    'taveint',
    'tavestart',
    'tc',
    'tmax',
    'tsnapint',
    'tsnapstart',
    'twrite',
    'W',
]

# Transform certain key coordinates
transformations = {
    'time': lambda x: np.array([x.t]),
    'lev': lambda x: np.arange(1,x.nz+1),
    'lev_mid': lambda x: np.arange(1.5,x.nz+.5),
    'x': lambda x: x.x[0,:],
    'y': lambda x: x.y[:,0]
}



#######################
def model_to_dataset(m):
    '''Convert outputs from model to an xarray dataset'''

    # Create a dictionary of variables
    variables = {}
    for vname in dim_database:
        if hasattr(m,vname):
            data = getattr(m,vname, None).copy()
            if 'time' in dim_database[vname]:
                variables[vname] = (dim_database[vname], data[np.newaxis,...], var_attr_database[vname])
            else:
                variables[vname] = (dim_database[vname], data, var_attr_database[vname])

    # Create a dictionary of coordinates
    coordinates = {}
    for cname in coord_database:
        if(cname=='time'):
            data = np.array([m.t]).copy()
        else:
            data = transformations[cname](m).copy()
        coordinates[cname] = (coord_database[cname], data, coord_attr_database[cname])

    # Create a dictionary of global attributes
    global_attrs = {}
    for aname in attribute_database:
        if hasattr(m, aname):
            data = getattr(m, aname)
            global_attrs[f"pyqg:{aname}"] = (data)
        
    diagnostics = {}
    
    
    ds = xr.Dataset(variables, coords=coordinates, attrs=global_attrs)
    ds.attrs['title'] = 'pyqg: Python Quasigeostrophic Model'
    ds.attrs['reference'] = 'https://pyqg.readthedocs.io/en/latest/index.html'
    
    return ds





################################################################
def create_dataset_gl(m, nt_total=None):
    '''Convert outputs from model to an xarray dataset'''

    if(nt_total==None):
       nt_total=1

    # Create a dictionary of coordinates
    coordinates = {}
    for cname in coord_database:
        if(cname=='time'):
            data = np.arange(nt_total+0.0) # to create a real
            coordinates[cname] = (coord_database[cname], data, coord_attr_database[cname])
        else:
            data = transformations[cname](m).copy()
            coordinates[cname] = (coord_database[cname], data, coord_attr_database[cname])


# Create a dictionary of variables
    variables = {}
    for vname in dim_database:
        if(vname=='times'):
                 data = np.arange(nt_total+0.0) # to create a real
                 variables[vname] = (dim_database[vname], data, var_attr_database[vname])
        else:
                 if hasattr(m,vname):
                        if 'time' in dim_database[vname]:
                            #data = getattr(m,vname, None).copy()
                            #print(vname,np.shape(data))
                            #variables[vname] = (dim_database[vname], data[np.newaxis,...], var_attr_database[vname])
                            data = np.zeros((
                            coordinates['time'][1].shape[0], coordinates['lev'][1].shape[0],
                                coordinates['y'][1].shape[0], coordinates['x'][1].shape[0]))
                            variables[vname] = (dim_database[vname], data, var_attr_database[vname])
                        else:
                            print('not here',vname)
                            pass
                            data = getattr(m,vname, None).copy()
                            variables[vname] = (dim_database[vname], data, var_attr_database[vname])

    # Create a dictionary of global attributes
    global_attrs = {}
    for aname in attribute_database:
        if hasattr(m, aname):
            data = getattr(m, aname)
            global_attrs[f"pyqg:{aname}"] = (data)
        
    
    ds = xr.Dataset(variables, coords=coordinates, attrs=global_attrs)
    ds.attrs['title'] = 'pyqg: Python Quasigeostrophic Model'
    ds.attrs['reference'] = 'https://pyqg.readthedocs.io/en/latest/index.html'
    
    return ds


    
def model_to_dataset_gl(ds, m, itime=None):
    '''Convert outputs from model to an xarray dataset'''

    if(itime==None):
        itime=0

    variables = {}
    for vname in dim_database:
        if (vname=='times'):
             ds[vname][itime]=m.t
             ds['time'].data[itime]=m.t
        else:
             ds[vname][itime,:,:,:] =getattr(m, vname, None).copy()


    return ds
