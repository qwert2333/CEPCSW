import numpy as np
import h5py

f = h5py.File('test.h5','r')
print( f.keys() )
keys = list(f.keys())[0]


print( 'data size: ', len(f[list(f.keys())[0]]) )
print( '%s shape: ' % list(f.keys())[0] , f[ list(f.keys())[0] ][0].shape)
print( '%s shape: ' % list(f.keys())[1] , f[ list(f.keys())[1] ][0].shape)
print( '%s shape: ' % list(f.keys())[2] , f[ list(f.keys())[2] ][0].shape)


