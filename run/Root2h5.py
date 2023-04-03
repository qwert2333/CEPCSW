import ROOT
import numpy as np
import matplotlib.pyplot as plt
import h5py


path = 'readDigi.root'

rfile = ROOT.TFile(path);
rtree = rfile.Get('DigiHits')
print(rtree.GetEntries())

Nevt = rtree.GetEntries()
data_dim = [Nevt]+[512,5] #(x, y, z, E, T) for at most 512 hits
htag_dim = [Nevt]+[512,1] #tag for each hit
label_dim = [Nevt]+[1] # event pid

data = np.zeros(data_dim)
htag = np.zeros(htag_dim)
pid = np.zeros(label_dim)

# Read from Root file: 
for evt in range(0, Nevt):
  rtree.GetEntry(evt)
  EcalEHit_x = np.array(getattr(rtree, 'EcalBHit_x'))
  EcalEHit_y = np.array(getattr(rtree, 'EcalBHit_y'))
  EcalEHit_z = np.array(getattr(rtree, 'EcalBHit_z'))
  EcalEHit_E = np.array(getattr(rtree, 'EcalBHit_E'))
  EcalEHit_T = np.array(getattr(rtree, 'EcalBHit_T'))
  EcalEHit_pid = np.array(getattr(rtree, 'EcalBHit_MCpid'))
  EcalEHit_Htag = np.array(getattr(rtree, 'EcalBHit_MCtag'))

  data[evt][:len(EcalEHit_x), 0] = EcalEHit_x
  data[evt][:len(EcalEHit_x), 1] = EcalEHit_y
  data[evt][:len(EcalEHit_x), 2] = EcalEHit_z
  data[evt][:len(EcalEHit_x), 3] = EcalEHit_E
  data[evt][:len(EcalEHit_x), 4] = EcalEHit_T
  htag[evt][:len(EcalEHit_x), 0] = EcalEHit_Htag
  pid[evt] = EcalEHit_pid[0]


#Write to hdf5:
hf = h5py.File('test.h5','w')
hf.create_dataset('data', data=data)
hf.create_dataset('htag', data=htag)
hf.create_dataset('pid', data=pid)


'''
wfile = h5py.File('dataset.h5', 'w')
wfile.create_dataset('HitCloud_x', data = EcalEHit_x)
wfile.create_dataset('HitCloud_y', data = EcalEHit_y)
wfile.create_dataset('HitCloud_z', data = EcalEHit_z)
wfile.create_dataset('HitCloud_E', data = EcalEHit_E)
wfile.create_dataset('HitCloud_T', data = EcalEHit_T)
wfile.create_dataset('HitCloud_pid', data = EcalEHit_pid)
'''


