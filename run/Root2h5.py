import ROOT
import numpy as np
import matplotlib.pyplot as plt
import h5py


path = 'readDigi.root'

rfile = ROOT.TFile(path);
rtree = rfile.Get('DigiHits')
print(rtree.GetEntries())

Nevt = rtree.GetEntries()
NhitMax = 512
data_dim = [Nevt]+[NhitMax,5] #(x, y, z, E, T) for at most 512 hits
htag_dim = [Nevt]+[NhitMax,1] #tag for each hit
label_dim = [Nevt]+[1] # event pid

data = np.zeros(data_dim)
htag = np.zeros(htag_dim)
pid = np.zeros(label_dim)

# Read from Root file: 
for evt in range(0, Nevt):
  rtree.GetEntry(evt)
  EcalBHit_x = np.array(getattr(rtree, 'EcalBHit_x'))
  EcalBHit_y = np.array(getattr(rtree, 'EcalBHit_y'))
  EcalBHit_z = np.array(getattr(rtree, 'EcalBHit_z'))
  EcalBHit_E = np.array(getattr(rtree, 'EcalBHit_E'))
  EcalBHit_T = np.array(getattr(rtree, 'EcalBHit_T'))
  EcalBHit_pid = np.array(getattr(rtree, 'EcalBHit_MCpid'))
  EcalBHit_Htag = np.array(getattr(rtree, 'EcalBHit_MCtag'))

  HcalBHit_x = np.array(getattr(rtree, 'HcalBHit_x'))
  HcalBHit_y = np.array(getattr(rtree, 'HcalBHit_y'))
  HcalBHit_z = np.array(getattr(rtree, 'HcalBHit_z'))
  HcalBHit_E = np.array(getattr(rtree, 'HcalBHit_E'))
  HcalBHit_T = np.array(getattr(rtree, 'HcalBHit_T'))
  HcalBHit_pid = np.array(getattr(rtree, 'HcalBHit_MCpid'))
  HcalBHit_Htag = np.array(getattr(rtree, 'HcalBHit_MCtag'))

  if (len(EcalBHit_x)+len(HcalBHit_x)>NhitMax): continue 
  if (len(EcalBHit_x)==0 and len(HcalBHit_x)==0): continue

  data[evt][:len(EcalBHit_x), 0] = EcalBHit_x
  data[evt][:len(EcalBHit_x), 1] = EcalBHit_y
  data[evt][:len(EcalBHit_x), 2] = EcalBHit_z
  data[evt][:len(EcalBHit_x), 3] = EcalBHit_E
  data[evt][:len(EcalBHit_x), 4] = EcalBHit_T
  htag[evt][:len(EcalBHit_x), 0] = EcalBHit_Htag


  data[evt][len(EcalBHit_x):len(EcalBHit_x)+len(HcalBHit_x), 0] = HcalBHit_x
  data[evt][len(EcalBHit_x):len(EcalBHit_x)+len(HcalBHit_x), 1] = HcalBHit_y
  data[evt][len(EcalBHit_x):len(EcalBHit_x)+len(HcalBHit_x), 2] = HcalBHit_z
  data[evt][len(EcalBHit_x):len(EcalBHit_x)+len(HcalBHit_x), 3] = HcalBHit_E
  data[evt][len(EcalBHit_x):len(EcalBHit_x)+len(HcalBHit_x), 4] = HcalBHit_T
  htag[evt][len(EcalBHit_x):len(EcalBHit_x)+len(HcalBHit_x), 0] = HcalBHit_Htag
  pid[evt] = HcalBHit_pid[0]



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


