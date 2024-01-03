import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import pyvista as pv
import sys
sys.path.append('/workfs2/bes/zyang/code-store/yang/python_functions')
import yyy_draw_options as yyy
import ECALParams as EP

#######################################
def Morandi_color_palette(i):
  color_list = ['#D89C7A', '#849B91', '#686789', '#B0B1B6', '#979771', 
                '#91AD9E', '#B77F70', '#BEB1A8', '#A79A89', '#8A95A9', 
                '#9AA690', '#99857E', '#7D7465', '#88878D', '#B4746B', 
                '#676662', '#AB545A', '#9A7549', '#B57C82', "#C9842E"
                ]
  return color_list[i%len(color_list)]

def color_list(i):
  color_list = ['navy', 'firebrick', 'darkgreen', 'darkorange', 'blue', 'cyan', 'magenta', 'chocolate']
  return color_list[i%len(color_list)]

def point_in_Tracker(x, y):
  xx = np.array([])
  yy = np.array([])
  for i in range(len(x)):
    if x[i]>1860 or x[i]<-1860:
      continue
    if y[i]>1860 or y[i]<-1860:
      continue
    if y[i]>-x[i]+1860*(2**0.5) or y[i]>x[i]+1860*(2**0.5):
      continue
    if y[i]<-x[i]-1860*(2**0.5) or y[i]<x[i]-1860*(2**0.5):
      continue
    xx = np.append(xx, x[i])
    yy = np.append(yy, y[i])
  return xx, yy


#######################################
# MC
def draw_MC_Track(rec_path, plotter, event, pTcut=0):
  recfile = uproot.open(rec_path)
  mcStatus = recfile['RecClusters']['mcStatus'].array()
  mcPx = recfile['RecClusters']['mcPx'].array()
  mcPy = recfile['RecClusters']['mcPy'].array()
  mcPz = recfile['RecClusters']['mcPz'].array()
  mcMass = recfile['RecClusters']['mcMass'].array()
  mcCharge = recfile['RecClusters']['mcCharge'].array()
  
  evt_status = mcStatus[event]
  evt_px = mcPx[event]
  evt_py = mcPy[event]
  evt_pz = mcPz[event]
  evt_mass = mcMass[event]
  evt_charge = mcCharge[event]

  px = evt_px[(evt_status==1) & (evt_charge!=0)]
  py = evt_py[(evt_status==1) & (evt_charge!=0)]
  pz = evt_pz[(evt_status==1) & (evt_charge!=0)]
  mass = evt_mass[(evt_status==1) & (evt_charge!=0)]
  charge = evt_charge[(evt_status==1) & (evt_charge!=0)]

  pt = (px**2 + py**2)**0.5
  rho = 10000/9*pt  # unit: mm
  phi0 = np.arctan2(py, px)  # initial momentum direction
  phi1 = phi0 - charge*np.pi/2
  # center point of the circle
  cx = rho * np.cos(phi1)
  cy = rho * np.sin(phi1)
  for i in range(len(pt)):
    if pt[i]<pTcut:
      continue
    # if rho[i] > 1860/2:
    dphi = np.linspace(0, np.pi, 1000)
    xx = ( cx[i] + rho[i] * np.cos(phi1[i] - np.pi + dphi) ) if charge[i]<0 else ( cx[i] + rho[i] * np.cos(phi1[i] - np.pi - dphi) )
    yy = ( cy[i] + rho[i] * np.sin(phi1[i] - np.pi + dphi) ) if charge[i]<0 else ( cy[i] + rho[i] * np.sin(phi1[i] - np.pi - dphi) )

    f_x, f_y = point_in_Tracker(xx, yy)
    f_z = (f_x**2 + f_y**2)**0.5 * pz[i] / pt[i]
    line = pv.MultipleLines(np.column_stack((f_x, f_y, f_z)))
    plotter.add_mesh(line, color='deepskyblue', opacity=0.2, line_width=10)

def draw_MC_photon(rec_path, plotter, event, Ecut=0):
  recfile = uproot.open(rec_path)
  mcStatus = recfile['RecClusters']['mcStatus'].array()
  mcPdgid = recfile['RecClusters']['mcPdgid'].array()
  mcPx = recfile['RecClusters']['mcPx'].array()
  mcPy = recfile['RecClusters']['mcPy'].array()
  mcPz = recfile['RecClusters']['mcPz'].array()
  
  evt_status = np.array(mcStatus[event])
  evt_pdg = np.array(mcPdgid[event])
  px = np.array(mcPx[event])[(evt_status==1) & (evt_pdg==22)][2:]
  py = np.array(mcPy[event])[(evt_status==1) & (evt_pdg==22)][2:]
  pz = np.array(mcPz[event])[(evt_status==1) & (evt_pdg==22)][2:]
  p = (px**2+py**2+pz**2)**0.5
  x = px[p>Ecut]/p[p>Ecut]*3300
  y = py[p>Ecut]/p[p>Ecut]*3300
  z = pz[p>Ecut]/p[p>Ecut]*3300
  for i in range(len(x)):
    line = pv.Line((0, 0, 0), (x[i], y[i], z[i]))
    plotter.add_mesh(line, color='orangered', opacity=0.7, line_width=2)

def draw_MC_h0(rec_path, plotter, event, Ecut=0):
  recfile = uproot.open(rec_path)
  mcStatus = recfile['RecClusters']['mcStatus'].array()
  mcPdgid = recfile['RecClusters']['mcPdgid'].array()
  mcCharge = recfile['RecClusters']['mcCharge'].array()
  mcMass = recfile['RecClusters']['mcMass'].array()
  mcPx = recfile['RecClusters']['mcPx'].array()
  mcPy = recfile['RecClusters']['mcPy'].array()
  mcPz = recfile['RecClusters']['mcPz'].array()
  
  evt_status = np.array(mcStatus[event])
  evt_pdg = np.array(mcPdgid[event])
  evt_charge = np.array(mcCharge[event])
  evt_mass = np.array(mcMass[event])
  evt_px = np.array(mcPx[event])
  evt_py = np.array(mcPy[event])
  evt_pz = np.array(mcPz[event])

  pdg = evt_pdg[(evt_status==1) & ((evt_pdg==130) | (evt_pdg==2112) | (evt_pdg==-2112))]
  result = np.any((pdg != 130) & (pdg != 2112) & (pdg != -2112))
  if result:
    print("Error!!!")

  px = evt_px[(evt_status==1) & ((evt_pdg==130) | (evt_pdg==2112) | (evt_pdg==-2112))]
  py = evt_py[(evt_status==1) & ((evt_pdg==130) | (evt_pdg==2112) | (evt_pdg==-2112))]
  pz = evt_pz[(evt_status==1) & ((evt_pdg==130) | (evt_pdg==2112) | (evt_pdg==-2112))]
  p = (px**2+py**2+pz**2)**0.5
  x = px[p>Ecut]/p[p>Ecut]*3300
  y = py[p>Ecut]/p[p>Ecut]*3300
  z = pz[p>Ecut]/p[p>Ecut]*3300
  for i in range(len(x)):
    line = pv.Line((0, 0, 0), (x[i], y[i], z[i]))
    plotter.add_mesh(line, color='limegreen', opacity=0.7, line_width=2)

def draw_endpoint(sim_path, plotter, event):
  simfile = uproot.open(sim_path)
  status = simfile['events']['MCParticle.generatorStatus'].array()
  end_x = simfile['events']['MCParticle.endpoint.x'].array()
  end_y = simfile['events']['MCParticle.endpoint.y'].array()
  end_z = simfile['events']['MCParticle.endpoint.z'].array()

  evt_status = np.array(status[event])
  evt_end_x = np.array(end_x[event])[(evt_status==1)]
  evt_end_y = np.array(end_y[event])[(evt_status==1)]
  evt_end_z = np.array(end_z[event])[(evt_status==1)]

  xx = evt_end_x[(np.abs(evt_end_x)<3500) & (np.abs(evt_end_y)<3500) & (np.abs(evt_end_z)<5000)]
  yy = evt_end_y[(np.abs(evt_end_x)<3500) & (np.abs(evt_end_y)<3500) & (np.abs(evt_end_z)<5000)]
  zz = evt_end_z[(np.abs(evt_end_x)<3500) & (np.abs(evt_end_y)<3500) & (np.abs(evt_end_z)<5000)]

  points = np.column_stack((xx, yy, zz))
  pts = pv.PolyData(points)
  plotter.add_mesh(pts, point_size=10, style='points', color='red')

def draw_vertex(sim_path, plotter, event):
  simfile = uproot.open(sim_path)
  status = simfile['events']['MCParticle.generatorStatus'].array()
  vtx_x = simfile['events']['MCParticle.vertex.x'].array()
  vtx_y = simfile['events']['MCParticle.vertex.y'].array()
  vtx_z = simfile['events']['MCParticle.vertex.z'].array()

  evt_status = np.array(status[event])
  evt_vtx_x = np.array(vtx_x[event])[(evt_status==1)]
  evt_vtx_y = np.array(vtx_y[event])[(evt_status==1)]
  evt_vtx_z = np.array(vtx_z[event])[(evt_status==1)]

  xx = evt_vtx_x[(np.abs(evt_vtx_x)<3500) & (np.abs(evt_vtx_y)<3500) & (np.abs(evt_vtx_z)<5000)]
  yy = evt_vtx_y[(np.abs(evt_vtx_x)<3500) & (np.abs(evt_vtx_y)<3500) & (np.abs(evt_vtx_z)<5000)]
  zz = evt_vtx_z[(np.abs(evt_vtx_x)<3500) & (np.abs(evt_vtx_y)<3500) & (np.abs(evt_vtx_z)<5000)]

  points = np.column_stack((xx, yy, zz))
  pts = pv.PolyData(points)
  plotter.add_mesh(pts, point_size=10, style='points', color='gold')
#######################################
# Tracks
def draw_rec_Track(rec_path, plotter, event):
  recfile = uproot.open(rec_path)
  m_trkstate_phi = recfile['RecTracks']['m_trkstate_phi'].array()
  m_trkstate_tanL = recfile['RecTracks']['m_trkstate_tanL'].array()
  m_trkstate_kappa = recfile['RecTracks']['m_trkstate_kappa'].array()
  m_trkstate_location = recfile['RecTracks']['m_trkstate_location'].array()
  
  evt_phi = np.array(m_trkstate_phi[event])
  evt_tanL = np.array(m_trkstate_tanL[event])
  evt_kappa = np.array(m_trkstate_kappa[event])
  evt_location = np.array(m_trkstate_location[event])

  loc1_phi = evt_phi[evt_location==1]
  loc1_tanL = evt_tanL[evt_location==1]
  loc1_kappa = evt_kappa[evt_location==1]

  charge = np.where(loc1_kappa>0, 1, -1)
  B = 3.0
  rho = 10000/3/B/np.abs(loc1_kappa)
  p_T = 9/10000*rho
  px = p_T * np.cos(loc1_phi)
  py = p_T * np.sin(loc1_phi)
  pz = p_T * loc1_tanL

  pt = (px**2 + py**2)**0.5
  rho = 10000/9*pt  # unit: mm
  phi0 = np.arctan2(py, px)  # initial momentum direction
  phi1 = phi0 - charge*np.pi/2
  # center point of the circle
  cx = rho * np.cos(phi1)
  cy = rho * np.sin(phi1)
  for i in range(len(px)):
    dphi = np.linspace(0, np.pi, 1000)
    xx = ( cx[i] + rho[i] * np.cos(phi1[i] - np.pi + dphi) ) if charge[i]<0 else ( cx[i] + rho[i] * np.cos(phi1[i] - np.pi - dphi) )
    yy = ( cy[i] + rho[i] * np.sin(phi1[i] - np.pi + dphi) ) if charge[i]<0 else ( cy[i] + rho[i] * np.sin(phi1[i] - np.pi - dphi) )

    f_x, f_y = point_in_Tracker(xx, yy)
    f_z = (f_x**2 + f_y**2)**0.5 * pz[i] / pt[i]
    line = pv.MultipleLines(np.column_stack((f_x, f_y, f_z)))
    plotter.add_mesh(line, color='dodgerblue', opacity=0.4, line_width=2)

def draw_tracker_hit(sim_parh, plotter, event):
  simfile = uproot.open(sim_parh)
  x_VXD_all = simfile['events']['VXDCollection.position.x'].array()
  y_VXD_all = simfile['events']['VXDCollection.position.y'].array()
  z_VXD_all = simfile['events']['VXDCollection.position.z'].array()
  x_SIT_all = simfile['events']['SITCollection.position.x'].array()
  y_SIT_all = simfile['events']['SITCollection.position.y'].array()
  z_SIT_all = simfile['events']['SITCollection.position.z'].array()
  x_TPC_all = simfile['events']['TPCCollection.position.x'].array()
  y_TPC_all = simfile['events']['TPCCollection.position.y'].array()
  z_TPC_all = simfile['events']['TPCCollection.position.z'].array()
  x_SET_all = simfile['events']['SETCollection.position.x'].array()
  y_SET_all = simfile['events']['SETCollection.position.y'].array()
  z_SET_all = simfile['events']['SETCollection.position.z'].array()

  x_VXD_evt = np.array(x_VXD_all[event])
  y_VXD_evt = np.array(y_VXD_all[event])
  z_VXD_evt = np.array(z_VXD_all[event])
  x_SIT_evt = np.array(x_SIT_all[event])
  y_SIT_evt = np.array(y_SIT_all[event])
  z_SIT_evt = np.array(z_SIT_all[event])
  x_TPC_evt = np.array(x_TPC_all[event])
  y_TPC_evt = np.array(y_TPC_all[event])
  z_TPC_evt = np.array(z_TPC_all[event])
  x_SET_evt = np.array(x_SET_all[event])
  y_SET_evt = np.array(y_SET_all[event])
  z_SET_evt = np.array(z_SET_all[event])

  xx = np.concatenate((x_VXD_evt, x_SIT_evt, x_TPC_evt, x_SET_evt))
  yy = np.concatenate((y_VXD_evt, y_SIT_evt, y_TPC_evt, y_SET_evt))
  zz = np.concatenate((z_VXD_evt, z_SIT_evt, z_TPC_evt, z_SET_evt))

  points = np.column_stack((xx, yy, zz))
  pts = pv.PolyData(points)
  plotter.add_mesh(pts, point_size=4, style='points', color='blue')
  

#######################################
# ECAL crystal bars
def get_module_dlayer_slayer(x, y, z):
  """
  given position of a bar,, return the module, dlayer and slayer of the bar
  """
  # norm vectors of all modules
  norm_vectors = [[0, 1], [-1/2**0.5, 1/2**0.5], [-1, 0], [-1/2**0.5, -1/2**0.5], [0, -1], [1/2**0.5, -1/2**0.5], [1, 0], [1/2**0.5, 1/2**0.5]]
  module = -1
  dlayer = -1
  slayer = -1

  if (y>1860) and (y<2140) and (y<1860*(2**0.5)-x): 
    module = 0
  elif (x>-2140) and (x<-1860) and (y<1860*(2**0.5)+x): 
    module = 2
  elif (y>-2140) and (y<-1860) and (y>-1860*(2**0.5)-x): 
    module = 4
  elif (x>1860) and (x<2140) and (y>-1860*(2**0.5)+x):  
    module = 6
  elif (x<0) and (y>0):
    module = 1
  elif (x<0) and (y<0):
    module = 3
  elif (x>0) and (y<0):
    module = 5
  elif (x>0) and (y>0):
    module = 7
  else:
    print('wrong module')
    return module, dlayer, slayer
  
  # verticle distance for (0, 0) to (x, y)
  distance = (x * norm_vectors[module][0]) + (y * norm_vectors[module][1])
  dlayer = int((distance-1860)/20) + 1
  if dlayer<1 or dlayer>14:
    print('wrong dlayer')
  # slayer = 0 if np.abs((distance-1865)%20)<1 else 1
  if np.abs((round(distance)-1865)%20)<1:
    slayer = 0
  elif np.abs((round(distance)-1865)%20 - 10)<1:
    slayer = 1
  else:
    print('wrong slayer', distance, (distance-1865)/20)
  return module, dlayer, slayer

def draw_bar(plotter, x, y, z, module, dlayer, slayer, color):
  """
  draw a bar using the given infomation
  """
  bar_width = 10
  # parallel to z axis
  if slayer==1:  
    bar_length = 600
    if module in [0, 2, 4, 6]:
      cube = pv.Cube(center=[x, y, z], 
                      x_length=bar_width, y_length=bar_width, z_length=bar_length)
      plotter.add_mesh(cube, color=color, opacity=0.2)
    elif module in [1, 3, 5, 7]:
      cube = pv.Cube(center=[x, y, z], 
                      x_length=bar_width, y_length=bar_width, z_length=bar_length)
      cube.rotate_z(45, (x, y, z), inplace=True)
      plotter.add_mesh(cube, color=color, opacity=0.2)
    else:
      print('..wrong module..')    
    
  # perpendicular to z axis
  elif slayer==0:  
    bar_length = 480 - dlayer*10
    rotate = 0
    
    if module in [0, 4]:
      rotate = 0
    elif module in [2, 6]:
      rotate = 90
    elif module in [1, 5]:
      rotate = 45
    elif module in [3, 7]:
      rotate = -45

    cube = pv.Cube(center=[x, y, z], 
                    x_length=bar_length, y_length=bar_width, z_length=bar_width)
    cube.rotate_z(rotate, (x, y, z), inplace=True)
    plotter.add_mesh(cube, color=color, opacity=0.2)
  else:
    print('..wrong slayer..')

def draw_hitbars(rec_path, plotter, event, Ecut=0.001, barCol='UV'):
  recfile = uproot.open(rec_path)
  simBar_x = recfile['SimBarHit']['simBar_x'].array()
  simBar_y = recfile['SimBarHit']['simBar_y'].array()
  simBar_z = recfile['SimBarHit']['simBar_z'].array()
  simBar_Q1 = recfile['SimBarHit']['simBar_Q1'].array()
  simBar_slayer = recfile['SimBarHit']['simBar_slayer'].array()
  bar_x = np.array(simBar_x[event])
  bar_y = np.array(simBar_y[event])
  bar_z = np.array(simBar_z[event])
  bar_E = np.array(simBar_Q1[event])
  bar_slayer = np.array(simBar_slayer[event])
  barV_x = bar_x[bar_slayer==1]
  barV_y = bar_y[bar_slayer==1]
  barV_z = bar_z[bar_slayer==1]
  barV_E = bar_E[bar_slayer==1]
  barU_x = bar_x[bar_slayer==0]
  barU_y = bar_y[bar_slayer==0]
  barU_z = bar_z[bar_slayer==0]
  barU_E = bar_E[bar_slayer==0]

  # color map set with enregy information
  cmap = cm.get_cmap('YlOrRd')
  norm = colors.LogNorm(vmin=max(min(bar_E), Ecut), vmax=max(bar_E))

  if ('U' in barCol) :
    for i in range(len(barU_x)):
      if barU_E[i]<Ecut:
        continue
      module, dlayer, slayer = get_module_dlayer_slayer(barU_x[i], barU_y[i], barU_z[i])
      # draw_bar(plotter, barU_x[i], barU_y[i], barU_z[i], module, dlayer, 0, cmap(norm(barU_E[i])))
      draw_bar(plotter, barU_x[i], barU_y[i], barU_z[i], module, dlayer, 0, cmap(norm(barU_E[i])))
  if ('V' in barCol):
    for i in range(len(barV_x)):
      if barV_E[i]<Ecut:
        continue
      module, dlayer, slayer = get_module_dlayer_slayer(barV_x[i], barV_y[i], barV_z[i])
      # draw_bar(plotter, barV_x[i], barV_y[i], barV_z[i], module, dlayer, 1, cmap(norm(barV_E[i])))
      draw_bar(plotter, barV_x[i], barV_y[i], barV_z[i], module, dlayer, 1, cmap(norm(barV_E[i])))

def draw_lm(rec_path, plotter, event, barCol='UV'):
  recfile = uproot.open(rec_path)
  barShowerV_x = recfile['RecLayers']['barShowerV_x'].array()
  barShowerV_y = recfile['RecLayers']['barShowerV_y'].array()
  barShowerV_z = recfile['RecLayers']['barShowerV_z'].array()
  barShowerV_E = recfile['RecLayers']['barShowerV_E'].array()
  barShowerU_x = recfile['RecLayers']['barShowerU_x'].array()
  barShowerU_y = recfile['RecLayers']['barShowerU_y'].array()
  barShowerU_z = recfile['RecLayers']['barShowerU_z'].array()
  barShowerU_E = recfile['RecLayers']['barShowerU_E'].array()
  barV_x = np.array(barShowerV_x[event])
  barV_y = np.array(barShowerV_y[event])
  barV_z = np.array(barShowerV_z[event])
  barV_E = np.array(barShowerV_E[event])
  barU_x = np.array(barShowerU_x[event])
  barU_y = np.array(barShowerU_y[event])
  barU_z = np.array(barShowerU_z[event])
  barU_E = np.array(barShowerU_E[event])

  # color map set with enregy information
  cmap = cm.get_cmap('YlOrRd')
  bar_E = np.append(barV_E, barU_E)
  norm = colors.LogNorm(vmin=min(bar_E), vmax=max(bar_E))

  if 'U' in barCol:
    for i in range(len(barU_x)):
      module, dlayer, slayer = get_module_dlayer_slayer(barU_x[i], barU_y[i], barU_z[i])
      draw_bar(plotter, barU_x[i], barU_y[i], barU_z[i], module, dlayer, 0, cmap(norm(barU_E[i])))
  if 'V' in barCol:
    for i in range(len(barV_x)):
      module, dlayer, slayer = get_module_dlayer_slayer(barV_x[i], barV_y[i], barV_z[i])
      draw_bar(plotter, barV_x[i], barV_y[i], barV_z[i], module, dlayer, 1, cmap(norm(barV_E[i])))

def draw_ECAL_shower_center(rec_path, plotter, event):
  recfile = uproot.open(rec_path)
  Clus_x = recfile['RecClusters']['Clus_x'].array()
  Clus_y = recfile['RecClusters']['Clus_y'].array()
  Clus_z = recfile['RecClusters']['Clus_z'].array()
  Clus_E = recfile['RecClusters']['Clus_E'].array()
  evt_clus_x = np.array(Clus_x[event])
  evt_clus_y = np.array(Clus_y[event])
  evt_clus_z = np.array(Clus_z[event])
  evt_clus_E = np.array(Clus_E[event])

  for i in range(len(evt_clus_E)):
    sphere = pv.Sphere(center=[evt_clus_x[i], evt_clus_y[i], evt_clus_z[i]], radius=50)
    plotter.add_mesh(sphere, color='red', opacity=0.2)

#######################################
# HCAL cell
def get_HCAL_cell_angle(x, y):
  """
  Given the (x, y) of a cell in HCAL, return its rotate angle as a cube.
  The default cube is set in the upper module
  """
  phi = np.arctan2(y, x)
  angle = 0
  if (phi>-np.pi/8 and phi<np.pi/8) or (phi>np.pi*7/8) or (phi<-np.pi*7/8):
    # The lefy and right module
    angle = 90
  elif (phi>np.pi*3/8 and phi<np.pi*5/8) or (phi<-np.pi*3/8 and phi>-np.pi*5/8):
    # The upper and bottom module
    angle = 0
  elif (phi>np.pi*5/8 and phi<np.pi*7/8) or (phi<-np.pi/8 and phi>-np.pi*3/8):
    # The upper left and bottom right module
    angle = 45
  elif (phi>np.pi/8 and phi<np.pi*3/8) or (phi<-np.pi*5/8 and phi>-np.pi*7/8):
    # The upper right and bottom left module
    angle = -45
  else:
    print("ERROR: There are other cell in HCAL, phi =", phi/np.pi*180)
  return angle

def draw_cell(plotter, x, y, z, color):
  """
  draw a cell in HCAL using the given infomation
  """
  cell_size = 40
  cell_thickness = 6.5
  rotate_angle = get_HCAL_cell_angle(x, y)
  cube = pv.Cube(center=[x, y, z], 
                  x_length=cell_size, y_length=cell_thickness, z_length=cell_size)
  cube.rotate_z(rotate_angle, (x, y, z), inplace=True)
  plotter.add_mesh(cube, color=color, opacity=0.3)

def draw_HCAL_cells(rec_path, plotter, event):
  recfile = uproot.open(rec_path)
  HcalHit_x = recfile['SimBarHit']['HcalHit_x'].array()
  HcalHit_y = recfile['SimBarHit']['HcalHit_y'].array()
  HcalHit_z = recfile['SimBarHit']['HcalHit_z'].array()
  HcalHit_E = recfile['SimBarHit']['HcalHit_E'].array()
  cell_x = np.array(HcalHit_x[event])
  cell_y = np.array(HcalHit_y[event])
  cell_z = np.array(HcalHit_z[event])
  cell_E = np.array(HcalHit_E[event])

  # color map set with enregy information
  cmap = cm.get_cmap('YlOrRd')
  norm = colors.LogNorm(vmin=min(cell_E), vmax=max(cell_E))

  for i in range(len(cell_x)):
    draw_cell(plotter, cell_x[i], cell_y[i], cell_z[i], cmap(norm(cell_E[i])))

def draw_HCAL_3Dcluster(rec_path, plotter, event):
  recfile = uproot.open(rec_path)
  Hcal_hit_tag = recfile['RecHcalClusters']['Hcal_hit_tag'].array()
  Hcal_hit_x = recfile['RecHcalClusters']['Hcal_hit_x'].array()
  Hcal_hit_y = recfile['RecHcalClusters']['Hcal_hit_y'].array()
  Hcal_hit_z = recfile['RecHcalClusters']['Hcal_hit_z'].array()
  evt_hcal_tag = Hcal_hit_tag[event]
  evt_hcal_x = Hcal_hit_x[event]
  evt_hcal_y = Hcal_hit_y[event]
  evt_hcal_z = Hcal_hit_z[event]
  set_hcal_tag = np.array(list(set(evt_hcal_tag)))
  for i, clusterE in enumerate(set_hcal_tag):
    xx = evt_hcal_x[evt_hcal_tag==clusterE]
    yy = evt_hcal_y[evt_hcal_tag==clusterE]
    zz = evt_hcal_z[evt_hcal_tag==clusterE]
    for bb in range(len(xx)):
      draw_cell(plotter, xx[bb], yy[bb], zz[bb], color_list(i))

def draw_HCAL_shower_center(rec_path, plotter, event):
  recfile = uproot.open(rec_path)
  Clus_x = recfile['RecHcalClusters']['Hcal_clus_x'].array()
  Clus_y = recfile['RecHcalClusters']['Hcal_clus_y'].array()
  Clus_z = recfile['RecHcalClusters']['Hcal_clus_z'].array()
  Clus_E = recfile['RecHcalClusters']['Hcal_clus_E'].array()
  evt_clus_x = np.array(Clus_x[event])
  evt_clus_y = np.array(Clus_y[event])
  evt_clus_z = np.array(Clus_z[event])
  evt_clus_E = np.array(Clus_E[event])

  for i in range(len(evt_clus_E)):
    sphere = pv.Sphere(center=[evt_clus_x[i], evt_clus_y[i], evt_clus_z[i]], radius=50)
    plotter.add_mesh(sphere, color='orangered', opacity=0.1)

#######################################
# Other
def draw_extrapolate_points(rec_path, plotter, event):
  recfile = uproot.open(rec_path)
  trkstate_x = recfile['RecTracks']['m_trkstate_refx'].array()
  trkstate_y = recfile['RecTracks']['m_trkstate_refy'].array()
  trkstate_z = recfile['RecTracks']['m_trkstate_refz'].array()
  evt_x = np.array(trkstate_x[event])
  evt_y = np.array(trkstate_y[event])
  evt_z = np.array(trkstate_z[event])

  # color map set with enregy information
  # cmap = cm.get_cmap('YlOrRd')
  # # norm = colors.LogNorm(vmin=0.00001, vmax=1.01)
  # norm = colors.Normalize(vmin=0.00001, vmax=1.01)

  for i in range(len(evt_x)):
    cyld = pv.Cylinder(center=[evt_x[i], evt_y[i], evt_z[i]], direction=[1, 1, 1], radius=10, height=10)
    plotter.add_mesh(cyld, color='dodgerblue', opacity=0.3)

def draw_ECAL_profile(transparency, plotter):
  for i in range(8):
    rotate = i*45
    ECAL_inner_plane = pv.Plane(center=(1860, 0.0, 0.0), direction=(1.0, 0.0, 0.0), i_size=6600, j_size=1540)
    ECAL_inner_plane.rotate_z(rotate, (0, 0, 0), inplace=True)
    plotter.add_mesh(ECAL_inner_plane, color='gray', opacity=transparency)

    ECAL_outer_plane = pv.Plane(center=(2140, 0.0, 0.0), direction=(1.0, 0.0, 0.0), i_size=6600, j_size=1773)
    ECAL_outer_plane.rotate_z(rotate, (0, 0, 0), inplace=True)
    plotter.add_mesh(ECAL_outer_plane, color='gray', opacity=transparency)

    ECAL_bot_plane1 = pv.Triangle([ [1860,  1540/2, -6600/2], 
                                    [1860, -1540/2, -6600/2], 
                                    [2140,  1773/2, -6600/2]])
    ECAL_bot_plane1.rotate_z(rotate, (0, 0, 0), inplace=True)
    plotter.add_mesh(ECAL_bot_plane1, color='gray', opacity=transparency)

    ECAL_bot_plane2 = pv.Triangle([ [1860, -1540/2, -6600/2], 
                                    [2140, -1773/2, -6600/2], 
                                    [2140,  1773/2, -6600/2]])
    ECAL_bot_plane2.rotate_z(rotate, (0, 0, 0), inplace=True)
    plotter.add_mesh(ECAL_bot_plane2, color='gray', opacity=transparency)

    ECAL_top_plane1 = pv.Triangle([ [1860,  1540/2, 6600/2], 
                                    [1860, -1540/2, 6600/2], 
                                    [2140,  1773/2, 6600/2]])
    ECAL_top_plane1.rotate_z(rotate, (0, 0, 0), inplace=True)
    plotter.add_mesh(ECAL_top_plane1, color='gray', opacity=transparency)

    ECAL_top_plane2 = pv.Triangle([ [1860, -1540/2, 6600/2], 
                                    [2140, -1773/2, 6600/2], 
                                    [2140,  1773/2, 6600/2]])
    ECAL_top_plane2.rotate_z(rotate, (0, 0, 0), inplace=True)
    plotter.add_mesh(ECAL_top_plane2, color='gray', opacity=transparency)

    ECAL_inner_line = pv.Line((1860, 1540/2, 6600/2), (1860, 1540/2, -6600/2))
    ECAL_innerBot_line = pv.Line((1860, -1540/2, -6600/2), (1860, 1540/2, -6600/2))
    ECAL_innerTop_line = pv.Line((1860, -1540/2, 6600/2), (1860, 1540/2, 6600/2))
    ECAL_inner_line.rotate_z(rotate, (0, 0, 0), inplace=True)
    ECAL_innerBot_line.rotate_z(rotate, (0, 0, 0), inplace=True)
    ECAL_innerTop_line.rotate_z(rotate, (0, 0, 0), inplace=True)
    plotter.add_mesh(ECAL_inner_line, color='gray', opacity=5*transparency)
    plotter.add_mesh(ECAL_innerBot_line, color='gray', opacity=5*transparency)
    plotter.add_mesh(ECAL_innerTop_line, color='gray', opacity=5*transparency)

    ECAL_outer_line = pv.Line((2140, 1773/2, 6600/2), (2140, 1773/2, -6600/2))
    ECAL_outerBot_line = pv.Line((2140, -1773/2, -6600/2), (2140, 1773/2, -6600/2))
    ECAL_outerTop_line = pv.Line((2140, -1773/2, 6600/2), (2140, 1773/2, 6600/2))
    ECAL_outerBot_line.rotate_z(rotate, (0, 0, 0), inplace=True)
    ECAL_outer_line.rotate_z(rotate, (0, 0, 0), inplace=True)
    ECAL_outerTop_line.rotate_z(rotate, (0, 0, 0), inplace=True)
    plotter.add_mesh(ECAL_outer_line, color='gray', opacity=5*transparency)
    plotter.add_mesh(ECAL_outerBot_line, color='gray', opacity=5*transparency)
    plotter.add_mesh(ECAL_outerTop_line, color='gray', opacity=5*transparency)


    HCAL_inner_plane = pv.Plane(center=(2150, 0.0, 0.0), direction=(1.0, 0.0, 0.0), i_size=8960, j_size=1781)
    HCAL_inner_plane.rotate_z(rotate, (0, 0, 0), inplace=True)
    plotter.add_mesh(HCAL_inner_plane, color='darkgray', opacity=transparency)

    HCAL_outer_plane = pv.Plane(center=(3225, 0.0, 0.0), direction=(1.0, 0.0, 0.0), i_size=8960, j_size=2672)
    HCAL_outer_plane.rotate_z(rotate, (0, 0, 0), inplace=True)
    plotter.add_mesh(HCAL_outer_plane, color='darkgray', opacity=transparency)
    
    HCAL_bot_plane1 = pv.Triangle([ [2150,  1781/2, -8960/2], 
                                    [2150, -1781/2, -8960/2], 
                                    [3225,  2672/2, -8960/2]])
    HCAL_bot_plane1.rotate_z(rotate, (0, 0, 0), inplace=True)
    plotter.add_mesh(HCAL_bot_plane1, color='darkgray', opacity=transparency)

    HCAL_bot_plane2 = pv.Triangle([ [2150, -1781/2, -8960/2], 
                                    [3225, -2672/2, -8960/2], 
                                    [3225,  2672/2, -8960/2]])
    HCAL_bot_plane2.rotate_z(rotate, (0, 0, 0), inplace=True)
    plotter.add_mesh(HCAL_bot_plane2, color='darkgray', opacity=transparency)

    HCAL_top_plane1 = pv.Triangle([ [2150,  1781/2, 8960/2], 
                                    [2150, -1781/2, 8960/2], 
                                    [3225,  2672/2, 8960/2]])
    HCAL_top_plane1.rotate_z(rotate, (0, 0, 0), inplace=True)
    plotter.add_mesh(HCAL_top_plane1, color='darkgray', opacity=transparency)

    HCAL_top_plane2 = pv.Triangle([ [2150, -1781/2, 8960/2], 
                                    [3225, -2672/2, 8960/2], 
                                    [3225,  2672/2, 8960/2]])
    HCAL_top_plane2.rotate_z(rotate, (0, 0, 0), inplace=True)
    plotter.add_mesh(HCAL_top_plane2, color='darkgray', opacity=transparency)

    HCAL_inner_line = pv.Line((2150, 1781/2, 8960/2), (2150, 1781/2, -8960/2))
    HCAL_innerBot_line = pv.Line((2150, -1781/2, -8960/2), (2150, 1781/2, -8960/2))
    HCAL_innerTop_line = pv.Line((2150, -1781/2, 8960/2), (2150, 1781/2, 8960/2))
    HCAL_inner_line.rotate_z(rotate, (0, 0, 0), inplace=True)
    HCAL_innerBot_line.rotate_z(rotate, (0, 0, 0), inplace=True)
    HCAL_innerTop_line.rotate_z(rotate, (0, 0, 0), inplace=True)
    plotter.add_mesh(HCAL_inner_line, color='darkgray', opacity=5*transparency)
    plotter.add_mesh(HCAL_innerBot_line, color='darkgray', opacity=5*transparency)
    plotter.add_mesh(HCAL_innerTop_line, color='darkgray', opacity=5*transparency)
    
    HCAL_outer_line = pv.Line((3225, 2672/2, 8960/2), (3225, 2672/2, -8960/2))
    HCAL_outerBot_line = pv.Line((3225, -2672/2, -8960/2), (3225, 2672/2, -8960/2))
    HCAL_outerTop_line = pv.Line((3225, -2672/2, 8960/2), (3225, 2672/2, 8960/2))
    HCAL_outer_line.rotate_z(rotate, (0, 0, 0), inplace=True)
    HCAL_outerBot_line.rotate_z(rotate, (0, 0, 0), inplace=True)
    HCAL_outerTop_line.rotate_z(rotate, (0, 0, 0), inplace=True)
    plotter.add_mesh(HCAL_outer_line, color='darkgray', opacity=5*transparency)
    plotter.add_mesh(HCAL_outerBot_line, color='darkgray', opacity=5*transparency)
    plotter.add_mesh(HCAL_outerTop_line, color='darkgray', opacity=5*transparency)

def draw_steps(digi_ECAL_path, digi_HCAL_path, event, plotter):
  ecal_file = uproot.open(digi_ECAL_path)
  ecal_x = ecal_file['SimStep']['step_x'].array()
  ecal_y = ecal_file['SimStep']['step_y'].array()
  ecal_z = ecal_file['SimStep']['step_z'].array()
  evt_ecal_x = np.array(ecal_x[event])
  evt_ecal_y = np.array(ecal_y[event])
  evt_ecal_z = np.array(ecal_z[event])

  hcal_file = uproot.open(digi_HCAL_path)
  hcal_x = hcal_file['SimStep']['step_x'].array()
  hcal_y = hcal_file['SimStep']['step_y'].array()
  hcal_z = hcal_file['SimStep']['step_z'].array()
  evt_hcal_x = np.array(hcal_x[event])
  evt_hcal_y = np.array(hcal_y[event])
  evt_hcal_z = np.array(hcal_z[event])

  ECAL_steps = np.column_stack((evt_ecal_x, evt_ecal_y, evt_ecal_z))
  HCAL_steps = np.column_stack((evt_hcal_x, evt_hcal_y, evt_hcal_z))

  ECAL_points = pv.PolyData(ECAL_steps)
  plotter.add_mesh(ECAL_points, point_size=5, style='points', color='deepskyblue')
  HCAL_points = pv.PolyData(HCAL_steps)
  plotter.add_mesh(HCAL_points, point_size=10, style='points', color='orangered')

def draw_pfo(rec_path, plotter, event):
  recfile = uproot.open(rec_path)
  pfo_tag = recfile['PFO']['pfo_tag'].array()
  n_track = recfile['PFO']['n_track'].array()
  ecal_pfo_tag = recfile['PFO']['ecal_pfo_tag'].array()
  ecal_clus_x = recfile['PFO']['ecal_clus_x'].array()
  ecal_clus_y = recfile['PFO']['ecal_clus_y'].array()
  ecal_clus_z = recfile['PFO']['ecal_clus_z'].array()
  hcal_pfo_tag = recfile['PFO']['hcal_pfo_tag'].array()
  hcal_clus_x = recfile['PFO']['hcal_clus_x'].array()
  hcal_clus_y = recfile['PFO']['hcal_clus_y'].array()
  hcal_clus_z = recfile['PFO']['hcal_clus_z'].array()

  evt_pfo_tag = pfo_tag[event]
  evt_n_track = n_track[event]
  evt_ecal_pfo_tag = ecal_pfo_tag[event]
  evt_ecal_clus_x = ecal_clus_x[event]
  evt_ecal_clus_y = ecal_clus_y[event]
  evt_ecal_clus_z = ecal_clus_z[event]
  evt_hcal_pfo_tag = hcal_pfo_tag[event]
  evt_hcal_clus_x = hcal_clus_x[event]
  evt_hcal_clus_y = hcal_clus_y[event]
  evt_hcal_clus_z = hcal_clus_z[event]

  for it in evt_pfo_tag:
    pfo_ecal_x = evt_ecal_clus_x[evt_ecal_pfo_tag==it]
    pfo_ecal_y = evt_ecal_clus_y[evt_ecal_pfo_tag==it]
    pfo_ecal_z = evt_ecal_clus_z[evt_ecal_pfo_tag==it]
    # print("  For pfo", it, ", N ecal =", len(pfo_ecal_x))
    for ip in range(len(pfo_ecal_x)):
      if evt_n_track[it]<0.5:
        sphere = pv.Sphere(center=[pfo_ecal_x[ip], pfo_ecal_y[ip], pfo_ecal_z[ip]], radius=50)
        plotter.add_mesh(sphere, color=color_list(it), opacity=0.4)
      else:
        tetrahedron = pv.Tetrahedron(center=[pfo_ecal_x[ip], pfo_ecal_y[ip], pfo_ecal_z[ip]], radius=50)
        plotter.add_mesh(tetrahedron, color=color_list(it), opacity=0.4)

    pfo_hcal_x = evt_hcal_clus_x[evt_hcal_pfo_tag==it]
    pfo_hcal_y = evt_hcal_clus_y[evt_hcal_pfo_tag==it]
    pfo_hcal_z = evt_hcal_clus_z[evt_hcal_pfo_tag==it]
    # print("             , N hcal =", len(pfo_hcal_x))
    for ip in range(len(pfo_hcal_x)):
      if evt_n_track[it]<0.5:
        sphere = pv.Sphere(center=[pfo_hcal_x[ip], pfo_hcal_y[ip], pfo_hcal_z[ip]], radius=10)
        plotter.add_mesh(sphere, color=color_list(it), opacity=0.4)
      else:
        tetrahedron = pv.Tetrahedron(center=[pfo_hcal_x[ip], pfo_hcal_y[ip], pfo_hcal_z[ip]], radius=10)
        plotter.add_mesh(tetrahedron, color=color_list(it), opacity=0.4)

    if len(pfo_ecal_x)>0:
      for ie in range(len(pfo_ecal_x)):
        for ip in range(len(pfo_hcal_x)):
          line = pv.Line((pfo_ecal_x[ie], pfo_ecal_y[ie], pfo_ecal_z[ie]), (pfo_hcal_x[ip], pfo_hcal_y[ip], pfo_hcal_z[ip]))
          plotter.add_mesh(line, color=color_list(it), line_width=5, opacity=0.3)

def draw_pfo_direction(rec_path, plotter, event):
  recfile = uproot.open(rec_path)
  pfo_tag = recfile['PFO']['pfo_tag'].array()
  n_track = recfile['PFO']['n_track'].array()
  ecal_pfo_tag = recfile['PFO']['ecal_pfo_tag'].array()
  ecal_clus_x = recfile['PFO']['ecal_clus_x'].array()
  ecal_clus_y = recfile['PFO']['ecal_clus_y'].array()
  ecal_clus_z = recfile['PFO']['ecal_clus_z'].array()
  ecal_clus_E = recfile['PFO']['ecal_clus_E'].array()
  hcal_pfo_tag = recfile['PFO']['hcal_pfo_tag'].array()
  hcal_clus_x = recfile['PFO']['hcal_clus_x'].array()
  hcal_clus_y = recfile['PFO']['hcal_clus_y'].array()
  hcal_clus_z = recfile['PFO']['hcal_clus_z'].array()
  hcal_clus_E = recfile['PFO']['hcal_clus_E'].array()

  evt_pfo_tag = np.array(pfo_tag[event])
  evt_n_track = np.array(n_track[event])
  evt_ecal_pfo_tag = np.array(ecal_pfo_tag[event])
  evt_ecal_clus_x = np.array(ecal_clus_x[event])
  evt_ecal_clus_y = np.array(ecal_clus_y[event])
  evt_ecal_clus_z = np.array(ecal_clus_z[event])
  evt_ecal_clus_E = np.array(ecal_clus_E[event])
  evt_hcal_pfo_tag = np.array(hcal_pfo_tag[event])
  evt_hcal_clus_x = np.array(hcal_clus_x[event])
  evt_hcal_clus_y = np.array(hcal_clus_y[event])
  evt_hcal_clus_z = np.array(hcal_clus_z[event])
  evt_hcal_clus_E = np.array(hcal_clus_E[event])

  for it in evt_pfo_tag:
    pfo_ecal_x = evt_ecal_clus_x[evt_ecal_pfo_tag==it]
    pfo_ecal_y = evt_ecal_clus_y[evt_ecal_pfo_tag==it]
    pfo_ecal_z = evt_ecal_clus_z[evt_ecal_pfo_tag==it]
    pfo_ecal_E = evt_ecal_clus_E[evt_ecal_pfo_tag==it]
    pfo_hcal_x = evt_hcal_clus_x[evt_hcal_pfo_tag==it]
    pfo_hcal_y = evt_hcal_clus_y[evt_hcal_pfo_tag==it]
    pfo_hcal_z = evt_hcal_clus_z[evt_hcal_pfo_tag==it]
    pfo_hcal_E = evt_hcal_clus_E[evt_hcal_pfo_tag==it]*44

    if len(pfo_ecal_x)==0 and len(pfo_hcal_x)==0:
      continue
    
    pfo_x = ( ( np.append(pfo_ecal_x, pfo_hcal_x)*np.append(pfo_ecal_E, pfo_hcal_E) ).sum() / 
              ( np.append(pfo_ecal_E, pfo_hcal_E) ).sum() )
    pfo_y = ( ( np.append(pfo_ecal_y, pfo_hcal_y)*np.append(pfo_ecal_E, pfo_hcal_E) ).sum() / 
              ( np.append(pfo_ecal_E, pfo_hcal_E) ).sum() )
    pfo_z = ( ( np.append(pfo_ecal_z, pfo_hcal_z)*np.append(pfo_ecal_E, pfo_hcal_E) ).sum() / 
              ( np.append(pfo_ecal_E, pfo_hcal_E) ).sum() )
    line = pv.Line((0, 0, 0), (pfo_x, pfo_y, pfo_z))
    plotter.add_mesh(line, color=color_list(it), opacity=0.7, line_width=3)


#######################################
# Print
def print_pfo(rec_path, event):
  recfile = uproot.open(rec_path)
  pfo_tag = recfile['PFO']['pfo_tag'].array()
  ecal_pfo_tag = recfile['PFO']['ecal_pfo_tag'].array()
  ecal_clus_x = recfile['PFO']['ecal_clus_x'].array()
  ecal_clus_y = recfile['PFO']['ecal_clus_y'].array()
  ecal_clus_z = recfile['PFO']['ecal_clus_z'].array()
  ecal_clus_E = recfile['PFO']['ecal_clus_E'].array()
  hcal_pfo_tag = recfile['PFO']['hcal_pfo_tag'].array()
  hcal_clus_x = recfile['PFO']['hcal_clus_x'].array()
  hcal_clus_y = recfile['PFO']['hcal_clus_y'].array()
  hcal_clus_z = recfile['PFO']['hcal_clus_z'].array()
  hcal_clus_E = recfile['PFO']['hcal_clus_E'].array()

  evt_pfo_tag = np.array(pfo_tag[event])
  evt_ecal_pfo_tag = np.array(ecal_pfo_tag[event])
  evt_ecal_clus_x = np.array(ecal_clus_x[event])
  evt_ecal_clus_y = np.array(ecal_clus_y[event])
  evt_ecal_clus_z = np.array(ecal_clus_z[event])
  evt_ecal_clus_E = np.array(ecal_clus_E[event])
  evt_hcal_pfo_tag = np.array(hcal_pfo_tag[event])
  evt_hcal_clus_x = np.array(hcal_clus_x[event])
  evt_hcal_clus_y = np.array(hcal_clus_y[event])
  evt_hcal_clus_z = np.array(hcal_clus_z[event])
  evt_hcal_clus_E = np.array(hcal_clus_E[event])

  print("Number of PFO:", len(evt_pfo_tag))
  print("Total energy = ", evt_ecal_clus_E.sum()+evt_hcal_clus_E.sum()*44)

  for it in evt_pfo_tag:
    pfo_ecal_x = evt_ecal_clus_x[evt_ecal_pfo_tag==it]
    pfo_ecal_y = evt_ecal_clus_y[evt_ecal_pfo_tag==it]
    pfo_ecal_z = evt_ecal_clus_z[evt_ecal_pfo_tag==it]
    pfo_ecal_E = evt_ecal_clus_E[evt_ecal_pfo_tag==it]
    pfo_hcal_x = evt_hcal_clus_x[evt_hcal_pfo_tag==it]
    pfo_hcal_y = evt_hcal_clus_y[evt_hcal_pfo_tag==it]
    pfo_hcal_z = evt_hcal_clus_z[evt_hcal_pfo_tag==it]
    pfo_hcal_E = evt_hcal_clus_E[evt_hcal_pfo_tag==it]
    
    print("  For pfo", it, )
    print("    Total energy = ", pfo_ecal_E.sum()+pfo_hcal_E.sum()*44)
    print("    N ecal =", len(pfo_ecal_x))

    for ip in range(len(pfo_ecal_x)):
      print("      ({:.2f}, {:.2f}, {:.2f}, {:.3f})".format(pfo_ecal_x[ip], pfo_ecal_y[ip], pfo_ecal_z[ip], pfo_ecal_E[ip]))

    print("    N hcal =", len(pfo_hcal_x))
    for ip in range(len(pfo_hcal_x)):
      print("      ({:.2f}, {:.2f}, {:.2f}, {:.3f})".format(pfo_hcal_x[ip], pfo_hcal_y[ip], pfo_hcal_z[ip], pfo_hcal_E[ip]*44))




#######################################
sim_path  = '/cefs/higgs/zyang/cepcsoft/CEPCSW_v2.1.5.alpha/yang/y_physics/simdir/nnH_gg_WO_ISR_00001.root'
rec_path  = '/cefs/higgs/zyang/cepcsoft/CEPCSW_v2.1.5.alpha/yang/y_physics/recdir/rec_nnHgg_all_pfo3.root'
# rec_path  = '/cefs/higgs/zyang/cepcsoft/CEPCSW_v2.1.5.alpha/yang/pion/recdir/rec_pi-_5GeV_pfo2.root'
# rec_path  = '/cefs/higgs/zyang/cepcsoft/CEPCSW_v2.1.5.alpha/yang/neutron/recdir/rec_n_5GeV_pfo2.root'
# rec_path = '/cefs/higgs/zyang/cepcsoft/CEPCSW_v2.1.5.alpha/yang/kaon/recdir/rec_kL_10GeV_pfo2.root'

event = 12


############################## Bars and cells ##############################
plotter = pv.Plotter()
# plotter.set_background(color='black')
##### sim #####
draw_MC_Track(rec_path, plotter, event)
# draw_MC_h0(rec_path, plotter, event)
# draw_MC_photon(rec_path, plotter, event)
draw_vertex(sim_path, plotter, event)
draw_endpoint(sim_path, plotter, event)
##### track #####
draw_tracker_hit(sim_path, plotter, event)
draw_rec_Track(rec_path, plotter, event)
##### ECAL #####
# draw_hitbars(rec_path, plotter, event, Ecut=0.001, barCol='UV')
# draw_lm(rec_path, plotter, event, barCol='UV')
##### HCAL #####
# draw_HCAL_cells(rec_path, plotter, event)
# draw_HCAL_3Dcluster(rec_path, plotter, event)
##### other #####
draw_ECAL_profile(0.08, plotter)
# draw_ECAL_shower_center(rec_path, plotter, event)
# draw_HCAL_shower_center(rec_path, plotter, event)
# draw_extrapolate_points(rec_path, plotter, event)
# draw_pfo(rec_path, plotter, event)
# draw_pfo_direction(rec_path, plotter, event)
##### print #####
# print_pfo(rec_path, event)


# plotter.show_grid(color='silver')
plotter.camera_position = [ (15000, 5000, 23000),  # position of the camera
                            (0, 0, 0),  # focal point of the camera
                            (0, 1, 0)   # view up of the camera
                          ]
plotter.show()




############################## steps ##############################
# digi_ECAL_path = '/cefs/higgs/zyang/cepcsoft/CEPCSW_v2.1.5.alpha/yang/y_physics/digidir/ECAL_nnHgg_00001.root'
# digi_HCAL_path = '/cefs/higgs/zyang/cepcsoft/CEPCSW_v2.1.5.alpha/yang/y_physics/digidir/HCAL_nnHgg_00001.root'
# digi_ECAL_path = '/cefs/higgs/zyang/cepcsoft/CEPCSW_v2.1.5.alpha/yang/pion/digidir/ECAL_pi-_5GeV.root'
# digi_HCAL_path = '/cefs/higgs/zyang/cepcsoft/CEPCSW_v2.1.5.alpha/yang/pion/digidir/HCAL_pi-_5GeV.root'
# digi_ECAL_path = '/cefs/higgs/zyang/cepcsoft/CEPCSW_v2.1.5.alpha/yang/neutron/digidir/ECAL_n_5GeV.root'
# digi_HCAL_path = '/cefs/higgs/zyang/cepcsoft/CEPCSW_v2.1.5.alpha/yang/neutron/digidir/HCAL_n_5GeV.root'


# plotter2 = pv.Plotter()
# # plotter2.set_background(color='black')
# draw_ECAL_profile(0.05, plotter2)
# draw_steps(digi_ECAL_path, digi_HCAL_path, event, plotter2)
# plotter2.show_grid(color='silver')
# plotter2.show()