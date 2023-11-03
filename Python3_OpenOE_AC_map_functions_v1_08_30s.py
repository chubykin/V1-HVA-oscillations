
# coding: utf-8

# In[1]:

# OpenOE_AC_map_v1_01, 06/30/15
#__version__=1.02
#__version__=1.03 # removing the initial gap (acquired when Open Ephys is initialized)
"""
__version__=1.04 # 07/03/15, adding saving the averaged file into an hdf5-file
__version__=1.05 # fixing the averaging, AAC
__version__=1.06 # fixed variable initial TTL-triggered datapoints removal, zeroing time
__version__ # CSD added CSD analysis
__version__=1.07 # file output is now into the report directory

__version__=1.06 - map functions file, 03/05/2016 - adding frequency analysis


"""
# OpenOE_AC_map_functions*.py

__version__=1.02 # new: 
"""
  df2_array=df2_array[:,:-(row_number_units-trial_to_divide)]
  times2=times1[:-(row_number_units-trial_to_divide)]
"""
__version__=1.03 # old: df2_array=df2_array[:,(row_number_units-trials_to_devide):]
__version__=1.04 # same order, fixed file_list sorting. Everything is working properly. 07/27/15
__version__=1.05 # 07/29/15 real python-based CSD, https://github.com/Junji110/iCSD
__version__=1.07 # 07/29/16 - df_avg - removing the trials which have more than 1500 uV std (movement artifact)


def get_probe_map(probe):
    if probe=='128A':
        # 128A probe mapping
        probe_map={16:0.00, 0:1.00, 46:2.00,
                   17:0.01, 63:1.01, 45:2.01,
                   18:0.02, 1:1.02, 44:2.02,
                   19:0.03, 62:1.03, 43:2.03,
                   20:0.04, 2:1.04, 42:2.04,
                   21:0.05, 61:1.05, 41:2.05,
                   22:0.06, 3:1.06, 40:2.06,
                   23:0.07, 60:1.07, 39:2.07,
                   24:0.08, 4:1.08, 38:2.08,
                   25:0.09, 59:1.09, 37:2.09,
                   26:0.10, 5:1.10, 36:2.10,
                   27:0.11, 58:1.11, 35:2.11,
                   28:0.12, 6:1.12, 34:2.12,
                   29:0.13, 57:1.13, 33:2.13,
                   30:0.14, 7:1.14, 32:2.14,
                   31:0.15, 56:1.15, 47:2.15,
                   15:0.16, 8:1.16, 48:2.16,
                   14:0.17, 55:1.17, 49:2.17,
                   13:0.18, 9:1.18, 50:2.18,
                   12:0.19, 54:1.19, 51:2.19,
                   11:0.20, 10:1.20, 52:2.20,
                            53:1.21,

                   80:3.00, 64:4.00, 110:5.00,
                   81:3.01, 127:4.01, 109:5.01,
                   82:3.02, 65:4.02, 108:5.02,
                   83:3.03, 126:4.03, 107:5.03,
                   84:3.04, 66:4.04, 106:5.04,
                   85:3.05, 125:4.05, 105:5.05,
                   86:3.06, 67:4.06, 104:5.06,
                   87:3.07, 124:4.07, 103:5.07,
                   88:3.08, 68:4.08, 102:5.08,
                   89:3.09, 123:4.09, 101:5.09,
                   90:3.10,69:4.10, 100:5.10,
                   91:3.11, 122:4.11, 99:5.11,
                   92:3.12, 70:4.12, 98:5.12,
                   93:3.13, 121:4.13, 97:5.13,
                   94:3.14, 71:4.14, 96:5.14,
                   95:3.15, 120:4.15, 111:5.15,
                   79:3.16, 72:4.16, 112:5.16,
                   78:3.17, 119:4.17, 113:5.17,
                   77:3.18, 73:4.18, 114:5.18,
                   76:3.19, 118:4.19, 115:5.19,
                   75:3.20, 74:4.20, 116:5.20,
                            117:4.21,

                   }
    elif probe=='64DA':

         # 64D probe mapping, channels face me
            probe_map={47:0.00, 63:1.00, 17:2.00,
                       46:0.01, 0:1.01, 18:2.01,
                       45:0.02, 62:1.02, 19:2.02,
                       44:0.03, 1:1.03, 20:2.03,
                       43:0.04, 61:1.04, 21:2.04,
                       42:0.05, 2:1.05, 22:2.05,
                       41:0.06, 60:1.06, 23:2.06,
                       40:0.07, 3:1.07, 24:2.07,
                       39:0.08, 59:1.08, 25:2.08,
                       38:0.09, 4:1.09, 26:2.09,
                       37:0.10, 58:1.10, 27:2.10,
                       36:0.11, 5:1.11, 28:2.11,
                       35:0.12, 57:1.12, 29:2.12,
                       34:0.13, 6:1.13, 30:2.13,
                       33:0.14, 56:1.14, 31:2.14,
                       32:0.15, 7:1.15, 16:2.15,
                       48:0.16, 55:1.16, 15:2.16,
                       49:0.17, 8:1.17, 14:2.17,
                       50:0.18, 54:1.18, 13:2.18,
                       51:0.19, 9:1.19, 12:2.19,
                       52:0.20, 53:1.20, 11:2.20,
                                10:1.21
                       }
    elif probe=='58DA': ### used to be 64DA-6

         # 64D probe mapping, channels face me
            probe_map={47:0.00, 63:1.00, 17:2.00,
                       46:0.01, 0:1.01, 18:2.01,
                       45:0.02, 62:1.02, 19:2.02,
                       44:0.03, 1:1.03, 20:2.03,
                       43:0.04, 61:1.04, 21:2.04,
                       42:0.05, 2:1.05, 22:2.05,
                       41:0.06, 60:1.06, 23:2.06,
                       40:0.07, 3:1.07, 24:2.07,
                       39:0.08, 59:1.08, 25:2.08,
                       38:0.09, 4:1.09, 26:2.09,
                       37:0.10, 58:1.10, 27:2.10,
                       36:0.11, 5:1.11, 28:2.11,
                       35:0.12, 57:1.12, 29:2.12,
                       34:0.13, 6:1.13, 30:2.13,
                       33:0.14, 56:1.14, 31:2.14,
                       32:0.15, 7:1.15, 16:2.15,
                       48:0.16, 55:1.16, 15:2.16,
                       49:0.17, 8:1.17, 14:2.17,
                       50:0.18, 54:1.18, 13:2.18,
                       51:0.19, 9:1.19, 12:2.19,
                       52:0.20, 53:1.20, 11:2.20,
                                10:1.21
                      }


    elif probe == '64DB':
        # 64D probe mapping, channels face monitor
            probe_map={16:0.00, 0:1.00, 46:2.00,
                       17:0.01, 63:1.01, 45:2.01,
                       18:0.02, 1:1.02, 44:2.02,
                       19:0.03, 62:1.03, 43:2.03,
                       20:0.04, 2:1.04, 42:2.04,
                       21:0.05, 61:1.05, 41:2.05,
                       22:0.06, 3:1.06, 40:2.06,
                       23:0.07, 60:1.07, 39:2.07,
                       24:0.08, 4:1.08, 38:2.08,
                       25:0.09, 59:1.09, 37:2.09,
                       26:0.10, 5:1.10, 36:2.10,
                       27:0.11, 58:1.11, 35:2.11,
                       28:0.12, 6:1.12, 34:2.12,
                       29:0.13, 57:1.13, 33:2.13,
                       30:0.14, 7:1.14, 32:2.14,
                       31:0.15, 56:1.15, 47:2.15,
                       15:0.16, 8:1.16, 48:2.16,
                       14:0.17, 55:1.17, 49:2.17,
                       13:0.18, 9:1.18, 50:2.18,
                       12:0.19, 54:1.19, 51:2.19,
                       11:0.20, 10:1.20, 52:2.20,
                                53:1.21
                       }
    elif probe == '58DB':
        # 64D probe mapping, channels face monitor
            probe_map={16:0.00, 0:1.00, 46:2.00,
                       17:0.01, 63:1.01, 45:2.01,
                       18:0.02, 1:1.02, 44:2.02,
                       19:0.03, 62:1.03, 43:2.03,
                       20:0.04, 2:1.04, 42:2.04,
                       21:0.05, 61:1.05, 41:2.05,
                       22:0.06, 3:1.06, 40:2.06,
                       23:0.07, 60:1.07, 39:2.07,
                       24:0.08, 4:1.08, 38:2.08,
                       25:0.09, 59:1.09, 37:2.09,
                       26:0.10, 5:1.10, 36:2.10,
                       27:0.11, 58:1.11, 35:2.11,
                       28:0.12, 6:1.12, 34:2.12,
                       29:0.13, 57:1.13, 33:2.13,
                       30:0.14, 7:1.14, 32:2.14,
                       31:0.15, 56:1.15, 47:2.15,
                       15:0.16, 8:1.16, 48:2.16,
                       14:0.17, 55:1.17, 49:2.17,
                       13:0.18, 9:1.18, 50:2.18,
                       12:0.19, 54:1.19, 51:2.19,
                       11:0.20, 10:1.20, 52:2.20,
                                53:1.21
                       }
    probe_map_update=probe_map
    return probe_map_update
# number_presentations = 3
# Sampling_Rate=30000 # in Hz
# si=1.0/Sampling_Rate # sampling interval in s

#probe_map={16:0.00, 0:1.00, 46:2.00,
#           17:0.01, 63:1.01, 45:2.01,
#           18:0.02, 1:1.02, 44:2.02,
#           19:0.03, 62:1.03, 43:2.03,
#           20:0.04, 2:1.04, 42:2.04,
#           21:0.05, 61:1.05, 41:2.05,
#           22:0.06, 3:1.06, 40:2.06,
#           23:0.07, 60:1.07, 39:2.07,
#           24:0.08, 4:1.08, 38:2.08,
#           25:0.09, 59:1.09, 37:2.09,
#           26:0.10, 5:1.10, 36:2.10,
#           27:0.11, 58:1.11, 35:2.11,
#           28:0.12, 6:1.12, 34:2.12,
#           29:0.13, 57:1.13, 33:2.13,
#           30:0.14, 7:1.14, 32:2.14,
#           31:0.15, 56:1.15, 47:2.15,
#           15:0.16, 8:1.16, 48:2.16,
#           14:0.17, 55:1.17, 49:2.17,
#           13:0.18, 9:1.18, 50:2.18,
#           12:0.19, 54:1.19, 51:2.19,
#           11:0.20, 10:1.20, 52:2.20,
#                    53:1.21
#           }


# In[2]:

#import h5py
# import resampy
from scipy.fftpack import fft, ifft
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import AxesGrid
#import OpenEphys2 as OE
#import OpenEphys_V12 as OE
# 06/30/15 - version 1_01
import Python3_OpenEphys_V14 as OE
import Python3_new_OpenEphys_orig as orig
import pandas as pd
from pandas import DataFrame
import numpy as np
import glob
import os
import re
import sys
import Python3_icsd as icsd
import scipy.signal as sg
import quantities as pq
from collections import Counter
from operator import itemgetter 
import itertools 
# from resampy import filters
"""
# In[3]:

app = QApplication(sys.argv)
dialog=QFileDialog()
fname=QFileDialog.getOpenFileName(dialog,"Open an Open Ephys .continuous file in the folder to process", "", ".continuous Files (*.continuous)") # returns a QString fname    
fname=str(fname)
#folder_names=xlrd_column(fname,0)
"""

# In[15]:

# http://stackoverflow.com/questions/18874135/how-to-plot-pcolor-colorbar-in-a-different-subplot-matplotlib
# http://stackoverflow.com/questions/16595138/standalone-colorbar-matplotlib
# haven't implemented yet. 07/30/15

def ax_template_VEP():
    ax=list(range(6))
    ax[0]=plt.subplot2grid((6,4),(0,0),rowspan=2,colspan=2) # LFP all channels
    #ax[0].yaxis.set_major_locator(MaxNLocator(2)) # number of vertical major intervals, one less than max number of ticks
    ax[1]=plt.subplot2grid((6,4),(2,0),rowspan=2,colspan=2) # single channel, all trials
    ax[2]=plt.subplot2grid((6,4),(4,0),rowspan=2,colspan=2) # single channel, averaged
    #ax[2].yaxis.set_major_locator(MaxNLocator(2))
    ax[3]=plt.subplot2grid((6,4),(0,2),rowspan=2,colspan=2) # all channels, all trials
    #ax[3].yaxis.set_major_locator(MaxNLocator(2))
    ax[4]=plt.subplot2grid((6,4),(2,2),rowspan=2,colspan=2) # LFP colormap
    ax[5]=plt.subplot2grid((6,4),(4,2),rowspan=2,colspan=2) # CSD colormap
    #ax[5].yaxis.set_major_locator(MaxNLocator(2))
    #ax[6]=plt.subplot2grid((7,6),(0,4),rowspan=2,colspan=2) # PSTH mean baseline 
    #ax[6].axis["right"].set_visible(False)
    #ax[7]=plt.subplot2grid((7,6),(5,4),rowspan=2,colspan=2, frameon=False,xticks=[], yticks=[]) # PSTH mean post-cond
    
    plt.tight_layout()
    return ax

def ax_template_VEP2():
    ax=list(range(8))
    ax[0]=plt.subplot2grid((8,4),(0,0),rowspan=2,colspan=2) # all trials
    #ax[0].yaxis.set_major_locator(MaxNLocator(2)) # number of vertical major intervals, one less than max number of ticks
    ax[1]=plt.subplot2grid((8,4),(2,0),rowspan=2,colspan=2) # single channel, all trials
    ax[2]=plt.subplot2grid((8,4),(4,0),rowspan=2,colspan=2) # single channel, averaged
    #ax[2].yaxis.set_major_locator(MaxNLocator(2))
    ax[3]=plt.subplot2grid((8,4),(0,2),rowspan=2,colspan=2) # all channels, averaged
    #ax[3].yaxis.set_major_locator(MaxNLocator(2))
    ax[4]=plt.subplot2grid((8,4),(2,2),rowspan=2,colspan=2) # LFP colormap
    ax[5]=plt.subplot2grid((8,4),(4,2),rowspan=2,colspan=2) # CSD colormap
    #ax[5].yaxis.set_major_locator(MaxNLocator(2))
    ax[6]=plt.subplot2grid((8,4),(6,0),rowspan=2,colspan=2) # colorbar for LFP 
    #ax[6].axis["right"].set_visible(False)
    #ax[7]=plt.subplot2grid((6,5),(5,4),rowspan=2,colspan=2, frameon=False,xticks=[], yticks=[]) # colorbar for CSD
    ax[7]=plt.subplot2grid((8,4),(6,2),rowspan=2,colspan=2) # colorbar for CSD

    plt.tight_layout()
    return ax

# not being used right now. 07/30/15
def ax_template_VEP_grid(fig):
    grid = AxesGrid(fig, 326, # similar to subplot(326)
                  nrows_ncols = (3, 2),
                  axes_pad = 0.1,
                  label_mode = "1",
                  share_all = False,
                  cbar_location="top",
                  cbar_mode="each",
                  cbar_size="7%",
                  cbar_pad="2%",
                  )
    return grid

def splitall(path):
    allparts = []
    while 1:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        elif parts[1] == path: # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])
    return allparts

def extract_channel_number(fname):
    searchObj=re.search(r'100_CH[\d]+', fname)
    return int(searchObj.group()[6:])

# In[17]:
"""
dirname=os.path.dirname(fname)
filename=os.path.basename(fname)
filename_only,extension=os.path.splitext(filename)
dirname_split=splitall(dirname)
searchObj=re.search(r'[\d.]+\_CH',filename_only)
file_list=glob.glob(dirname+'/*'+searchObj.group()+'*'+extension)
file_list.sort()
print len(file_list)
print dirname_split


# In[7]:

get_ipython().magic(u'matplotlib notebook')
"""

def load_continuous_aligned_pad(file_list):
  print ("\nfile_list:", file_list,'\n\n')
  dirname=os.path.dirname(file_list[0])
  filename=os.path.basename(file_list[0])
  filename_only,extension=os.path.splitext(filename)
  dirname_split=splitall(dirname)
  
  """
  # file_list is sorted before it's passed to the function, 07/27/15
  # added on 7/24/15 - to filter 10*_CH*.continuous only (not AUX)
  searchObj=re.search(r'[\d.]+\_CH',filename_only)
  if searchObj:
    file1=glob.glob(dirname+'/*'+searchObj.group()+'*'+extension)
    file_list=sorted(file1, key=lambda x: extract_channel_number(x))
  else:
    print 'glob didn\'t find anything'
  """

 

  a = orig.loadContinuous(file_list[0]) 
  Sampling_Rate=float(a['header']['sampleRate'])
  si=1.0/Sampling_Rate # sampling interval in s

  samples_per_record = 1024
  

  data_np=np.array(a['data'])
  ts = a['timestamps']/samples_per_record
  

  
  
  # data_aligned = align_timestamps_pad(data_np)

  for i in file_list[1:]:
      
      a = orig.loadContinuous(i)
      
      data_np=np.vstack((data_np, a['data']))


  ls_block = find_trial_length(ts) 
  block_length = max(ls_block)

  for i in range(len(ls_block)+1): 
      if ls_block[i]>(block_length - 10): 
          gap = sum(ls_block[:i])*1024
          gap_idx = i
          break
  print (gap)

  data_ar = data_np[:,gap:]
  exp = {}
  #fs = 3e4
  #exp['duration'] = len(ts)*1024//fs
  #exp['trial_length'] = max(ls_block)*1024//fs 
  #exp['n_trials'] = exp['duration']/exp['trial_length']
  #print ('experiment info',exp)
  fs = Sampling_Rate
  exp['duration'] = len(ts)*1024//fs
  exp['trial_length'] = max(ls_block)*1024//fs 
  exp['n_trials'] = len(set(a['recordingNumber']))
#   exp['n_trials'] = exp['duration']/exp['trial_length']
  print ('experiment info',exp)

  trial_length_samples = int(exp['trial_length']*Sampling_Rate)
  number_presentations = int(exp['n_trials'])

  

  npad = ((0,0), (0,((block_length-ls_block[gap_idx])*1024)))
  data_aln = data_ar[:,:ls_block[gap_idx]*1024]
  data_aln = (np.pad(data_aln, pad_width=npad, mode='constant', constant_values=0))
  data_aln = data_aln[:,:trial_length_samples]
  start = ls_block[gap_idx]*1024
  for ii in range(len(ls_block[gap_idx+1:])):


      end = ls_block[gap_idx+ii+1]*1024 + start

      if ls_block[gap_idx+ii+1]>(block_length-10):
          npad = ((0,0), (0,((block_length-ls_block[gap_idx+ii+1])*1024)))
          tmp = data_ar[:,start:end]
          tmp = (np.pad(tmp, pad_width=npad, mode='constant', constant_values=0))
          tmp = tmp[:,:trial_length_samples]
          data_aln = np.hstack((data_aln, tmp) )
          start = end 
      else:
          start = end

  # down_fs = Sampling_Rate/30.0
  # data_aln = data_aln[:,::30] # resample to 1kHz
  # si = 1/down_fs
  df=DataFrame(data_aln)
  print ("Data dimensions:", np.shape(data_aln))


  # In[9]:

  #a=OE.loadContinuous(file_list[0])
  times=np.linspace(0,np.shape(df)[1]*si,np.shape(df)[1])
  
  return df, times, si, dirname, dirname_split, filename_only, number_presentations # last edit, 11:49, 7/24/15



def find_trial_length(ts): 
    ls=[]
    for k,g in itertools.groupby(list(enumerate(ts)),key=lambda x:x[0]-x[1]):
        ls.append(len(list(map(itemgetter(1), g)))) 
    return ls 


def df_probe_map(df,times,probe_map,n_presentations, ax, probe_column=2 ):

  df.index=df.index.map(probe_map.get) # v1_01, 06/30/15
  #df["mapping"]=df.index.map(probe_map.get) # 06/30/15. This works with drifting gratings test file, but not with S1-S2 stim
  # another alternative not optimal solution is to replace one of the data columns with mapped indices
  #df[0].put(probe_map.keys(), probe_map.values())
  #print df.head()
  

 

  #result = df.sort('mapping',ascending=False)
  result=df.sort_index(ascending=False)
  print (result.head())




  df_2=result[(result.index>=probe_column)&(result.index<(probe_column+1))]
  #plt.figure() # switching to object-oriented form for subplot reports for each recording session
  """
  ax.set_title('All channels and trials')
  for i in df_2.index:
      ax.plot(times1,(df_2.ix[i]+800*(i-2)*100)/800) #These 4 lines to plot all channels all Trials
      print i
    """

  ax.set_xlabel('Time (s)')
  ax.set_ylabel('Channel #')

  # row_number_units=df_2.ix[probe_column].count()
  # trial_units=np.fix((row_number_units-1)/number_presentations)
  # trial_to_divide=trial_units*number_presentations
  number_presentations = n_presentations
  row_number_units = np.shape(df_2)[1]
  trial_units=int(np.fix((row_number_units)/number_presentations))
  trial_to_divide=trial_units*number_presentations


  print (trial_units)
  print (row_number_units)

  print (len(df_2.index))
  print (trial_units*number_presentations)
  print ("Dividing total number of units (minus one)/number of presentations: ",(row_number_units-1)/number_presentations)

  
  df2_array=np.array(df_2)
  
  
  df2_array=df2_array[:,int(row_number_units-trial_to_divide):]
  times2=times[int(row_number_units-trial_to_divide):]
  
  #df2_array=df2_array[:,:-(row_number_units-trial_to_divide)]
  #times2=times1[:-(row_number_units-trial_to_divide)]
  return df2_array,times2, trial_units

def select_trace(df_2,times1,time0,time1, ax):
  
  times1_s=pd.Series(times1)
  #print times1_s
  ix_selected=times1_s[(times1_s>154.5)&(times1_s<158)].index
  print (ix_selected)

  df2_selected=df_2.ix[:,ix_selected]

  #plt.figure()
  
  ax.set_title('Select Trace:'+time0+'-'+time1)
  for i in df2_selected.index:
      ax.plot((df2_selected.ix[i]+800*(i-2)*100)/800)
      print (i)
  return df2_selected



def reshape_df2_array(df2_array, number_presentations, trial_units):
  print ("df2_array shape:",np.shape(df2_array))
  df3_array=np.reshape(df2_array,(np.shape(df2_array)[0],number_presentations, int(trial_units)))
  print ("df3 array shape:",np.shape(df3_array))
  return df3_array



def df_avg(df2_array, si, number_presentations, trial_units, ax):

  print ("df2_array shape:",np.shape(df2_array))
  df3_array=np.reshape(df2_array,(np.shape(df2_array)[0],number_presentations, int(trial_units)))
  print ("df3 array shape:",np.shape(df3_array))
  #df3_avg=np.mean(df3_array,0).transpose()
  #df3_avg=np.mean(df3_array,1)
  
  # DataFrames don't work with 3-dimensional arrays

  #mask = np.all(np.greater(df3_array, -1500), axis=1)
  # df3_array=df3_array[np.all(np.greater(df3_array, -1500), axis=1)]  # running our of memory with this one
  # ixs=np.argwhere(df3_array<-1500)
  # df3_array=np.delete(df3_array, list(np.unique(ixs[:,1])),axis=1)



  """

  ### remove bad trials untested, _AP
    # ls = []
    # for i in range(np.shape(df3_array)[0]):


    #     for j in range(np.shape(df3_array)[1]):

    #         if np.max(df3_array[i][j])>2000:
    #             print 'channel', i , 'trial', j, 'val',np.max(df3_array[i][j])
    #             ls.append(j) 
    #             continue
    #         elif np.min(df3_array[i][j])<-2000: 
    #             print 'channel ', i , 'trial ', j
    #             ls.append(j)

    # if len(ls)>20: # if lots of bad trials found , remove top 3
    #     count = Counter(ls)


    #     most_common = count.most_common(3)


    #     ls = [elem[0] for elem in most_common]
                
    # if len(ls)>0:
        
    #     df3_array = np.delete(df3_array, ls, axis = 1)
    #     print ls, 'trials removed'

  """









  df3_avg=np.mean(df3_array,1)

  print ("df3 array shape 2:",np.shape(df3_array))

  """
  df3_std=np.std(df3_array,1)
  print "max std:",np.max(df3_std)
  trials_movement=np.argwhere(df3_std>-1500)

  if len(trials_movement)>0:
    print "trials with movement:",trials_movement
    df3_avg=np.mean(np.delete(df3_array, trials_movement))
  """

  print ("df3_avg shape:",np.shape(df3_avg))
  #df3=DataFrame(df3_avg)
  #print 'df3.ix[0] length:',len(df3.ix[0])

  #times3=np.arange(0,np.shape(df3_avg)[1]*si,si) # times in seconds. 02/05/16 - wrong, leads to a bug, see correction with linspace
  times3=np.linspace(0,np.shape(df3_avg)[1]*si,np.shape(df3_avg)[1]) # times in second
  print ('times3 length (created from arange, df3):',len(times3))

  ax.set_title("All channels, averaged")
  for i in range(np.shape(df3_avg)[0]):
   
      #print i, "shapes time3",len(times3), "; df3_avg", np.shape(df3_avg) # troubleshooting

      #plt.plot(times2_avg,(df3.ix[i]+1000*i)/1000)
      #plt.plot(times3, (df3.ix[i]+1000*i)/1000)
      ax.plot(times3, (df3_avg[i,:]+400*i)/400)

  ax.set_xlabel('Time (s)')
  ax.set_ylabel('Channel #')

  return df3_avg, times3

def df_one_channel_analysis(df2_array, times3, number_presentations, trial_units, ax,ax1,channel_selected=7):
  # In[23]:

  # simplifying, doing it only for one channel:
  ddf=df2_array[channel_selected,:]
  print (np.shape(ddf))
  ddf2=np.reshape(ddf,(number_presentations, trial_units))
  print (np.shape(ddf2))


  # In[26]:

  for i in range(np.shape(ddf2)[0]):
      #plt.plot(times2_avg,(df3.ix[i]+1000*i)/1000)
      ax.plot(times3, (ddf2[i,:]+1000*i)/1000)
  ax.set_title('One channel, all trials')
  ax.set_ylabel('Trials')
  ax.set_xlabel('Time (s)')

  # It seems like there is not synchronous and not consistent stimulus onset/recording!!!
  # Look into this!
  # 

  # In[28]:
  ddf_avg=np.mean(ddf2,0)
  ax1.plot(times3,ddf_avg)
  ax1.set_title('One channel, averaged')
  ax1.set_ylabel(r'$\mu$V')
  ax1.set_xlabel('Time (s)')
  return ddf_avg

# In[29]:
# these two functions are added in ver 1_06 of map functions, 03/05/16

def df_one_channel_plotting(ddf2, times3, ax,ax1):
  for i in range(np.shape(ddf2)[0]):
      #plt.plot(times2_avg,(df3.ix[i]+1000*i)/1000)
      ax.plot(times3, (ddf2[i,:]+1000*i)/1000)
  ax.set_title('One channel, all trials')
  ax.set_ylabel('Trials')
  ax.set_xlabel('Time (s)')

  ddf_avg=np.mean(ddf2,0)
  ax1.plot(times3,ddf_avg)
  ax1.set_title('One channel, averaged')
  ax1.set_ylabel(r'$\mu$V')
  ax1.set_xlabel('Time (s)')
  return ddf_avg

def plot_spectrogram(x, nfft=256, Sampling_Rate=30000, Freq_top=400):#original: nfft = 256
  x_freq_top=[x[i] for i in range(0,len(x),Sampling_Rate//Freq_top)]
  
  return plt.specgram(x_freq_top, NFFT=nfft, Fs=Freq_top, window = sg.get_window(('kaiser', 14 ), nfft),
                                  cmap='jet', vmin = -10, vmax = 60,  mode = 'psd', noverlap=nfft/2)

def plot_FFT(x, ax, Sampling_Rate=30000):
  fft_out = np.fft.rfft(x)
  fft_mag = [np.sqrt(i.real**2 + i.imag**2)/len(fft_out) for i in fft_out]
  num_samples = len(x)
  rfreqs = [(i*1.0/num_samples)*Sampling_Rate for i in range(num_samples//2+1)]
  ax.plot(rfreqs[0:500], fft_mag[0:500])
  ax.set_title('FFT')
  ax.set_xlabel('FFT Bins')
  ax.set_ylabel('Magnitude')
  return ax

def downsample(trace,down):
    downsampled = sg.resample(trace,np.shape(trace)[0]/down)
    return downsampled


def pcolor(dd, th_bin, iscolorbar=1):
    dd=dd.T
    ttr=np.arange(0,np.shape(dd)[0]*th_bin,th_bin)  
    x=ttr
    y=np.arange(np.shape(dd)[1])
    X, Y = np.meshgrid(x,y)
    C=np.asarray(dd)
    extent=(min(x),max(x),min(y),max(y))
    imgplot=plt.imshow(C.T,cmap='jet',origin="lower",aspect='auto',vmax=50,vmin=-50, extent=extent,interpolation='bilinear') # added vmax=50 (Hz) to set up the max value of the colormap 
    ax = plt.gca()
    #ax.set_ylim(ax.get_ylim()[::-1]) # this and origin="lower" in imshow() is crucial for the correct orientation and order   
    if iscolorbar==1:
        plt.colorbar()
    #plt.xlim(0,2500)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlabel('Time (ms)')
    plt.ylabel('Channels')
    #plt.show()
    #plt.savefig(fname[:-4]+'_pcolor.pdf', format='pdf')
    return imgplot # returns the image for further post-processing

def df_CSD_analysis2(df3_avg, times3, ax, ax1):
  # ## CSD Analysis

  # In[31]:

  df4=np.array(df3_avg[1,:]-df3_avg[0,:])
  for i in range(np.shape(df3_avg)[0]-1):
      print (np.shape(df4))
      df4=np.vstack((df4,df3_avg[i+1,:]-df3_avg[i,:]))

  print ("df4 shape:",np.shape(df4))

  for i in range(np.shape(df3_avg)[0]-1):
      ax.plot(times3, (df4[i,:]+400*i)/400)
  #df4=np.array(df4)
  ax.set_ylabel('Trials')
  ax.set_xlabel('Time (s)')
  ax.set_title('CSD analysis, deltaX')

  print (np.shape(df4)) 
  ax1.set_title('CSD analysis, colorplot')
  #plt.imshow(df4)
  pcolor(df4,1.0/30000.0)
  return df4


def df_CSD_analysis(data,Channel_Number=21,colormap='jet',si=1.0/1000, show=True):

  #prepare lfp data for use, by changing the units to SI and append quantities,
  #along with electrode geometry and conductivities
  #lfp_data = test_data['pot1'] * 1E-3 * pq.V        # [mV] -> [V]
  lfp_data = data * 1E-3 * pq.V        # [mV] -> [V]
  #z_data = np.linspace(100E-6, 2300E-6, 23) * pq.m  # [m]
  z_data = np.linspace(100E-6, 1100E-6, Channel_Number) * pq.m  # [m]
  diam = 500E-6 * pq.m                              # [m]
  sigma = 0.3 * pq.S / pq.m                         # [S/m] or [1/(ohm*m)]
  sigma_top = 0. * pq.S / pq.m                      # [S/m] or [1/(ohm*m)]

  # Input dictionaries for each method
  delta_input = {
      'lfp' : lfp_data,
      'coord_electrode' : z_data,
      'diam' : diam,        # source diameter
      'sigma' : sigma,           # extracellular conductivity
      'sigma_top' : sigma,       # conductivity on top of cortex
      'f_type' : 'gaussian',  # gaussian filter
      'f_order' : (3, 1),     # 3-point filter, sigma = 1.
  }
  step_input = {
      'lfp' : lfp_data,
      'coord_electrode' : z_data,
      'diam' : diam,
      'sigma' : sigma,
      'sigma_top' : sigma,
      'tol' : 1E-12,          # Tolerance in numerical integration
      'f_type' : 'gaussian',
      'f_order' : (3, 1),
  }
  spline_input = {
      'lfp' : lfp_data,
      'coord_electrode' : z_data,
      'diam' : diam,
      'sigma' : sigma,
      'sigma_top' : sigma,
      'num_steps' : 201,      # Spatial CSD upsampling to N steps
      'tol' : 1E-12,
      'f_type' : 'gaussian',
      'f_order' : (15, 5),
  }
  std_input = {
      'lfp' : lfp_data,
      'coord_electrode' : z_data,
      'f_type' : 'gaussian',
      'f_order' : (3, 1),
  }


  #Create the different CSD-method class instances. We use the class methods
  #get_csd() and filter_csd() below to get the raw and spatially filtered
  #versions of the current-source density estimates.
  csd_dict = dict(
      #delta_icsd = icsd.DeltaiCSD(**delta_input),
      #step_icsd = icsd.StepiCSD(**step_input),
      spline_icsd = icsd.SplineiCSD(**spline_input),
      #std_csd = icsd.StandardCSD(**std_input), 
  )


  for method, csd_obj in csd_dict.items():
      #fig, axes = plt.subplots(3,1, figsize=(8,8))
      
      #plot LFP signal
      if show:
          f,axes=plt.subplots(1,2,figsize=(10,4))
          ax = axes[0]

          """
          im = ax.imshow(np.array(lfp_data), origin='upper', vmin=-abs(lfp_data).max(), \
                    #vmax=abs(lfp_data).max(), cmap='jet_r', interpolation='nearest')
                     vmax=abs(lfp_data).max(), cmap='jet', interpolation='nearest')
          """

          im = ax.imshow(np.array(lfp_data), origin='upper', vmin=-abs(lfp_data).max(), \
          #im = ax.imshow(np.array(lfp_data), origin='upper', vmin=-1e7, \
                    #vmax=abs(lfp_data).max(), cmap='jet_r', interpolation='nearest')
                    #vmax=abs(lfp_data).max(), cmap='jet', interpolation='nearest')
                    vmax=abs(lfp_data).max(), cmap=colormap, interpolation='nearest')
                    #vmax=1e7, cmap=colormap, interpolation='nearest')
          ax.axis(ax.axis('tight'))
          cb = plt.colorbar(im, ax=ax)
          cb.set_label('LFP (%s)' % lfp_data.dimensionality.string)
          ax.set_xticklabels([])
          ax.set_title('LFP')
          ax.set_ylabel('Channel #')
     

      #plot raw csd estimate
      
      csd = csd_obj.get_csd()
      
      """
      ax = axes[1]

      im = ax.imshow(np.array(csd), origin='upper', vmin=-abs(csd).max(), \
            #vmax=abs(csd).max(), cmap='jet_r', interpolation='nearest')
             vmax=abs(csd).max(), cmap='jet', interpolation='nearest')
      ax.axis(ax.axis('tight'))
      ax.set_title(csd_obj.name)
      cb = plt.colorbar(im, ax=ax)
      cb.set_label('CSD (%s)' % csd.dimensionality.string)
      ax.set_xticklabels([])
      ax.set_ylabel('ch #')
      """

      #plot spatially filtered csd estimate
      if show:
          ax = axes[1]
      csd = csd_obj.filter_csd(csd)
      #si=1.0/30000
      ttr=np.arange(0,np.shape(csd)[1]*si,si)
      x=ttr
      
      y=np.arange(0,np.shape(csd)[0]/10,0.1)
      #y=np.arange(np.shape(csd)[0])
      #y=np.arange(np.shape(csd)[0]-1,0,-1) # doesn't make any difference. 07/30/15

      X, Y = np.meshgrid(x,y)
      #extent=(min(x),max(x),min(y),max(y))
      extent=(min(x),max(x),max(y),min(y))


      if show:
          im = ax.imshow(np.array(csd), origin='upper', 
                     extent=extent, cmap=colormap, interpolation='nearest')
          ax.axis(ax.axis('tight'))
          ax.set_title(csd_obj.name + ', filtered')
          cb = plt.colorbar(im, ax=ax)
          cb.formatter.set_powerlimits((0, 0))
          cb.update_ticks()
          cb.set_label('CSD (%s)' % csd.dimensionality.string)
          ax.set_ylabel('Channel #')
          ax.set_xlabel('Time (s)')
          plt.show()
      #ax.set_yticks(y_ticks)
      return csd



def save_df3_avg_h5(ddf2, df3_avg,times3,si,fname):

  # ## need to save the averaged df3 and times3:
  # 
  # 

  # In[34]:

  #f = pd.HDFStore(dirname+filename_only+'_df_avg.h5')
  f = pd.HDFStore(fname)
  f['df_avg'] =  DataFrame(df3_avg)
  f['times']= pd.Series(times3)
  f['si']=pd.Series(si)
  f['layer4']=DataFrame(ddf2)

  f.close()

  # trying to determine the initial gap to remove. 06/30/15, v1_03
  """
  data1=a['data']
  t0=573600
  t1=579490
  """


  # In[11]:

  """"
  # Test that we can read the file: - works
  df_now=pd.read_hdf(dirname+filename_only+'_df_avg.h5','df_avg')
  times_now=pd.read_hdf(dirname+filename_only+'_df_avg.h5','times')
  si_now=pd.read_hdf(dirname+filename_only+'_df_avg.h5','si')
  df_now
  """


  # In[ ]:
  return None
def tf_cmw(df_res, num_freq = 40, method = 3):
    

    """
	performs complex Morlet wavelet convolution
    adapted fro Mike Cohen, AnalyzingNeuralTimeSeries Matlab code 
    (gives exactly the same result as Matlab, validated by Alex Pak)
    input: pd dataframe, n x m (trials x time)
    output: tf  n x m x k, method (total, phase/non-phase locked, itpc) x frequency x time

    range_cycles: cycled for Morlte wavelets, larger n of cycles -> better frequency precision, 
                                              smaller n of cycles -> beter temporal precision
    frex: range of frequencies
    """

    # global vars for wavelets

    base_idx = [0, 400]
    min_freq = 2
    max_freq = 40
    num_frex = num_freq
    range_cycles = [3, 10]

    # data info
    Sampling_Rate = 30000
    down_fs = Sampling_Rate/30.0
    down_sample=False
    if down_sample==True:
        base_idx = [0, 400]
        df_res = resampy.resample(df_res, Sampling_Rate, down_fs , axis=1)	        
        #df_res = resampy.resample(df_res.values, Sampling_Rate, down_fs , axis=1)	
    num_trials= df_res.shape[0]
    num_samples = df_res.shape[1]

    #frequencies vector
    frex = np.logspace(np.log10(min_freq),np.log10(max_freq),num_frex)
    #frex = np.linspace(min_freq,max_freq,num_frex)
    #time = np.linspace(0, num_samples, int(num_samples))
    time = np.linspace(0, num_samples/down_fs, int(num_samples) )

    #wavelet parameters
    s = np.divide(np.logspace(np.log10(range_cycles[0]),np.log10(range_cycles[-1]),num_frex), 2*np.pi*frex)
#     wavtime = np.linspace(-1, 1, 4*int(down_fs)+1)
    wavtime = np.linspace(-2, 2, 4*int(down_fs)+1) ## default
    half_wave = (len(wavtime)-1)/2

    #FFT parameters
    nWave = len(wavtime)
    nData = num_trials * num_samples 
    # .shape[0] - trials, shape[1] - time


    # len of convolutions
    nConv = [nWave+nData-1, nWave+nData-1 ,  nWave+num_samples-1 ]
    # container for the output
    tf = np.zeros((4, len(frex),num_samples) )


    phase_loc = df_res.mean(axis=0)
    # should be equal to 0, while averaged
    non_phase_loc = df_res-phase_loc.T

    dataX = {}
    #FFT of total data
    dataX[0] = fft( np.array(df_res).flatten(), nConv[0])
    #FFT of nonphase-locked data
    dataX[1] = fft( np.array(non_phase_loc).flatten(), nConv[1])
    #FFT of ERP (phase-locked data)
    dataX[2] = fft( phase_loc ,nConv[2])
    
    tf3d = []
    
    #main loop
    for fi in range(len(frex)):


        # create wavelet and get its FFT
        # the wavelet doesn't change on each trial...
        wavelet  = np.exp(2*1j*np.pi*frex[fi]*wavtime) * np.exp(-wavtime**2/(2*s[fi]**2))


         # run convolution for each of total, induced, and evoked
        for methodi in range(method):

            # need separate FFT 
            waveletX = fft(wavelet,nConv[methodi])
            waveletX = waveletX / max(waveletX)

            # notice that the fft_EEG cell changes on each iteration
            a_sig = ifft(waveletX*dataX[int(methodi)],nConv[int(methodi)])

            a_sig = a_sig[int(half_wave): int(len(a_sig)-half_wave)]


            if methodi<2:
                a_sig = np.reshape(a_sig, (num_trials, num_samples ))
                
                # compute power
                temppow = np.mean(abs(a_sig)**2,0)
            else:
            	#non-phase locked response
                temppow = abs(a_sig)**2
            if methodi==0:

                tf3d.append(abs(a_sig)**2) #TOTAL POWER(0), non-phase lock (1), phase lock (2)

            # db correct power, baseline +1 to make the same results as matlab range (1,501)
#             tf[methodi,fi,:] = 10*np.log10( temppow / np.nanmean(temppow[base_idx[0]:base_idx[1]+1]) )
            tf[methodi,fi,:] =  temppow / np.nanmean(temppow[base_idx[0]:base_idx[1]+1]) 

            # inter-trial phase consistency on total EEG
            if methodi==0:
                tf[3, fi, :] = abs(np.mean(np.exp(1j*np.angle(a_sig)),0))
    #contour_levels = np.arange(-10, 15, 1)
    contour_levels = np.arange(-10, 25, 1)
        

    return tf, time, frex, tf3d

def tf_cmw1(ax, df_res, num_freq = 40, method = 3, show=True):
    

    """
	performs complex Morlet wavelet convolution
    adapted fro Mike Cohen, AnalyzingNeuralTimeSeries Matlab code 
    (gives exactly the same result as Matlab, validated by Alex Pak)
    input: pd dataframe, n x m (trials x time)
    output: tf  n x m x k, method (total, phase/non-phase locked, itpc) x frequency x time

    range_cycles: cycled for Morlte wavelets, larger n of cycles -> better frequency precision, 
                                              smaller n of cycles -> beter temporal precision
    frex: range of frequencies
    """

    # global vars for wavelets

    base_idx = [0, 400]
    min_freq = .5
    max_freq = 80
    num_frex = num_freq
    range_cycles = [3, 10]

    # data info
    Sampling_Rate = 30000
    down_fs = Sampling_Rate/30.0
    down_sample=False
    if down_sample==True:
        base_idx = [0, 400]
        #df_res = resampy.resample(df_res, Sampling_Rate, down_fs , axis=1)	        
        df_res = resampy.resample(df_res.values, Sampling_Rate, down_fs , axis=1)	
    num_trials= df_res.shape[0]
    num_samples = df_res.shape[1]

    #frequencies vector
    frex = np.logspace(np.log10(min_freq),np.log10(max_freq),num_frex)
    #frex = np.linspace(min_freq,max_freq,num_frex)
    #time = np.linspace(0, num_samples, int(num_samples))
    time = np.linspace(0, num_samples/down_fs, int(num_samples) )

    #wavelet parameters
    s = np.divide(np.logspace(np.log10(range_cycles[0]),np.log10(range_cycles[-1]),num_frex), 2*np.pi*frex)
    wavtime = np.linspace(-2, 2, int(4*down_fs+1))
    half_wave = (len(wavtime)-1)/2

    #FFT parameters
    nWave = len(wavtime)
    nData = num_trials * num_samples 
    # .shape[0] - trials, shape[1] - time


    # len of convolutions
    nConv = [nWave+nData-1, nWave+nData-1 ,  nWave+num_samples-1 ]
    # container for the output
    tf = np.zeros((4, len(frex),num_samples) )


    phase_loc = df_res.mean(axis=0)
    # should be equal to 0, while averaged
    non_phase_loc = df_res-phase_loc.T

    dataX = {}
    #FFT of total data
    dataX[0] = fft( np.array(df_res).flatten(), nConv[0])
    #FFT of nonphase-locked data
    dataX[1] = fft( np.array(non_phase_loc).flatten(), nConv[1])
    #FFT of ERP (phase-locked data)
    dataX[2] = fft( phase_loc ,nConv[2])
    
    tf3d = []
    
    #main loop
    for fi in range(len(frex)):


        # create wavelet and get its FFT
        # the wavelet doesn't change on each trial...
        wavelet  = np.exp(2*1j*np.pi*frex[fi]*wavtime) * np.exp(-wavtime**2/(2*s[fi]**2))


         # run convolution for each of total, induced, and evoked
        for methodi in range(method):

            # need separate FFT 
            waveletX = fft(wavelet,nConv[methodi])
            waveletX = waveletX / max(waveletX)

            # notice that the fft_EEG cell changes on each iteration
            a_sig = ifft(waveletX*dataX[methodi],nConv[methodi])

            a_sig = a_sig[int(half_wave): int(len(a_sig)-half_wave)]


            if methodi<2:
                a_sig = np.reshape(a_sig, (num_trials, num_samples ))
                
                # compute power
                temppow = np.mean(abs(a_sig)**2,0)
            else:
            	#non-phase locked response
                temppow = abs(a_sig)**2
            if methodi==1:

                tf3d.append(abs(a_sig)**2) #TOTAL POWER(0), non-phase lock (1), phase lock (2)

            # db correct power, baseline +1 to make the same results as matlab range (1,501)
            tf[methodi,fi,:] = 10*np.log10( temppow / np.mean(temppow[base_idx[0]:base_idx[1]+1]) )
            #tf[methodi,fi,:] =  temppow / np.mean(temppow[base_idx[0]:base_idx[1]+1]) 

            # inter-trial phase consistency on total EEG
            if methodi==0:
                tf[3, fi, :] = abs(np.mean(np.exp(1j*np.angle(a_sig)),0))
    #contour_levels = np.arange(-10, 25, 1)
    contour_levels = np.arange(-10, 18, 1)
    if show==True:
	    ax.set_yscale('log')
	    ax.set_yticks(np.logspace(np.log10(min_freq),np.log10(max_freq),6))
	    ax.set_yticklabels(np.round(np.logspace(np.log10(min_freq),np.log10(max_freq),6)))
	    tf_plot = ax.contourf(time, frex, tf[0,:,:], contour_levels,  cmap = 'jet' , extend = 'both')

	    ax.set_ylabel('Frequency')
	    cb_tf = plt.colorbar(tf_plot, ax=ax)
	    cb_tf.set_label('dB from baseline')

        


    return tf, time, frex, tf3d


def bandpass_default(x, f_range, Fs, rmv_edge = True, w = 3):
    """
    Default bandpass filter
    
    Parameters
    ----------
    x : array-like 1d
        voltage time series
    f_range : (low, high), Hz
        frequency range for narrowband signal of interest
    Fs : float
        The sampling rate
    rmv_edge : bool
        if True, remove edge artifacts
    w : float
        Length of filter order, in cycles. Filter order = ceil(Fs * w / f_range[0])
        
    Returns
    -------
    x_filt : array-like 1d
        filtered time series
    taps : array-like 1d
        filter kernel
    """
    
    # Default Ntaps as w if not provided
    Ntaps = np.ceil(Fs*w/f_range[0])
    
    # Force Ntaps to be odd
    if Ntaps % 2 == 0:
        Ntaps = Ntaps + 1
    
    # Compute filter
    taps = sg.firwin(int(Ntaps), np.array(f_range) / (Fs/2.), pass_zero=False)
    
    # Apply filter
    x_filt = np.convolve(taps,x,'same')
    
    # Remove edge artifacts
    N_rmv = int(Ntaps/2.)
    if rmv_edge:
        return x_filt[N_rmv:-N_rmv], Ntaps
    else:
        return x_filt, taps
#filtering functions
def butter_bandpass(lowcut, highcut, fs, order=4):
    #lowcut is the lower bound of the frequency that we want to isolate
    #hicut is the upper bound of the frequency that we want to isolate
    #fs is the sampling rate of our data
    nyq = 0.5 * fs #nyquist frequency - see http://www.dspguide.com/ if you want more info
    low = float(lowcut) / nyq
    high = float(highcut) / nyq
    b, a = sg.butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(mydata, lowcut, highcut, fs, order=4):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = sg.filtfilt(b, a, mydata)
    return y

def circCorr(ang,line):
    n = len(ang)
    rxs = scipy.stats.pearsonr(line,np.sin(ang))
    rxs = rxs[0]
    rxc = scipy.stats.pearsonr(line,np.cos(ang))
    rxc = rxc[0]
    rcs = scipy.stats.pearsonr(np.sin(ang),np.cos(ang))
    rcs = rcs[0]
    rho = np.sqrt((rxc**2 + rxs**2 - 2*rxc*rxs*rcs)/(1-rcs**2)) #r
    r_2 = rho**2 #r squared
    pval = 1- scipy.stats.chi2.cdf(n*(rho**2),1)
    standard_error = np.sqrt((1-r_2)/(n-2))

    return rho, pval, r_2,standard_error