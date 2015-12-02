#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from numpy.random import *

def get_text_positions(x_data, y_data, txt_width, txt_height):
  """ define text positions for plotting
  Original author and descriptinons are given here:
  http://stackoverflow.com/questions/8850142/matplotlib-overlapping-annotations
  """
  a = zip(y_data, x_data)
  text_positions = y_data.copy()
  for index, (y, x) in enumerate(a):
    local_text_positions = [i for i in a if i[0] > (y - txt_height) 
	        and (abs(i[1] - x) < txt_width * 2) and i != (y,x)]
    if local_text_positions:
      sorted_ltp = sorted(local_text_positions)
      if abs(sorted_ltp[0][0] - y) < txt_height: #True == collision
        differ = np.diff(sorted_ltp, axis=0)
        a[index] = (sorted_ltp[-1][0] + txt_height, a[index][1])
        text_positions[index] = sorted_ltp[-1][0] + txt_height
        for k, (j, m) in enumerate(differ):
          #j is the vertical distance between words
          if j > txt_height * 2: #if True then room to fit a word in
            a[index] = (sorted_ltp[k][0] + txt_height, a[index][1])
            text_positions[index] = sorted_ltp[k][0] + txt_height
            break
  return text_positions

def text_plotter(x_data, y_data, text_content, text_positions, axis,txt_width,txt_height, color):
    for x,y,s,t in zip(x_data, y_data, text_content, text_positions):
        axis.text(x - txt_width, 1.01*t, s, rotation=0, color=color)
        if y != t:
            axis.arrow(x, t,0,y-t, color='grey',alpha=0.3, width=0.01, 
                       head_width=0.2, head_length=txt_width*0.5, 
                       zorder=0,length_includes_head=True)

def makePlot(x_data, y_data, tags, filename, plottitle, labelx, labely, functionlabel, textcolor):
  #random test data:
  #x_data = random_sample(100)
  #y_data = random_integers(10,50,(100))

  #GOOD PLOT:
  fig2 = plt.figure()
  ax2 = fig2.add_subplot(111)
  ax2.plot(x_data, y_data, 's'+textcolor, markersize=8, label=functionlabel)
  plt.xlabel(labelx)
  plt.ylabel(labely)
  ax2.grid()
  plt.title(plottitle)
  plt.legend()
  #set the bbox for the text. Increase txt_width for wider text.
  txt_height = 0.08*(plt.ylim()[1] - plt.ylim()[0])
  txt_width = 0.01*(plt.xlim()[1] - plt.xlim()[0])
  #Get the corrected text positions, then write the text.
  text_positions = get_text_positions(x_data, y_data, txt_width, txt_height)
  text_plotter(x_data, y_data, tags, text_positions, ax2, txt_width, txt_height, textcolor)

  plt.ylim(0,max(text_positions)+2*txt_height)
  #plt.xlim(-0.1,1.1)
  
  #plt.show()
  plt.savefig(filename)
  return
  
# Main Program

#makePlot()