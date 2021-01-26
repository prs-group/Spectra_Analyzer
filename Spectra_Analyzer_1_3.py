#!/usr/bin/python
# -*- coding: utf-8 -*-



import time

import csv
import matplotlib.pyplot as plt
from matplotlib import animation

import numpy as np
import lmfit
import re
from matplotlib.font_manager import FontProperties
import scipy
from scipy import integrate
import pandas as pd

#Version 1.3 26.01.21
#author = Marcel Ruth
#simps integration got replaced by trapz integration due to an integration error when using simps
#    _._     _, -'""`-._
#   (,-.`._,'(  MR   |\`-/
#       `-.-' \ )-`( , o o)
#             `-    \`_`"'-
#files have to be in same folder as this script 


#-------------------------
#Change Variables here 
how_often = 20 #number of time_vars, files, imports = last number + 1 so if BB_137.0270.csv is the last file, then 271
name_sample = "BB_137" #Samplenumber for "Samplenumber.xxxx.csv"
Spectra_cut = False  #cut in spectra data?
break_time_var = 19 #how long was the break time_var in h
cut_number = 124 #Last number of first data set +1, e.g. BB_137.0123.csv is the last, then choose 124
want_legend = True
want_animation = False #want animation? Gives Animation of both bands
second_plot = False
interval_time_var = 1 #time_var between measurements in h 
want_fit = False #just to test of the data is corrupted, if so turn this to false

#modelfit_parameters change them only if NaN error occours 
#model1: a*np.exp(x*k) + b
param_a = 1
param_k = -1
param_b = 1
#model2: A*(1-np.exp(-K*X))+B
param_A = 1
param_B = 0
param_K = 0.01

base_line_left = 1248.2   #left has to be the smaller number
base_line_right = 1249.5  #right has to be the smaller number

base_line_left_2 = 735.0    #left has to be the smaller number
base_line_right_2 = 738.0   #right has to be the smaller number

# name of the files, file data like .eps, .gif ... has to be given
file_format = "eps" #has to be set to the same as the format of the plot
name_plot_1 = "Ethynyhydroxycarbene_1249_Ald_1688_integrated.eps"
name_anim_1 = 'peak_decay_carbene.gif'
name_anim_2 = 'peak_decay_aldehyde.gif'

#colors best to use with color hex values, something like https://htmlcolorcodes.com
#PRS COLORS: Green = #008F00 Red = #8B2E2E Grey = #CFCFCF LightBlue = #C3E8FE Blue = #89ABC8 DarkBlue = #6388AC
color_axis_plot = "#6388AC"
color_label_plot = "#6388AC"
color_plot_data_points = "#6388AC"
color_plot_fit_curve = "#8B2E2E"

color_axis_plot_2 = "#008F00"
color_label_plot_2 = "#008F00"
color_plot_data_points_2 = "#008F00"
color_plot_fit_curve_2 = "#CFCFCF"

color_animation = "#6388AC"
color_animation_2 = "#6388AC"

#name in legend for data points
label_data_points = "Experimental Data at 855 cm$^{-1}$"
label_data_points_2 = "Experimental Data at 736 cm$^{-1}$"
size_legend = "xx-small" # x-small, normal, or numerical values like 20 can be used
#name for fit has to be changed in line 384 and 385

#axis label plot
ax1_label = "Peak-area at 855 cm$^{-1}$ / a.u. "
ax2_label = "Peak-area at 736 cm$^{-1}$ / a.u. "
x_axis_label = "Time / h"

#axis label animation
x_axis_label_animation = "wavenumber / cm$^{–1}$"
y_axis_label_animation = "intensity / a.u. "

#values for animation
low_y_animation = 0.04
high_y_animation = 0.3
low_y_animation_2 = 0.05
high_y_animation_2 = 0.42

#-------------------------
deltax = base_line_right -base_line_left

#create time_var_list
#-------------------------
interv = 0 #counters

def time_var_steps(interv):
    time_var_list = []
    if Spectra_cut == True:
        for i in range(how_often):
            if i < cut_number:
                time_var_list.append(interv)
                interv += interval_time_var
            elif i == cut_number:
                time_var_list.append(interv+break_time_var)
                interv += interval_time_var*break_time_var
            elif i > cut_number:
                time_var_list.append(interv)
                interv += interval_time_var
    else:
        for i in range(how_often):
            time_var_list.append(interv)
            interv += interval_time_var
    return time_var_list

time_var = time_var_steps(interv)
#-------------------------

#create integrated lists
#-------------------------

area_list = []
area_list_2 = []

# Print iterations progress
want_progress = True #turn off if you don't want to see the progress
if want_progress == True:
    def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
            printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
        """
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
        # Print New Line on Complete
        if iteration == total: 
            print()

if want_progress == True:       
    printProgressBar(0, how_often, prefix = 'Progress:', suffix = 'Complete', length = 50)

for i in range(how_often):
    if i < 10:
        file_name = name_sample + "." + "000" + str(i) + ".csv"
        with open(file_name, "r", encoding="utf-8") as csvfile:
            reader = csv.reader(csvfile)
            x_list = []
            y_list = []
            for row in reader:
                x_val = round(float((row[0])), 1)
                y_val = float(row[1])
                x_list.append(x_val)
                y_list.append(y_val)
        #creating values for the straight line
        x_index_for_left_side = x_list.index(base_line_left)
        x_index_for_right_side = x_list.index(base_line_right)
        if second_plot == True:
            x_index_for_left_side_2 = x_list.index(base_line_left_2)
            x_index_for_right_side_2 = x_list.index(base_line_right_2)

        y_value_for_left_side = y_list[x_index_for_left_side]
        y_value_for_right_side = y_list[x_index_for_right_side]
        if second_plot == True:
            y_value_for_left_side_2 = y_list[x_index_for_left_side_2]
            y_value_for_right_side_2 = y_list[x_index_for_right_side_2]

        straight_line_y = [y_value_for_left_side, y_value_for_right_side]
        straight_line_x = [base_line_left, base_line_right]
        if second_plot == True:
            straight_line_y_2 = [y_value_for_left_side_2, y_value_for_right_side_2]
            straight_line_x_2 = [base_line_left_2, base_line_right_2]

        area_straight = np.trapz(straight_line_y, x=straight_line_x)
        if second_plot == True:
            area_straight_2 = np.trapz(straight_line_y_2, x=straight_line_x_2)

        #creating region of intrest
        x_values_peak = x_list[x_index_for_right_side: x_index_for_left_side + 1]
        y_values_peak = y_list[x_index_for_right_side: x_index_for_left_side + 1]
        if second_plot == True:
            x_values_peak_2 = x_list[x_index_for_right_side_2: x_index_for_left_side_2 + 1]
            y_values_peak_2 = y_list[x_index_for_right_side_2: x_index_for_left_side_2 + 1]
        area_peak = np.trapz(y_values_peak, x=x_values_peak)

        if second_plot == True:
            area_peak_2 = np.trapz(y_values_peak_2, x=x_values_peak_2)

        #calculating difference
        area_effective = area_peak*(-1) - area_straight
        area_list.append(area_effective)
        if second_plot == True:
            area_effective_2 = area_peak_2*(-1) - area_straight_2
            area_list_2.append(area_effective_2)
        if i == 0:
            plot_1_high_y_limit = area_effective*1.001
            if second_plot == True:
                plot_2_low_y_limit = area_effective_2*0.975
        #clearance
        x_list.clear()
        y_list.clear()
        
    elif i >= 10 and i < 100:
        file_name = name_sample + "." + "00" + str(i) + ".csv"
        with open(file_name, "r") as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                x_val = round(float((row[0])), 1)
                y_val = float(row[1])
                x_list.append(x_val)
                y_list.append(y_val + 10)
        #creating values for the straight line
        x_index_for_left_side = x_list.index(base_line_left)
        x_index_for_right_side = x_list.index(base_line_right)
        if second_plot == True:
            x_index_for_left_side_2 = x_list.index(base_line_left_2)
            x_index_for_right_side_2 = x_list.index(base_line_right_2)

        y_value_for_left_side = y_list[x_index_for_left_side]
        y_value_for_right_side = y_list[x_index_for_right_side]
        if second_plot == True:
            y_value_for_left_side_2 = y_list[x_index_for_left_side_2]
            y_value_for_right_side_2 = y_list[x_index_for_right_side_2]

        straight_line_y = [y_value_for_left_side, y_value_for_right_side]
        straight_line_x = [base_line_left, base_line_right]
        if second_plot == True:
            straight_line_y_2 = [y_value_for_left_side_2, y_value_for_right_side_2]
            straight_line_x_2 = [base_line_left_2, base_line_right_2]

        area_straight = np.trapz(straight_line_y, x=straight_line_x)
        if second_plot == True:
            area_straight_2 = np.trapz(straight_line_y_2, x=straight_line_x_2)

        #creating region of intrest
        x_values_peak = x_list[x_index_for_right_side: x_index_for_left_side + 1]
        y_values_peak = y_list[x_index_for_right_side: x_index_for_left_side + 1]
        if second_plot == True:
            x_values_peak_2 = x_list[x_index_for_right_side_2: x_index_for_left_side_2 + 1]
            y_values_peak_2 = y_list[x_index_for_right_side_2: x_index_for_left_side_2 + 1]
        area_peak = np.trapz(y_values_peak, x=x_values_peak)
        
        if second_plot == True:
            area_peak_2 = np.trapz(y_values_peak_2, x=x_values_peak_2)

        #calculating difference
        area_effective = area_peak*(-1) - area_straight
        area_list.append(area_effective)
        if second_plot == True:
            area_effective_2 = area_peak_2*(-1) - area_straight_2
            area_list_2.append(area_effective_2)
        if i == how_often - 1:
            plot_1_low_y_limit = area_effective*0.975
            if second_plot == True:
                plot_2_high_y_limit = area_effective_2*1.01

        #clearance
        x_list.clear()
        y_list.clear()
        
    elif i >= 100 and i < 1000:
        file_name = name_sample + "." + "0" + str(i) + ".csv"
        with open(file_name, "r") as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                x_val = round(float((row[0])), 1)
                y_val = float(row[1])
                x_list.append(x_val)
                y_list.append(y_val)
        #creating values for the straight line
        x_index_for_left_side = x_list.index(base_line_left)
        x_index_for_right_side = x_list.index(base_line_right)
        if second_plot == True:
            x_index_for_left_side_2 = x_list.index(base_line_left_2)
            x_index_for_right_side_2 = x_list.index(base_line_right_2)

        y_value_for_left_side = y_list[x_index_for_left_side]
        y_value_for_right_side = y_list[x_index_for_right_side]
        if second_plot == True:
            y_value_for_left_side_2 = y_list[x_index_for_left_side_2]
            y_value_for_right_side_2 = y_list[x_index_for_right_side_2]

        straight_line_y = [y_value_for_left_side, y_value_for_right_side]
        straight_line_x = [base_line_left, base_line_right]
        if second_plot == True:
            straight_line_y_2 = [y_value_for_left_side_2, y_value_for_right_side_2]
            straight_line_x_2 = [base_line_left_2, base_line_right_2]

        area_straight = np.trapz(straight_line_y, x=straight_line_x)
        if second_plot == True:
            area_straight_2 = np.trapz(straight_line_y_2, x=straight_line_x_2)

        #creating region of intrest
        x_values_peak = x_list[x_index_for_right_side: x_index_for_left_side + 1]
        y_values_peak = y_list[x_index_for_right_side: x_index_for_left_side + 1]
        if second_plot == True:
            x_values_peak_2 = x_list[x_index_for_right_side_2: x_index_for_left_side_2 + 1]
            y_values_peak_2 = y_list[x_index_for_right_side_2: x_index_for_left_side_2 + 1]
        area_peak = np.trapz(y_values_peak, x=x_values_peak)
        
        if second_plot == True:
            area_peak_2 = np.trapz(y_values_peak_2, x=x_values_peak_2)

        #calculating difference
        area_effective = area_peak*(-1) - area_straight
        area_list.append(area_effective)
        if second_plot == True:
            area_effective_2 = area_peak_2*(-1) - area_straight_2
            area_list_2.append(area_effective_2)
        if i == how_often - 1:
            plot_1_low_y_limit = area_effective*0.975
            if second_plot == True:
                plot_2_high_y_limit = area_effective_2*1.01

        #clearance
        x_list.clear()
        y_list.clear()
              
    elif i >= 1000:   
        file_name = name_sample + "." + str(i) + ".csv"
        with open(file_name, "r") as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                x_val = round(float((row[0])), 1)
                y_val = float(row[1])
                x_list.append(x_val)
                y_list.append(y_val)
        #creating values for the straight line
        x_index_for_left_side = x_list.index(base_line_left)
        x_index_for_right_side = x_list.index(base_line_right)
        if second_plot == True:
            x_index_for_left_side_2 = x_list.index(base_line_left_2)
            x_index_for_right_side_2 = x_list.index(base_line_right_2)

        y_value_for_left_side = y_list[x_index_for_left_side]
        y_value_for_right_side = y_list[x_index_for_right_side]
        if second_plot == True:
            y_value_for_left_side_2 = y_list[x_index_for_left_side_2]
            y_value_for_right_side_2 = y_list[x_index_for_right_side_2]

        straight_line_y = [y_value_for_left_side, y_value_for_right_side]
        straight_line_x = [base_line_left, base_line_right]
        if second_plot == True:
            straight_line_y_2 = [y_value_for_left_side_2, y_value_for_right_side_2]
            straight_line_x_2 = [base_line_left_2, base_line_right_2]

        area_straight = np.trapz(straight_line_y, x=straight_line_x)
        if second_plot == True:
            area_straight_2 = np.trapz(straight_line_y_2, x=straight_line_x_2)

        #creating region of intrest
        x_values_peak = x_list[x_index_for_right_side: x_index_for_left_side + 1]
        y_values_peak = y_list[x_index_for_right_side: x_index_for_left_side + 1]
        if second_plot == True:
            x_values_peak_2 = x_list[x_index_for_right_side_2: x_index_for_left_side_2 + 1]
            y_values_peak_2 = y_list[x_index_for_right_side_2: x_index_for_left_side_2 + 1]
        area_peak = np.trapz(y_values_peak, x=x_values_peak)
        
        if second_plot == True:
            area_peak_2 = np.trapz(y_values_peak_2, x=x_values_peak_2)

        #calculating difference
        area_effective = area_peak*(-1) - area_straight
        area_list.append(area_effective)
        if second_plot == True:
            area_effective_2 = area_peak_2*(-1) - area_straight_2
            area_list_2.append(area_effective_2)
        if i == how_often - 1:
            plot_1_low_y_limit = area_effective*0.975
            if second_plot == True:
                plot_2_high_y_limit = area_effective_2*1.01

        #clearance
        x_list.clear()
        y_list.clear()
    else:
        continue
    if want_progress == True:
        time.sleep(0.1)
        # Update Progress Bar
        printProgressBar(i + 1, how_often, prefix = 'Progress:', suffix = 'Complete', length = 50)
#-------------------------          


#Fit function
#-------------------------  
if want_fit == True:
    def exponential(x, a, k, b):
        return a*np.exp(x*k) + b

    ExponentialModel = lmfit.Model(exponential)
    result = ExponentialModel.fit(area_list, x=time_var, a=param_a, k=param_k, b=param_b)
    #-------------------------

    #Fit function2
    #-------------------------
    if second_plot == True:
        def exponential_2(X, A, K, B):
            return  A*(1-np.exp(-K*X))+B

        exponentialModel_2 = lmfit.Model(exponential_2)
        result_2 = exponentialModel_2.fit(area_list_2, X=time_var, A=param_A, B=param_B, K=param_K)


    #search for parameters
    #-------------------------
    input_text = result.fit_report()
    k_crude = re.findall(r"k:[\n ]?[+-]?[\n ]?[0-9]*[\n ]*[0-9]*[\n ]*[.]?[\n ]*[0-9]+[\n ]*[0-9]+", input_text)
    k_calcul = float(k_crude[-1].replace("k: ", ""))
    k_const = round(k_calcul, 4)
    a_crude = re.findall(r"a:[\n ]?[+-]?[\n ]?[0-9]*[\n ]*[0-9]*[\n ]*[.]?[\n ]*[0-9]+[\n ]*[0-9]+", input_text)
    a_const = round(float(a_crude[-1].replace("a: ", "")), 4)
    b_crude = re.findall(r"b:[\n ]?[+-]?[\n ]?[0-9]*[\n ]*[0-9]*[\n ]*[.]?[\n ]*[0-9]+[\n ]*[0-9]+", input_text)
    b_const = round(float(b_crude[-1].replace("b: ", "")), 4)
    t_half = np.log(2)/(-k_const)
    t_half_round = round(t_half, 1)
    print(result.fit_report())

    if second_plot == True:
        input_text2 = result_2.fit_report() #search for parameters
        K_crude = re.findall(r"K:[\n ]?[+-]?[\n ]?[0-9]*[\n ]*[0-9]*[\n ]*[.]?[\n ]*[0-9]+[\n ]*[0-9]+", input_text2)
        K_calcul = float(K_crude[-1].replace("K: ", ""))
        K_const = round(K_calcul, 4)
        A_crude = re.findall(r"A:[\n ]?[+-]?[\n ]?[0-9]*[\n ]*[0-9]*[\n ]*[.]?[\n ]*[0-9]+[\n ]*[0-9]+", input_text2)
        A_const = round(float(A_crude[-1].replace("A: ", "")), 4)
        B_crude = re.findall(r"B:[\n ]?[+-]?[\n ]?[0-9]*[\n ]*[0-9]*[\n ]*[.]?[\n ]*[0-9]+[\n ]*[0-9]+", input_text2)
        B_const = round(float(B_crude[-1].replace("B: ", "")), 4)
        t_half2 = np.log(2)/(-K_const)
        t_half_round2 = round(-t_half2, 1)
        print(result_2.fit_report())


#Plot
#-------------------------
if want_fit == True:
    label_fit = "Best Fit: " + str(a_const) +" · exp(" + str(k_const) + "0 · $t$) + " + str(b_const) + "\nHalf-life = " + str(t_half_round) + " h"
    if second_plot == True:
        label_fit_2 = "Best Fit: " + str(A_const) +" · (1 - exp(-" + str(K_const) + " · $t$) + " + str(B_const) + "\nHalf-life = " + str(t_half_round2) + " h"

max_time_var = time_var[-1]
#-------------------------

fig, ax1 = plt.subplots()
color = color_axis_plot
ax1.set_xlabel(x_axis_label)
ax1.set_ylabel(ax1_label, color=color_label_plot)
plt1 = ax1.plot(time_var, area_list, 'o', markersize=1, color=color_plot_data_points, label=label_data_points)
if want_fit == True:
    plt2 = ax1.plot(time_var, result.best_fit, color=color_plot_fit_curve,  linewidth=2, linestyle=':', label=label_fit)
ax1.tick_params(axis='y', labelcolor=color)
ax1.axis([0, max_time_var, plot_1_low_y_limit, plot_1_high_y_limit])

if second_plot == True:
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = color_axis_plot_2
    ax2.set_ylabel(ax2_label, color=color_label_plot_2)  # we already handled the x-label with ax1
    plt3 = ax2.plot(time_var, area_list_2,'o', markersize=1, color=color_plot_data_points_2, label=label_data_points_2)
    if want_fit == True:
        plt4 = ax2.plot(time_var, result_2.best_fit, color=color_plot_fit_curve_2,  linewidth=2, linestyle=':', label=label_fit_2)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.axis([0, max_time_var, plot_2_low_y_limit, plot_2_high_y_limit])
    if want_legend == True and want_fit == True:
        fontP = FontProperties() #legend if 4 plots needed
        fontP.set_size(size_legend)
        plts = plt1 + plt2 + plt3 + plt4
        labs = [l.get_label() for l in plts]
        plt.legend(plts, labs, loc=7, prop=fontP)
else:
    if want_legend == True and want_fit == True: #legend if only 2 plots needed
        fontP = FontProperties()
        fontP.set_size(size_legend)
        plts = plt1 + plt2
        labs = [l.get_label() for l in plts]
        plt.legend(plts, labs, loc=7, prop=fontP)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig(name_plot_1, format=file_format)

#-------------------------
if want_fit == True:
    print("t_half = " +  str(t_half))
    if second_plot == True:
        print("t_half = " +  str(-t_half2))
#-------------------------
#ANIMATIONS

if want_animation == True:
    fig = plt.figure()
    ax = plt.axes(xlim=(base_line_right, base_line_left), ylim=(low_y_animation, high_y_animation))
    line, = ax.plot([], [], lw=2, color=color_animation)

    # initialization function: plot the background of each frame
    def init():
        line.set_data([], [])
        return line,

    # animation function.  This is called sequentially #get the values for the plot from csv files
    def animate(i):
        if i < 10:
            file_name = name_sample + "." + "000" + str(i) + ".csv"
            with open(file_name, "r") as csvfile:
                reader = csv.reader(csvfile)
                x_list = []
                y_list = []
                for row in reader:
                    x_val = float((row[0]))
                    y_val = float(row[1])
                    x_list.append(x_val)
                    y_list.append(y_val)
            del x_val
            del y_val

        elif i >= 10 and i < 100:
            file_name = name_sample + "." + "00" + str(i) + ".csv"
            with open(file_name, "r") as csvfile:
                reader = csv.reader(csvfile)
                x_list = []
                y_list = []
                for row in reader:
                    x_val = float((row[0]))
                    y_val = float(row[1])
                    x_list.append(x_val)
                    y_list.append(y_val)
            del x_val
            del y_val

        elif i >= 100 and i < 1000:
            file_name = name_sample + "." + "0" + str(i) + ".csv"
            with open(file_name, "r") as csvfile:
                reader = csv.reader(csvfile)
                x_list = []
                y_list = []
                for row in reader:
                    x_val = float((row[0]))
                    y_val = float(row[1])
                    x_list.append(x_val)
                    y_list.append(y_val)
            del x_val
            del y_val

        elif i >= 1000:
            file_name = name_sample + "." + str(i) + ".csv"
            with open(file_name, "r") as csvfile:
                reader = csv.reader(csvfile)
                x_list = []
                y_list = []
                for row in reader:
                    x_val = float((row[0]))
                    y_val = float(row[1])
                    x_list.append(x_val)
                    y_list.append(y_val)
            del x_val
            del y_val

        x = x_list
        y = y_list
        ax.set_xlabel(x_axis_label_animation)
        ax.set_ylabel(y_axis_label_animation)
        line.set_data(x, y)
        return line,

    #plot1
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=how_often, interval=20, blit=True)

    anim.save(name_anim_1, fps=30)

    #plot 2
    fig = plt.figure()
    ax = plt.axes(xlim=(base_line_right_2, base_line_left_2), ylim=(low_y_animation_2, high_y_animation_2))
    line, = ax.plot([], [], lw=2, color=color_animation_2)
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=how_often, interval=20, blit=True)
    anim.save(name_anim_2, fps=30)


df = pd.DataFrame(data={"time": time_var, "integral_band1": area_list})
df.to_csv("integral_values_band1.csv", sep=",", index=False)

if second_plot == True:
    df = pd.DataFrame(data={"time": time_var, "integral_band2": area_list_2})
    df.to_csv("integral_values_band2.csv", sep=",", index=False)

#END OF SCRIPT     
print("""
    _._     _,-'""`-._
   (,-.`._,'(  MR   |\`-/|
       `-.-' \ )-`( , o o)
             `-    \`_`"'-
""")


  