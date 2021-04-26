#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import sys
import os
import csv

#Getting the bond length and bond engle Written by zyl
'''
N O T E !

This script needs an additional file named 'atom_list.csv' to work

That csv file contains your selected atom names and numbers, and should be written like

C1,1
C2,5
N1,13
  .
  .
  .

Which contains name and number, seperated by comma, and has no title lines.


'''

'''
HOW TO USE:

Put this script in the same folder as CONTCAR

or POSCAR, together with your 'atom_list.csv' file,

then run this script.

'''

"""
----------------------------------------------------------
Please modify the bonding threshold in the code before use,

     and only supports POSCAR with "direct" format!
----------------------------------------------------------
"""
bonding_threshold = 2.3

#Calculate bond lenth
def calc_bond_lenth(atom1,atom2):
     return np.sqrt(np.sum(np.square(atom1 - atom2)))

#Calculate bond angle, centre at atom1
def calc_bond_angle(atom1_C,atom2_C,atom3_C):
     vector_2_1=np.array(atom2_C-atom1_C)
     vector_3_1=np.array(atom3_C-atom1_C)
     cos_angle=vector_2_1.dot(vector_3_1)/(calc_bond_lenth(atom1_C,atom2_C)*calc_bond_lenth(atom1_C,atom3_C))  #Bond angle among atom 2,atom 1, and atom 3. And atom 1 is the top of the angle!
     angle=np.arccos(cos_angle)*360/2/np.pi  #Unit conversion
     return angle

#Evaluate whether 2 bonds are connected, if so, return atom name to calculate bond angle
def isBonded(atom1,atom2,atom3,atom4):
     if atom1 == atom3:
          return [atom1,atom2,atom4]
     elif atom1 == atom4:
          return [atom1,atom2,atom3]
     elif atom2 == atom3:
          return [atom2,atom1,atom4]
     elif atom2 == atom4:
          return [atom2,atom1,atom3]
     else:
          return False


#This part reads CONTCAR and get cell vector
#Please do not change this part
file_poscar=open('CONTCAR','r')
line=file_poscar.readlines()
a1 = float(line[2].split()[0])
a2 = float(line[3].split()[0])
a3 = float(line[4].split()[0])
b1 = float(line[2].split()[1])
b2 = float(line[3].split()[1])
b3 = float(line[4].split()[1])
z1 = float(line[2].split()[2])
z2 = float(line[3].split()[2])
z3 = float(line[4].split()[2])
vector_a=np.array([a1,b1,z1])
vector_b=np.array([a2,b2,z2])
vector_c=np.array([a3,b3,z3])
vector_all=np.array([vector_a,vector_b,vector_c])

#Read atom to tackle from your file
with open('atom_list.csv','r') as csvfile:
     lines = csvfile.readlines()

#Two lists to represent both coordinates and name in a single number
atom_list = []
atom_name_list = []

#Pick your selected atoms from CONTCAR, including name and coordinates
for i in range(len(lines)):
     atom_name_list.append(lines[i].split(',')[0])
     line_number_raw = lines[i].split(',')[1]
     line_number = int(line_number_raw.rstrip())
     cell_coordinate = np.array([float(line[int(line_number)+7].split()[0]),float(line[int(line_number)+7].split()[1]),float(line[int(line_number)+7].split()[2])])
     for k in range(3):
          if cell_coordinate[k] > 0.5:
               cell_coordinate[k] = cell_coordinate[k] - 1    #This works for Periodic Boundary Conditions
     cartesian_coordinate = np.dot(cell_coordinate,vector_all)   #This turns cell coordinate to cartesian coordinate
     atom_list.append(cartesian_coordinate)


bond_list = []

#Write bond-lenth file
with open('bond_lenth.csv','w',newline='') as csvfile:
     writer = csv.writer(csvfile)
     writer.writerow(["atom1","atom2","bond lenth"])   #Title

     for i in range(len(atom_list)):
          for k in range(len(atom_list)):
               if k > i :
                    if 0 < calc_bond_lenth(atom_list[i],atom_list[k]) < bonding_threshold:    #Evaluate whether forms a bond. can be modified due to different systems
                         writer.writerow([atom_name_list[i],atom_name_list[k],calc_bond_lenth(atom_list[i],atom_list[k])])
                         bond_list.append([atom_name_list[i],atom_name_list[k]])     #This line is for calculating bond angles

#Write bond-angle file
with open('bond_angle.csv','w',newline='') as csvfile:
     writer = csv.writer(csvfile)
     writer.writerow(['atom1(centre)','atom2','atom3','bond angle'])

     for i in range(len(bond_list)):
          for k in range(len(bond_list)):
               if k > i :
                    if isBonded(bond_list[i][0],bond_list[i][1],bond_list[k][0],bond_list[k][1]) != False :
                         atom_name_1,atom_name_2,atom_name_3 = isBonded(bond_list[i][0],bond_list[i][1],bond_list[k][0],bond_list[k][1])
                         for m in range(len(atom_list)):
                              if atom_name_list[m] == atom_name_1:
                                   atom1 = atom_list[m]
                              if atom_name_list[m] == atom_name_2:
                                   atom2 = atom_list[m]
                              if atom_name_list[m] == atom_name_3:
                                   atom3 = atom_list[m]
                         writer.writerow([atom_name_1,atom_name_2,atom_name_3,calc_bond_angle(atom1,atom2,atom3)])
                         #If bonded, convert return name to coordinates, then calculate bond angle and write to file

file_poscar.close()

