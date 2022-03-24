""
#Script for writing blueprint files from a string defining the protein topology
#H[n-m]L[n-m]
#Length: xxx
""
# USE: python build_blueprints.v2.py -xml template_bb+design.xml -blueresfile

import itertools
import re
import sys
import os
import copy
from Blueprint import Blueprint
from argparse import ArgumentParser
import numpy as np
from BuildBP import *

#==============================
# INPUT PARAMETERS
#==============================
parser = ArgumentParser()
parser.add_argument('-xml', type=str, help="xml template")
parser.add_argument('-pdb', type=str, help="input pdb")
parser.add_argument('-wts', type=str, help="wts file for design")
parser.add_argument('-blueresfile', action='store_true')
parser.add_argument('-resfile', type=str, default='resfile')
parser.add_argument('-prefix', type=str)
parser.add_argument('-filt_topols', type=str)
args = parser.parse_args()

template_xml = args.xml


######################################
#------------------------------
# fundamental parameters
#------------------------------

# STRAND LENGTHS
s2 = [6,8] # EVEN
s3 = [7] # ODD
s6 = [6,7,8] # EVEN/ODD
s5 = [6,7,8] # EVEN/ODD
# REGISTER SHIFTS
shift_s3_s6=[0,1,2]
shift_s5_s2=[0,1,2]
# beta-arcade loops
l5l3 = [ ['BBG','ABABB'] ]
loop_pairs={'L5':'L3'} # key: G, value: A
# hairpins
l2_abego='GG'
l4_abego='GG'
l8_abego='GG'
l2 = len(l2_abego)
l4 = len(l4_abego)
l8 = len(l8_abego)
# intra-strand twist constraint
twist_value=40 # to remove this constraint set the value to zero.
twist_value_tol=20
# S5-S6 Helix connection
l6_loops = ['B','XB','XXB'] # Abego B as the last loop residue preceding the helix is a good motif for N-ter helix capping
h1 = [4,5,6] # 4,5
l7_loops = ['X','XX'] # 1,2

##################################

prefix=args.prefix
cwd = os.getcwd()
if args.filt_topols:
	filein = open( args.filt_topols )
	selected_topols = []
	for line in filein:
		topol = line.split()[0].strip()
		top = topol[topol.index('-')+1:]
		selected_topols.append(top)
		print top

		
dic_top={}
##################################

total_tops=0 ; rejected=0
param_combinations = list(itertools.product(s2,s3,s5,s6,shift_s5_s2,shift_s3_s6,l5l3,l6_loops,h1,l7_loops))
for param_comb in param_combinations:
	s2,s3,s5,s6,shift_s5_s2,shift_s3_s6,l5l3,l6_loop,h1,l7_loop = param_comb
        # dependent paramters.
        s1=s2
        s7=s6
        s4=s3

	l5_abego, l3_abego = l5l3
	l5 = len(l5_abego) ; l3 = len(l3_abego)
	l6_abego = l6_loop
	l6 = len(l6_abego)
        l7_abego = l7_loop
        l7 = len(l7_abego)
	
        # for each combination.
        ss_len={'E1':s1,'E2':s2,'E3':s3,'E4':s4,'E5':s5,'E6':s6,'E7':s7}

	############################
        # map sidechain directionality. Initialize
        sidech={}
        for key in ss_len:
                sidech[key] = [ 0 for k in range(ss_len[key]) ]


	### set sidechain directionalities ###
        # define loop edges orientation based on abegos
	OrientSegmentsFromLoops('E2','E3',l3_abego,sidech)
        OrientSegmentsFromLoops('E4','E5',l5_abego,sidech)
	# now hairpins
	sidech['E1'][-1]=-1 ; sidech['E2'][0]=-1
	sidech['E3'][-1]=1 ; sidech['E4'][0]=1
	sidech['E6'][-1]=1 ; sidech['E7'][0]=1
	#########################################

        # alternate
        for strand in ss_len.keys():
                Alternate(sidech[strand])

	#########################################

	# Check compatiblity between non loop pairs that are face-to-face: L4//L8
	if ( sidech['E2'][-1] == sidech['E5'][0] ) and ( sidech['E3'][0] == sidech['E4'][-1] ):
		compatible=True

	else:
		compatible=False

	if not compatible:
		print 'L3 and L5 are incompatible'
		rejected+=1
		continue

	# alternate
	#for strand in ss_len.keys():
	#	Alternate(sidech[strand])


        # check compatiblity of s4s7_shift
        if sidech['E6'][0+shift_s3_s6] != sidech['E3'][-1]:
                print  'shift_s3_s6 incompatible'
		rejected+=1
                continue

	if sidech['E2'][0+shift_s5_s2] != sidech['E5'][-1]:
		print  'shift_s5_s2 incompatible'
		rejected+=1
		continue



	#if len(l6_loops_list) == 0:
	#	rejected+=1

	#for l6_abego in l6_loops_list:
	if l6 > 0:
          l6 = len(l6_abego)

          #if not ( l6_abego == 'ABABB' ) :
          #      continue
	
	  bulges={}
	  e_bulges={}

	  # overall topology
	  topol = "L[1-1]E[%i-%i]L[%i-%i]E[%i-%i]L[%i-%i]E[%i-%i]L[%i-%i]E[%i-%i]L[%i-%i]E[%i-%i]L[%i-%i]H[%i-%i]L[%i-%i]E[%i-%i]L[%i-%i]E[%i-%i]L[1-1]" %(s1,s1,l2,l2,s2,s2,l3,l3,s3,s3,l4,l4,s4,s4,l5,l5,s5,s5,l6,l6,h1,h1,l7,l7,s6,s6,l8,l8,s7,s7)


	  ss,combinations = GetCombinations(topol)

	  print topol


	  # Make directory and bluprints for each combination along with bulge positions
	  for comb in combinations:
	  	# Make the directory name
		filename = '%s-%s-H%s-%s-' %(prefix,l6_abego,h1,l7_abego); strand=0 ; resnum=0
	  for i,s in enumerate(ss):
		filename+='%s%i' %(ss[i],comb[i])		


          filename+='-shE5E2.%i-shE3E6.%i' %(shift_s5_s2,shift_s3_s6)

	  noprefix_filename = filename[filename.index('-')+1:]
	  if args.filt_topols:
		  flag=False
		  for topology in selected_topols:
			  if noprefix_filename in topology:
				flag=True
				dic_top.setdefault(topology,0)
				dic_top[topology]+=1
	  	  if not flag:
			continue


	  MakePlainBlueprint(ss,comb,'bp')
	  blue = Blueprint('bp')

	  # Bulges
	  keys = bulges.keys() # bulged strand names
	  keys.sort()
	  bpos_dic = {} # all positions considered for each bulged strand
	  for key in keys:
		for seg in blue.segments:
			if seg.id == key:
				st_len = len(seg.bp_data)
				if bulges[key] == 'n-center':		
					if st_len %2 == 0: # even strand length. Bulge must be at odd position
						bposs = range(1,st_len,2)[1:-1]
					else:
						bposs = range(2,st_len,2)[1:-1]
				elif bulges[key] == 'c-center': # it depends on the orientation of the first residue. In this case...
					bposs = range(1,st_len,2)[1:-1]
				else:
					bposs=[bulges[key]]

				bpos_dic[key] = bposs

	  # Take all bulge combinations
	  bcomb=[]
	  for key in keys:
	        bcomb.append( bpos_dic[key] )

	  bcombinations = list(itertools.product(*bcomb))

	  st=''
	  for key in e_bulges.keys():
		st+= '-bE.%s.%i' %(key,e_bulges[key])

	  for j,bulcomb in enumerate(bcombinations):
			picked_bulges={}
			pathname=filename
			for k,key in enumerate(keys):
				pathname += '-b%s.%i' %(key,bulcomb[k]) # this is the new filename
				picked_bulges[key]=bulcomb[k]

			pathname += st
			if not os.path.exists(pathname):
				os.mkdir(pathname)
			os.chdir(pathname)

			## Build blueprints
			MakeRefBlueprint(ss,comb,picked_bulges,e_bulges=e_bulges,refblue = 'bp')

			#---------------------------------
                        MakeFirstBlueprint(refblue= 'bp', segments = ['E1','L2','E2','L3','E3','L4','E4','L5','E5','L6','H1','L7','E6','L8','E7'], newblue = 'bp1', ss_pairing={'E1':['E2.A'],'E2':['E5.A'],'E3':['E4.A'],'E6':['E3.A'],'E7':['E6.A']},seg_abego={'L3':l3_abego,'L4':l4_abego,'L5':l5_abego,'L2':l2_abego,'L6':l6_abego,'L7':l7_abego,'L8':l8_abego})
			#---------------------------------

			if args.pdb:
				os.system('cp ../%s input.pdb' %(args.pdb))
			else:
				write_dummy_pdb('input.pdb')

			if args.wts:
				os.system('cp ../%s .' %(args.wts))

			# XML
		    	os.system('cp ../%s foo.xml' %(template_xml))
			xml_lines = open('../%s' %(template_xml),'r').readlines()
   			# names of blueprints must be consistent between here and the template.xml
			# Move above dir for another topology

			################################################

			# step1
			blue = Blueprint('bp1.b') ; blue.reindex_blueprint(start=1)
			fileout = open('cst1','w')

			n1 = len( blue.segment_dict['E1'].bp_data )
			n2 = len( blue.segment_dict['E2'].bp_data )
			n3 = len( blue.segment_dict['E3'].bp_data )


                        ################################################
			# Define OFFSETS based on sidechain directionality:
                        ################################################
                        if sidech['E2'][0+shift_s5_s2] == -1:
                                offset_s5_s2=1
                        else:
                                offset_s5_s2=0

                        if sidech['E6'][0+shift_s3_s6] == 1:
                                offset_s3_s6=1
                        else:
                                offset_s3_s6=0

			

                        ################################################

			# strand pairs WIHOUT SHIFT

			sthb = HbondsRegularHairpin(strand1="E1",strand2="E2",blueprint=blue) ; fileout.write(sthb) # By definition there is no offset.
			sthb = HbondsRegularHairpin(strand1="E6",strand2="E7",blueprint=blue) ; fileout.write(sthb)
			sthb = HbondsRegularHairpin(strand1="E3",strand2="E4",blueprint=blue) ; fileout.write(sthb)

			# strand pairs WITH SHIFT
                        sthb = HbondsRegularHairpin(strand1="E3",strand2="E6",blueprint=blue,shift=shift_s3_s6,offset_s1=offset_s3_s6) ; fileout.write(sthb)
                        sthb = HbondsRegularHairpin(strand1="E5",strand2="E2",blueprint=blue,shift=shift_s5_s2,offset_s1=offset_s5_s2,exclude_hbs=[-1]) ; fileout.write(sthb)

                        #for strand in ['E1','E2','E3','E4','E5','E6','E7']:
                        for strand in ['E2','E3','E4','E5','E6']:
                                seg = blue.segment_dict[strand]
                                if len(seg.bp_data) >=6 and twist_value > 0:
                                        curv_st = RegularStrandCurvature(strand=strand,level=2,blueprint=blue,global_twist=twist_value,global_twist_value=twist_value_tol)

                                        fileout.write(curv_st)


                        ################################################
			total_tops+=1
			################################################

			# proper hbonding between loop pairs
			
                        for loop1 in loop_pairs:
				seg1 = blue.segment_dict[loop1]
				abego1 = seg1.abego()
                                seg2 = blue.segment_dict[loop_pairs[loop1]]
                                abego2 = seg2.abego()
				if ('A' in abego1) and ('G' in abego2):
					idxA = abego1.index('A')
					idxG = abego2.index('G')
        	                        posA = seg1.bp_data[idxA][0]
	                                posG = seg2.bp_data[idxG][0]
					loopG = abego2


                                elif ('G' in abego1) and ('A' in abego2):
                                        idxA = abego2.index('A')
                                        idxG = abego1.index('G')
                                	posA = seg2.bp_data[idxA][0]
	                                posG = seg1.bp_data[idxG][0]
					loopG = abego1


				st = CircularHBondConstraints(posG,posA-1) ; fileout.write(st)
	    			st = CircularHBondConstraints(posG+1,posG-2) ; fileout.write(st)			


			# add proline phi psi constraints where we expect to have proline
                        arch_loops=['L2','L3','L4','L6','L7','L8']
			pro_cstfile = open('pro_cst','w')
                        for loop in arch_loops:
                                abego=''
                                for res in blue.segment_dict[loop].bp_data:
                                        abego+=res[2][1]
				if (abego == 'ABB' or abego=='ABABB' or abego=='BABB'):
					pos = blue.segment_dict[loop].bp_data[-1][0]
					st = ProlineAbegoB_PhiPsiConstraints(pos) ; fileout.write(st) ; pro_cstfile.write(st)


			fileout.close()	

			pro_cstfile.close()
		
			#############################################
			# XML REPLACEMENTS
			#############################################
			# Secondary structure selectors from blueprint SS
			sec_struct=''
			for res in blue.bp_data:
				ss = res[2][0]
				sec_struct+=ss
			XMLReplaceTagsValues(xml_lines=xml_lines,identifier='SecondaryStructure',tags=['xxx'],values=[sec_struct])

                        # RESFILE
                        if args.resfile:
                                pikaa_list = WriteResfile( blueprint='bp1.b', resfile_name=args.resfile )

				# Resfile residues in string for proper residue selection.
				XMLReplaceTagsValues(xml_lines=xml_lines,identifier='resfile_residues',tags=['xxx'],values=[pikaa_list])	
			#############################################

			#-----------------------
			# Write modified xml
			#-----------------------
			xml_out = open('input.xml','w')
			for line in xml_lines:
				xml_out.write(line)

			os.chdir('../')




print 'Total number of accepted topologies: %i' %total_tops 
print 'Total number of rejected topologies: %i' %rejected

