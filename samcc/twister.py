import pandas as pd
import sys
import os
import math
import time
import glob
import itertools
import numpy as np
from Bio.PDB.PDBParser import PDBParser
import matplotlib.pyplot as plt
import seaborn as sns
if sys.version_info[0] == 3:
	from io import StringIO
else:
	import StringIO

np.seterr(all='raise')

def parse(fname, drop=[], silent=False):
	try:
		tfile = open(fname).read()

		seq = tfile[tfile.find("********* Sequence and heptad assignment *********"):tfile.find("********* Coiled-coil parameters per residue *********")]

		s=[]
		i=3
		seq = seq.split('\n')
		while i<len(seq):
			#print i, s
			s+= list(seq[i].strip().replace(' ', ''))
			i+=4

		s = s[1:-1]

		tfile = StringIO.StringIO(tfile[tfile.find("Res cc_phase"):tfile.find("Ave")].replace(" R ", "   "))
		df = pd.read_csv(tfile, error_bad_lines=False, delim_whitespace=True, header=0, skiprows=0, index_col=False, engine='c', dtype={"total_score": np.float64})

		#print fname

		df['seq'] = s

	except pd.io.common.EmptyDataError:
		if not silent: print('Error reading %s. Empty file?' % fname)
		return None
	except IOError:
		if not silent: print('Error reading %s. File does not exist.' % fname)
		return None

	else:
		for d in drop:
			df.drop(d, axis=1, inplace=True)

		# skip log files that contain nan values
		if df.isnull().values.any():
			if not silent: print('Error reading %s. File constains NaN.' % fname)
			return None
		#df['P'] = df['res/tur']/(1-(df['res/tur']*df.cc_dang/360.0))
		df['P'] = 3.627/(1-(3.627*df.cc_dang/360.0))
		try:
			df['A'] = np.degrees(np.arcsin(np.radians(df['cc_dang']) * df['cc_rad']) / df['a_ris'])
		except FloatingPointError:
			if not silent: print('Can''t calculate A %s' % fname)
			#return None

		return df

def average(df, fname="", allCrick=False):

	cr = df.Cr_ang.tolist()
	if not allCrick:
		cr = cr[0]

	fit = {'tw_P':df.P.mean(),\
		   'tw_p':df['res/tur'].mean(),\
		   'tw_w0':df.cc_dang.mean(),\
		   'tw_R0':df.cc_rad.mean(),\
		   'tw_R1':df.a_rad.mean(),\
		   'tw_d':df.a_ris.mean(),\
		   'tw_ccrise':df.cc_rise.mean(),\
		   'tw_ccpit':df.cc_pit.mean(),\
		   'tw_w1':df.a_dang.mean(),\
		   #'tw_ph1':float(df.iloc[[0]]['Cr_ang']),\
		   'tw_ph1':cr,\
		   }

	fit.update({'P_std':df.P.std(),\
		   'tw_p_std':df['res/tur'].std(),\
		   'tw_w0_std':df.cc_dang.std(),\
		   'tw_R0_std':df.cc_rad.std(),\
		   'tw_R1_std':df.a_rad.std(),\
		   'tw_d_std':df.a_ris.std(),\
		   'tw_ccrise_std':df.cc_rise.std(),\
		   })

	if fname!="":
		fit['id'] = fname

	return fit

def getrange(pdbfile):
	p = PDBParser(PERMISSIVE=1)
	s = p.get_structure("", pdbfile)
	chnames=[]
	pf, pt= 0, 1000000
	for c in s[0]: # iterate over chains

		cl = list(c) # list of all residues in a chain

		#if len(cl)<25: continue # skip short chains

		chnames.append(c.get_id())
		c = cl

		pos = 0

		assert (c[0].get_full_id()[3][0])==' ' # first residue is not an HET atom

		while c[pos].get_full_id()[3][0]==' ' and pos<len(c)-1:
			pos+=1

		pos-=1

		f,t = c[0].get_id()[1], c[pos].get_id()[1]

		#print f,t

		if f>pf: pf=f
		if t<pt: pt=t



	return pf, pt, ' '.join(chnames)


def run_twister_local(pdbfile, force=False):
	logfile = pdbfile.replace('.pdb', '.log')
	f, t, chnames = getrange(pdbfile)

	if os.path.isfile(logfile) and not force:
		return logfile

	os.system("~/apps/twister/runscript %s %s %s %s" % (pdbfile, f, t))
	return logfile

def run_twister_remote(pdbfile, force=False, f=None, t=None, chnames=""):
	logfile = pdbfile.replace('.pdb', '.log')
	resfile = pdbfile.replace('.pdb', '.res')
	#os.system('rm %s' % logfile)

	if f==None and t==None:
		f, t, chnames = getrange(pdbfile)
	else:
		if chnames=="":
			_, _, chnames = getrange(pdbfile)


	if os.path.isfile(logfile) and not force:
		return f,t
		#return logfile


	# edi
	#os.system('scp %s sdunin@edi01:/home/users/sdunin/apps/twister' % pdbfile)
	#cmd='ssh sdunin@edi01 "cd /home/users/sdunin/apps/twister/; ./runscript %s %s %s %s"' % (pdbfile, f, t, chnames)
	#print cmd
	#os.system(cmd)
	#os.system('scp sdunin@edi01:/home/users/sdunin/apps/twister/%s .' % logfile)

	os.system('scp %s sdh@hades.cent.uw.edu.pl:/home/sdh/apps/twister' % pdbfile)
	cmd='ssh sdh@hades.cent.uw.edu.pl "cd /home/sdh/apps/twister/; ./runscript %s %s %s %s"' % (pdbfile, f, t, chnames)
	print(cmd)
	os.system(cmd)
	os.system('scp sdh@hades.cent.uw.edu.pl:/home//sdh/apps/twister/%s .' % logfile)


	time.sleep(1)
	#return logfile
	# Edi 01
	#os.system('scp -q %s sdunin@edi01:/home/users/sdunin/apps/twister/' % pdbfile)
	#os.system('ssh sdunin@edi01 "cd /home/users/sdunin/apps/twister/; ./runscript %s"' % pdbfile)
	#os.system('scp -q sdunin@edi01:/home/users/sdunin/apps/twister/%s .' % logfile)
	#os.system('ssh sdunin@edi01 "cd /home/users/sdunin/apps/twister/; rm %s; rm %s; rm %s"' % (pdbfile, logfile, resfile))

	return f,t

def adj(ang):
	ang += 180.0
	if ang>360:
		c = int(ang / 360.0)
		return ang - (360*c) - 180.0
	else:
		return ang - 180.0

def adj2(ang):
	ang = adj(ang)
	if ang < -180:
		ang += 360
	return ang

def calc_crick_ang_dev(twister, exp_helix_crick, firstpos, lastpos, force_start_Cr_ang_pos=None):

	if firstpos == lastpos:
		assert lastpos == None


	# define Crick angle of the starting position
	if force_start_Cr_ang_pos==None:
		if lastpos==None:
			firstpos = 0
		if type(twister) is list:
			start_Cr_ang = twister[firstpos]
		else:
			start_Cr_ang = twister.iloc[firstpos]['Cr_ang']
		start_Cr_ang_pos, name, _ = crick_to_pos(start_Cr_ang, exp_helix_crick)
	else:
		start_Cr_ang_pos = force_start_Cr_ang_pos

	cr = itertools.cycle(exp_helix_crick)

	# skip to the first repeat pos
	for n in range(start_Cr_ang_pos):
		next(cr)

	if lastpos==None:
		laspos = len(twister)
	if type(twister) is list:
		data = twister[firstpos:lastpos]
	else:
		data = twister['Cr_ang'].iloc[firstpos:lastpos]

	Cr_ang_dev = [adj2(c - next(cr)) for c in data]

	m,b = np.polyfit(range(len(Cr_ang_dev)), Cr_ang_dev, 1)

	return Cr_ang_dev, m, b

def gen_expected_crick_angles(P, rep_len, start_ph1, ap=False):
	step = 360 / P
	if ap:
		sign=-1
	else:
		sign=1
	return [adj(start_ph1+(sign * i * float(step))) for i in range(rep_len)]

def crick_to_pos(start_Cr_ang, exp_helix_crick):

	diff = [abs(adj(start_Cr_ang-i)) for i in exp_helix_crick]
	mindiff = min(diff)

	start_Cr_ang_pos = diff.index(mindiff)

	try:
		name = chr(97+start_Cr_ang_pos)
	except ValueError:
		name = '?'

	return start_Cr_ang_pos, name, mindiff

	#print "region %s-%s, starting position %s (%s degrees)" % (firstpos, firstpos+fragmentlen, chr(97+start_Cr_ang_pos), exp_helix_crick[start_Cr_ang_pos])

	#print start_Cr_ang_pos
	#print chr(97+start_Cr_ang_pos)

	#print 'measured starting ph1', start_Cr_ang
	#print 'closest ideal ph1', exp_helix_crick[start_Cr_ang_pos]



def measure_axial_rotation(logname, P=3.5, replen=7, start_ph1=19.5):

	exp_ahelix_crick = gen_expected_crick_angles(P, replen, start_ph1)

	id = logname[:-4].split('/')[-1]
	df = parse(logname)
	if not type(df) is pd.DataFrame:
		return None, None, None

	#print df

	start_Cr_ang = df.Cr_ang.iloc[0]
	start_Cr_ang_pos, name, _ = crick_to_pos(start_Cr_ang, exp_ahelix_crick)
	crdev = calc_crick_ang_dev(df, exp_ahelix_crick, 0, len(df))[0]

	df['crdev'] = crdev

	return id, df, name



def measure_axial_rotations_in_dir(dir, P=3.5, replen=7, start_ph1=19.5):

	"""
	parses all twister log files in a directory and caluclates helix axial rotations
	"""
	res = {}

	for logname in glob.glob(dir+'/*.log'):
		id, df, name = measure_axial_rotation(logname, P=P, replen=replen, start_ph1=start_ph1)
		if not id==None:

			twister_std_crdev = np.std(df.crdev)
			twister_mean_crdev = np.mean(df.crdev)

			res[id] = {'twister_std_crdev': twister_std_crdev, 'twister_mean_crdev': twister_mean_crdev, 'heptad':name}

	df = pd.DataFrame.from_dict(res, orient='index')
	df['description'] = df.index

	return df


# plotting

def do_plots(dfs, colors):
	TICKFONTSIZE = 10
	AXESFONTSIZE = 15

	fig, axes = plt.subplots(nrows=1, sharey=False, sharex=False, figsize=(8, 5))
	for df, color in zip(dfs, colors):

		#print df

		#ax1 = df.crdev.rolling(window=5,center=False).mean()
		print(df.crdev.std())
		axes.plot(df.crdev.tolist(), c=color, linewidth=4.0)
		#axes.axhline(ls='--', color='0.5', y=df.crdev.mean(), linewidth=2.0)

		#for i in [17/5., 7/2., 25/7., 18/5., 11/3., 15/4., 19/5.]:
		#	axes.axhline(ls='--', color='0.5', y=i, linewidth=2.0)

		axes.set_xlabel(r'Position in the sequence', fontsize=AXESFONTSIZE)
		axes.set_ylabel(r'Helix axial rotation (degrees)', fontsize=AXESFONTSIZE)



		axes.set_ylim(-26, 26)
		#plt.yticks(np.arange(-26, 26, 2.0))

		for tick in axes.xaxis.get_major_ticks():
			tick.label.set_fontsize(TICKFONTSIZE)

		for tick in axes.yaxis.get_major_ticks():
			tick.label.set_fontsize(TICKFONTSIZE)
	plt.show()
	#plt.savefig('/Users/sdh/Desktop/twister.png')





if __name__ == "__main__":


	pdbfile = sys.argv[1]

	"""
	try:
		import renumber_PDB
		renumber_PDB.write_structure(pdbfile)
	except:
		pass
	"""

	run_twister_remote(pdbfile)

	logfile = pdbfile.replace('.pdb', '.log')
	#df = parse(logfile)
	id, df, name=measure_axial_rotation(logfile)


	do_plots([df], ['red'])
