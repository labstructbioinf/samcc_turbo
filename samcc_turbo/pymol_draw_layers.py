def save_layers_to_pymol(pdbpath, layer_points, savepath, suffix, pymol_version=2.0, helices_axis=None, color_selection=False, helix_order=None, ppo=None):
	''' save each layer as a set of distances between point determining the layer '''

	#FIXME clean from devel code, add docs

	import __main__
	__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
	import sys, time, os, itertools
	import pymol

	pymol.finish_launching()
	pymol.cmd.reinitialize()

	pdbname = pdbpath.split('/')[-1].split('.')[0]

	pymol.cmd.load(pdbpath, pdbname)

	# this draws points that are used to find helix order
	#print('PPO', ppo)
	for p in enumerate(ppo):
		pymol.cmd.pseudoatom('ppo_' + str(p[0]), pos=list(p[1]))

	### axis drawing code ###
	if helices_axis:
		for helix in enumerate(helices_axis):

			h_start  = list(helix[1][0])
			h_end    = list(helix[1][1])
			h_name_s = 'axp' + str(helix[0]) + '_start'
			h_name_e = 'axp' + str(helix[0]) + '_end'

			pymol.cmd.pseudoatom(h_name_s, pos=h_start)
			pymol.cmd.pseudoatom(h_name_e, pos=h_end)

			### this version works for pymol 1.7
			# pymol.cmd.distance('dist' + str(helix[0]), '/' + h_name_s + '/////1', '/' + h_name_e + '/////1')
			### this version works for pymol 2.0
			pymol.cmd.distance('axis' + str(helix[0]), '/' + h_name_s + '///1', '/' + h_name_e + '///1')
	### === ###

	### only show selection and color it code ###
	if color_selection:
		pymol.cmd.hide('cartoon', pdbname)
		pymol.cmd.show('cartoon', ' or '.join(color_selection))
		pymol.cmd.color('red', ' or '.join(color_selection))
	### === ###

	### layer creation code
	for layer in enumerate(layer_points):
		for point in enumerate(layer[1]):
			pymol.cmd.pseudoatom('l_' + str(layer[0]) + 'point' + str(point[0]), pos=point[1])

		for dist in itertools.combinations(range(len(layer[1])), 2):
			# draw only interactions between neighbouring helices (if helix_order not None)
			if (helix_order and (dist[0], dist[1]) not in helix_order):
				continue
			if pymol_version == 2.0:
				pymol.cmd.distance('l_' + str(layer[0]) + 'dist_' + str(dist[0]) + '_' + str(dist[1]), '/l_' + str(layer[0]) + 'point' + str(dist[0]) + '///1', '/l_' + str(layer[0]) + 'point' + str(dist[1]) + '///1')
			else:
				pymol.cmd.distance('l_' + str(layer[0]) + 'dist_' + str(dist[0]) + '_' + str(dist[1]), '/l_' + str(layer[0]) + 'point' + str(dist[0]) + '/////1', '/l_' + str(layer[0]) + 'point' + str(dist[1]) + '/////1')
	### end layer creation


	### alternate layer creation code - less entities in pymol
	# for layer in enumerate(layer_points):
	#     for point in enumerate(layer[1]):
	#         pymol.cmd.pseudoatom('layer' + str(layer[0]) + 'pts', pos=point[1])
	#
	#     #pymol.cmd.distance('layer' + str(layer[0]) + 'dist', '/layer' + str(layer[0]) + 'pts' + '/////1', '/layer' + str(layer[0]) + 'pts' + '/////1')
	#     pymol.cmd.distance('layer' + str(layer[0]) + 'dist', '/layer' + str(layer[0]) + 'pts' + '///1', '/layer' + str(layer[0]) + 'pts' + '///1')

	### end layer creation

	pymol.cmd.set('dash_gap', 0)
	pymol.cmd.set('dash_radius', 0.2)
	pymol.cmd.hide('labels')
	pymol.cmd.hide('lines')
	#pymol.cmd.show('cartoon')
	pymol.cmd.remove('resn HOH')

	print('Saving pymol session with layers...')
	pymol.cmd.save(savepath + pdbname + '_' + suffix + '.pse', format='pse')

	#pymol.cmd.quit()
	print(pdbname + ' saved.')
