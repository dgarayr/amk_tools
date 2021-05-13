'''April 2021, Diego Garay-Ruiz, Institute of Chemical Research of Catalonia
Module to read reaction networks from AutoMeKin, generating a NetworkX.Graph object
which can be then translated into a energy profile.
In contrast with the original implementation, here everything is read from files in
the FINAL_LL_/FINAL_HL_ folders, which have self-consistent indexing.
To keep it straightforward, only RXNet file & min.db, ts.db and prod.db files shall be read'''
import numpy as np
import sqlite3
import matplotlib.pyplot as plt
import matplotlib.collections as collections
import matplotlib.text
import networkx as nx

# Constant definitions (from scipy.constants, redefined to avoid additional dependencies)
hartree_J = 4.3597447222071e-18
Navogadro = 6.02214076e+23
hartree_kcal = hartree_J*Navogadro/4184

# Helper functions
def dict_adder(indict,key,value):
	# Check if an entry is in the dictionary before adding it
	if (key not in indict.keys()):
		indict[key] = value
	return None

def xyz_line_parse(xyzline):
	split_line = xyzline.split()
	at = split_line[0]
	xyz = [float(item) for item in split_line[1:]]
	return [at] + xyz 

def molden_vibration_parser(moldenfile):
	'''Parse a MOLDEN file containing vibrational structures'''
	with open(moldenfile,"r") as fmold:
		dump = fmold.read()
	sel_list = dump.split("]")[2:]
	proc_list = [item.split('[')[0].strip() for item in sel_list]
	# Now we can process every item
	freqs = [float(entry) for entry in proc_list[0].split()]
	coords = [xyz_line_parse(line) for line in proc_list[1].split("\n")]
	# Get no. of atoms from coordinates and use this to handle vibrations as a list of lists
	displ_block = proc_list[2].split("\n")
	nlines = len(coords) + 1
	displacements = [displ_block[ii+1:ii+nlines] for ii in range(0,len(displ_block),nlines)]
	for jj,displ in enumerate(displacements):
		displacements[jj] = [line.split() for line in displ]
	return freqs,coords,displacements

# Querying functions
def query_energy(dbcursor,filterstruc,tablename,add_zpe=False):
	'''Get energy and ZPEs from a SQL table and sum them'''
	qtemp = "SELECT energy,zpe FROM %s " % tablename
	if (filterstruc):
		qtemp += filterstruc
	matches = dbcursor.execute(qtemp).fetchall()
	if (add_zpe):
		energies = [sum(match) for match in matches]
	else:
		energies = [match[0] for match in matches]
	return energies

def query_all(dbcursor,filterstruc,tablename,add_zpe=False):
	'''Get energy and ZPEs from a SQL table and sum them
	Input:
	- dbcursor. sqlite3 cursor to a connected database
	- filterstruc. String, WHERE statement to limit the query
	- tablename. Name of the table to be accessed: min, ts or prod
	- add_zpe. Boolean, if True the energies will include zero-point 
	vibrational energy'''
	qtemp = "SELECT energy,zpe,geom,freq FROM %s " % tablename
	if (filterstruc):
		qtemp += filterstruc
	matches = dbcursor.execute(qtemp).fetchall()
	if (add_zpe):
		energies = [sum(match[0:2]) for match in matches]
	else:
		energies = [match[0] for match in matches]
	geometry = [match[2] for match in matches]
	frequencies = [match[3] for match in matches]
	return energies,geometry,frequencies

# Reaction network reading
def RX_parser(workfolder,rxnfile="RXNet"):
	# Handle the main file
	froute = "%s/%s" % (workfolder,rxnfile)
	with open(froute,"r") as frxn:
		# Structure of the file: TS, DE value and rxn path
		# Four possible path heuristics:
		# MIN x <--> MIN y
		# MIN x <--> PRa: A + B
		# PRa: A + B <--> MIN x
		# PRa: A + B <--> PRc: C + D
		tags = []
		indices = []
		dump = [item.strip() for item in frxn.readlines()[2:]] # skipping header
		# use --- to split the reaction path separator
		dump_proc = [item.split("---") for item in dump]
		for line in dump_proc:
			leftside,rightside = [item.strip("<>-").split() for item in line]
			ts = int(leftside[0])
			# Check whether we have a product or a minimum
			if (leftside[2] == "MIN"):
				m1 = int(leftside[3])
				t1 = "MIN"
			else:
				# We have PRa: and we have to extract the integer a
				m1 = int(leftside[2].strip(":").replace("PR",""))
				t1 = "PROD"
			# do the same for the other side
			if (rightside[0] == "MIN"):
				m2 = int(rightside[1])
				t2 = "MIN"
			else:
				m2 = int(rightside[0].strip(":").replace("PR",""))
				t2 = "PROD"
			tags.append(["TS",t1,t2])
			indices.append([ts,m1,m2])
		# fetch tags MIN/TS/PROD and the corresponding indices for each 
		data = [tags,indices]
	return data

def RX_builder(workfolder,data):
	network_info = {}
	# Establish DB connections and fetch all required information to generate a NetworkX graph
	# Prepare DB connections, as read-only
	dbnames = ["file:" + workfolder + "/%s.db?mode=ro" % entity for entity in ["min","ts","prod"]]
	dbconnections = [sqlite3.connect(db,uri=True) for db in dbnames]
	dbmin,dbts,dbprod = [dbcon.cursor() for dbcon in dbconnections]
	
	# We first need to have a energy reference: for consistency with the original implementation,
	# we will be using the middle column element with the minimum index, which due to ordering has the minimum energy
	smin = min([entry[1] for entry in data[1]])
	e_ref = query_energy(dbmin,"WHERE id==%s" % smin,"min")[0]

	# the units of this energy are kcal/mol from MOPAC (LL) and hartree from G09 (HL)
	if ("LL" in workfolder):
		network_info["e_units"] = "kcal"
	else:
		network_info["e_units"] = "hartree"

	# save in output dict
	network_info.update({"ref_struc":"MIN%d" % smin, "ref_energy":e_ref})
	# Build the graph
	nodelist = []
	edgelist = []
	node_build_dict = {}
	for tag,ndx in zip(data[0],data[1]):
		# Skip cases starting with PROD for consistency
		if (tag[1] == "PROD"):
			continue
		# Prepare queries and analyze the right side which can be PROD or MIN, treated differently
		ts_ii,side1_ii,side2_ii = ndx
		qts,qm1,qm2 = ["WHERE id==%s" % ival for ival in ndx]
		e_ts,geom_ts,freq_ts = [elem[0] for elem in query_all(dbts,qts,"ts")]
		e_m1,geom_m1,freq_m1 = [elem[0] for elem in query_all(dbmin,qm1,"min")]
		if (tag[2] == "PROD"):
			e_m2,geom_m2,freq_m2 = [elem[0] for elem in query_all(dbprod,qm2,"prod")]
		else:
			# check for self-loops from MINx to MINx and remove them
			if (side1_ii == side2_ii):
				continue
			e_m2,geom_m2,freq_m2 = [elem[0] for elem in query_all(dbmin,qm2,"min")]

		# Compute relative energies and generate consistent labels
		relvals = [e - e_ref for e in [e_ts,e_m1,e_m2]]
		labels = [name + str(ii) for name,ii in zip(tag,ndx)]
		
		# Construct energy dictionary and a full edge constructor with connections, name and energy parameters
		nn1,nn2 = labels[1:3]
		# Dictionary updaters
		if (nn1 not in node_build_dict.keys()):
			nodelist.append((nn1,{"name":nn1,"energy":relvals[1],"geometry":geom_m1,"frequencies":freq_m1}))
		if (nn2 not in node_build_dict.keys()):
			nodelist.append((nn2,{"name":nn2,"energy":relvals[2],"geometry":geom_m2,"frequencies":freq_m2}))
		edgelist.append((nn1,nn2,{"name":labels[0],"energy":relvals[0],"geometry":geom_ts,"frequencies":freq_ts}))

	# Now generate the graph and then add the corresponding node energies
	G = nx.Graph()
	G.add_edges_from(edgelist)
	G.add_nodes_from(nodelist)
	# Add the network_info dict as a Graph.Graph property
	G.graph.update(network_info)
	return G,network_info

def vibr_displ_parser(workfolder,G):
	'''Add vibrational displacements to a existing graph, for further visualization'''
	nm_folder = workfolder + "/normal_modes"

	# Nodes first, then edges, remembering that there is no info for PROD species
	for nd in G.nodes(data=True):
		if ("MIN" in nd[0]):
			ndid = int(nd[0].replace("MIN",""))
			nm_file = nm_folder + "/MIN%04d.molden" % ndid
			frq,coords,displ = molden_vibration_parser(nm_file)
			nd[1]["vibr_displace"] = displ

	for ed in G.edges(data=True):
		tsid = int(ed[2]["name"].replace("TS",""))
		nm_file = nm_folder + "/TS%04d.molden" % tsid
		frq,coords,displ = molden_vibration_parser(nm_file)
		ed[2]["vibr_displace"] = displ
	# And the graph is modified in-place
	return None

def graph_plotter(G):
	# Do a basic plot for the graph generated from the RXNet and the databases
	fig = plt.figure()
	ax = fig.add_subplot(111)
	# We need position layout & node energies for coloring
	posx = nx.shell_layout(G)
	energvals = [entry[1] for entry in nx.get_node_attributes(G,"energy").items()]
	nx.draw(G,posx,with_labels=False,node_color=energvals,alpha=0.5,ax=ax)
	nx.draw_networkx_labels(G,posx,nx.get_node_attributes(G,"name"),font_weight="bold",ax=ax)
	edgelabs = nx.get_edge_attributes(G,"name")
	nx.draw_networkx_edge_labels(G,posx,edgelabs,font_color="red",ax=ax)
	return fig,ax

def graph_plotter_interact(G,figsize=(10,10)):
	'''Generates an interactive plot in which nodes can be clicked to get more information:
	by now only energy in kcal/mol'''
	# Set nested event response functions so they can freely access all params in the function
	# check units
	unit_change = G.graph["e_units"] == "hartree"

	note_track = {}
	
	def annotator(axis,text,loc):
		note = matplotlib.text.Annotation(text,loc,backgroundcolor="#FF1493",fontsize=14)
		axis.add_artist(note)
		return note
	
	def onpickx(event):
		# Recall that what we have is a PathCollection (nodes, plt.scatter) or a LineCollection (edges, ArrowPaths)
		if isinstance(event.artist,collections.PathCollection):
			current_artist = event.artist
			element_type = "NODE"
		elif isinstance(event.artist,collections.LineCollection):
			current_artist = event.artist
			element_type = "EDGE"
		ind = event.ind[0]
		# Change reactivity for nodes or edges
		if (element_type == "NODE"):
			xx,yy = current_artist.get_offsets()[ind]
			idname = map_nodes[ind]
			
		elif (element_type == "EDGE"):
			v0,v1 = current_artist.get_paths()[ind].vertices
			# we want the middle point to place the label
			xx,yy = 0.5*(np.array(v0) + np.array(v1))
			idname = map_edges[ind]
			
		if (idname not in note_track.keys()):
			if (element_type == "NODE"):
				energy,objname = [G.nodes[idname][prop] for prop in ["energy","name"]]
			elif (element_type == "EDGE"):
				energy,objname = [G.edges[idname[0],idname[1]][prop] for prop in ["energy","name"]]
			if (unit_change):
				energy = energy*hartree_kcal
			estring = "%.2f" % energy
			notelabel = "%s \n (%.2f)" % (objname,energy)
			noteobj = annotator(ax,notelabel,(xx,yy))
			note_track[idname] = noteobj
		else:
			current_note = note_track[idname]
			#ax.remove_artist(current_note)
			current_note.remove()
			del note_track[idname]
		ax.figure.canvas.draw_idle()
		return ind
	fig = plt.figure(figsize=figsize)
	ax = fig.add_subplot(111)
	ax.set_frame_on(False)
	posx = nx.shell_layout(G)
	# map indices to keys for nodes and edges (nodes are taken from posx)
	map_nodes = {ii:name for ii,name in enumerate(posx.keys())}
	map_edges = {ie:ed for ie,ed in enumerate(G.edges)}
	# fetch energies
	energvals = [entry[1] for entry in nx.get_node_attributes(G,"energy").items()]
	nodecol = nx.draw_networkx_nodes(G,posx,node_size=500,ax=ax,alpha=0.7,node_color=energvals)
	nx.draw_networkx_labels(G,posx,font_size=16,font_weight="bold")
	edgecol = nx.draw_networkx_edges(G,posx)
	nx.draw_networkx_edge_labels(G,posx,nx.get_edge_attributes(G,"name"),font_size=16,font_color="red")
	# Make the plot reactive
	nodecol.set_picker(True)
	edgecol.set_picker(True)
	fig.canvas.mpl_connect('pick_event',onpickx)

	# axis rescaling avoids border clipping for nodes
	rescaling = 1.2
	ax.set_xlim([rescaling * x for x in ax.get_xlim()])
	ax.set_ylim([rescaling * y for y in ax.get_ylim()])
	return fig,ax

# Profile generation functions
def node_position(edge,pos_dict,x0,cntr):
	# Check whether the positions are already assigned
	n1,n2 = edge
	try:
		n1_x0 = pos_dict[n1]
	except:
		n1_x0 = x0
	# In the second node, only keep track for minima
	try:
		n2_x0 = pos_dict[n2]
		sign = np.sign(n2_x0 - n1_x0)
	except:
		# decide the sign according to the counter
		parity = cntr % 2
		if (parity == 0):
			sign = 1
		else:
			sign = -1
		n2_x0 = n1_x0 + sign*8
	# Check that positions don't overlap: if they do, move the second but don't track it
	if (n1_x0 == n2_x0):
		print("Overlapping positions:",n1_e0,n2_x0,edge)
	# now provide all X-positions: n1_x0, n1_x1, n2_x0, n2_x1 and the sign
	n1_x1,n2_x1 = [n_x0 + 2 for n_x0 in [n1_x0,n2_x0]]
	return [n1_x0,n1_x1],[n2_x0,n2_x1],sign

def line_iterator(graph,node,pos_dict,used_ts,x0=0):
	list_arrs = []
	list_labels = []
	# Check positioning
	neighborhood = list(graph.neighbors(node))
	neigh_edges = [(node,neigh) for neigh in neighborhood]
	cnt = 0
	for ie,ed in enumerate(neigh_edges):     
		tsname = graph.edges[ed]["name"]
		if (tsname in used_ts):
			continue
		used_ts.append(tsname)
		# horizontal positioning 
		n1pos,n2pos,sgn = node_position(ed,pos_dict,x0,cnt)
		tspos = [n1x+sgn*4 for n1x in n1pos]
		# Don't keep track of products: these can be repeated along the plot
		for nd,npos in zip(ed,[n1pos,n2pos]):
			if ("MIN" in nd):
				pos_dict[nd] = npos[0]
		# vertical positioning: energies
		n1_e = graph.nodes[ed[0]]["energy"]
		n2_e = graph.nodes[ed[1]]["energy"]
		ts_e = graph.edges[ed]["energy"]
		# prepare XY points & sort them by the X column for consistent plot generation
		x_vals = n1pos + tspos + n2pos
		y_vals = [n1_e]*2 + [ts_e]*2 + [n2_e]*2
		xyarray = np.array([x_vals,y_vals]).T
		mask = xyarray[:,0].argsort()
		xyarray = xyarray[mask]
		list_arrs.append(xyarray)
		# manage label list, /w duplicities, & reorder it with the same mask as before
		lablist = [ed[0]]*2 + [tsname]*2 + [ed[1]]*2
		lablist = np.array(lablist)[mask]
		list_labels += [lablist]
		cnt += 1
	return list_arrs,list_labels,neighborhood

def full_profile(G,startnode):
	# Here we can iterate along the whole network
	# define the intermediate lists & dicts
	pos_dict = {}
	tslist = []
	x0 = 0
	next_gen = []
	list_profiles,list_labels,next_nodes = line_iterator(G,startnode,pos_dict,tslist,x0)	
	update_flag = bool(list_profiles)
	while (update_flag):
		current_state = []
		nw_profiles = []
		nw_labels = []
		for nwstart in next_nodes:
			add_profiles,add_labels,nw_nodes = line_iterator(G,nwstart,pos_dict,tslist,x0)
			next_gen += [item for item in nw_nodes if item not in next_gen]
			nw_profiles += add_profiles
			nw_labels += add_labels
			current_state.append(bool(add_profiles))
		list_profiles += [prof for prof in nw_profiles]
		list_labels += nw_labels
		update_flag = np.any(current_state)    
		next_nodes = next_gen
	# Unit handling
	unit_change = G.graph["e_units"] == "hartree"
	if (unit_change):
		for prof in list_profiles:
			prof[:,1] *= hartree_kcal
	return list_profiles,list_labels

def profplotter(arrlist,lablist,put_energy=False,figsize=(10,10)):
	# Pure plotting function
	# shifts for plotting
	xshift = 1
	yshift = 2
	cmap = plt.get_cmap("tab10")
	fig,ax = plt.subplots(nrows=1,ncols=1,figsize=figsize)
	ax.set_frame_on(False)
	ax.axes.get_xaxis().set_visible(False)
	ax.get_yaxis().tick_left()
	ax.axes.grid(b=True,which='both',linestyle='--',c="silver")
	for ii,(arr,labs) in enumerate(zip(arrlist,lablist)):
		ax.plot(arr[:,0],arr[:,1],"--",c=cmap(ii))
		# wider lines?
		for jj in range(3):
			single_line = arr[2*jj:2*jj+2]
			ax.plot(single_line[:,0],single_line[:,1],lw=4,c=cmap(ii))
		# prepare labels
		sel_labs = labs[::2]
		sel_locs = arr[::2,:]
		for loc,text in zip(sel_locs,sel_labs):
			loc[0] += xshift
			ax.annotate(text,loc,horizontalalignment='center',verticalalignment='bottom')
			evalue = "%.2f" % loc[1]
			if (put_energy):
				ax.annotate(evalue,[loc[0],loc[1]-yshift],horizontalalignment='center',verticalalignment='top').draggable()
	return fig,ax

# Functions to fetch geometries
def xyz_generator(item_select):
	# from a node or edge, extract geometry in a xyz file, indicating name & energy in comment line
	# get no. of atoms
	geoblock = item_select["geometry"]
	Nat = geoblock.count("\n") + 1
	xyzblock = "%d\n" % Nat
	xyzblock += "%s    energy=%.2f \n" % (item_select["name"],item_select["energy"])
	xyzblock += item_select["geometry"]
	return xyzblock

def profile_geometrist(G,label_list,out_folder="xyz_profiles"):
	'''From the list of labels generated by full_profile(), fetch all geometries and join
	them in XYZ files for the 3-structure reaction stage'''
	for lab in label_list:
		px = lab[::2] # as these are duplicated
		out_name = "_".join(px)
		sel_n1 = G.nodes[px[0]]
		sel_n2 = G.nodes[px[2]]
		sel_ts = G.edges[(px[0],px[2])]
		join_xyz = "\n".join([xyz_generator(elem) for elem in [sel_n1,sel_ts,sel_n2]])
		froute = "%s/%s.xyz" % (out_folder,out_name)
		with open(froute,"w+") as fxyz:
			fxyz.write(join_xyz)
	return None