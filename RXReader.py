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
def simple_prop_fetch(G,key,prop):
	'''Function to fetch properties from edges or nodes in a NetworkX Graph
	Input: 
	- G. NetworkX graph to query.
	- key. Item of the graph to locate, either string (nodes) or tuple of strings (edge, defined by node names)
	- prop. String, name of the property to be fetched.
	'''
	if (isinstance(key,tuple)):
		return G.edges[key][prop]
	else:
		return G.nodes[key][prop]

def xyz_line_parse(xyzline):
	'''Direct parsing of a raw XYZ-format line to format it as string/float/float/float
	Input:
	- xyzline. String, directly read from file.
	Output:
	- List with x,y,z coordinates formatted as numbers'''
	split_line = xyzline.split()
	at = split_line[0]
	xyz = [float(item) for item in split_line[1:]]
	return [at] + xyz 

def molden_vibration_parser(moldenfile):
	'''Parse a MOLDEN file containing vibrational structures.
	Input:
	- moldenfile. String, filename to be read.
	Output:
	- freqs. List of floats containing frequencies parsed from the file.
	- coords. List of lists containing XYZ coordinates.
	- displacement. List of lists of lists containing XYZ coordinates for each atomic displacement.'''
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
# Using a dictionary allows to change the queried fields for ts, min or prod
query_dict = {"ts":"SELECT energy,zpe,geom,freq,lname FROM ts ",
			"min":"SELECT energy,zpe,geom,freq,lname FROM min ",
			"prod":"SELECT energy,zpe,geom,freq,name,formula FROM prod "}

def query_energy(dbcursor,filterstruc,tablename,add_zpe=False):
	'''Get energy and ZPEs from a SQL table and sum them.
	Input:
	- dbcursor. SQLite3 cursor for a DB connection.
	- filterstruc. String, SQL statement of the type WHERE to filter a given entry
	- tablename. String, name of the table in the SQL file.
	- add_zpe. Boolean, if True, sum zero-point vibrational energy to the energy results
	Output:
	- energies. List of floats containing queried energies. In principle it is expected to have 
	  one-element lists, but it is returned like this for safety.'''
	qtemp = "SELECT energy,zpe FROM %s " % tablename
	if (filterstruc):
		qtemp += filterstruc
	matches = dbcursor.execute(qtemp).fetchall()
	if (add_zpe):
		energies = [sum(match) for match in matches]
	else:
		energies = [match[0] for match in matches]
	return energies

def query_all(dbcursor,filterstruc,tablename,add_zpe=False,multi_query=False):
	'''Get energy and ZPEs from a SQL table and sum them
	Input:
	- dbcursor. SQLite3 cursor for a DB connection.
	- filterstruc. String, SQL statement of the type WHERE to filter a given entry
	- tablename. String, name of the table in the SQL file.
	- add_zpe. Boolean, if True, sum zero-point vibrational energy to the energy results.
	- multi_query. Boolean, if True, several elements are expected to be queried, and will be returned as lists. Else,
	a single element will be returned, corresponding to the first match of the query.
	Output:
	- output. Dictionary containing queried elements (as lists if multi_query == True). Keys:
	(In principle it is expected to have one-element queries, but lists are returned for safety. To keep only)
	+ energy. Float, queried energy. 
	+ geometry. String containing newline-separated blocks for Cartesian coordinates of a molecule.
	+ frequencies. Strings containing newline-separated blocks for all frequency values of a molecule
	+ fnames. String, name for the current molecule.
	+ *formula. For products only, string containing the molecular formula of the structure.
	'''
	# Name field depends on the situation
	qtemp = query_dict[tablename]
	if (filterstruc):
		qtemp += filterstruc
	matches = dbcursor.execute(qtemp).fetchall()
	if (add_zpe):
		energies = [sum(match[0:2]) for match in matches]
	else:
		energies = [match[0] for match in matches]
	geometry = [match[2] for match in matches]
	frequencies = [match[3] for match in matches]
	names = [match[4] for match in matches]
	if (multi_query):
		output = {"energy":energies,"geometry":geometry,"frequencies":frequencies,"fname":names}
	else:
		output = {"energy":energies[0],"geometry":geometry[0],"frequencies":frequencies[0],"fname":names[0]}
	# For products, the formula field shall be included too
	if (tablename == "prod"):
		formulae = [match[5] for match in matches]
		if (multi_query):
			output["formula"] = formulae
		else:
			output["formula"] = formulae[0]
	return output


# Reaction network reading
def RX_parser(workfolder,rxnfile="RXNet"):
	'''Direct parsing of RXNet files generated by AutoMeKin, preparing the general structure (connectivity) to later
	instantiate the network
	Input:
	- workfolder. String, folder where the RXNet file and the SQL databases will be read from
	- rxnfile. String, specific name of the RXNet file to be used
	Output:
	- data. List of lists, containing the tags and indices of the node - edge - node pairs in the reaction
	network required to build the corresponding graph'''
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

def dict_updater(indict,efield,namefield):
	'''Helper function for network building in RX_builder(). Takes an input dict as generated by fetch_all(),
	adapts frequencies (to a semicolon-separated format), replaces the energy by the value in efield (usually relative energy)
	and adds a name value for the corresponding node/edge
	Input:
	- indict. Dictionary, generated by fetch_all() with the multi_query=False switch
	- efield. Float, energy value for the output dict (e.g. the pre-computed relative energy).
	- name. String, name of the node or edge in the graph.
	Output:
	- outdict. Dictionary, analogous to indict but containing reformatted frequencies, modified energy and name field.'''
	nw_freq = indict["frequencies"].replace("\n",";")
	outdict = {"name":namefield,"energy":efield,"frequencies":nw_freq}
	outdict.update({k:v for k,v in indict.items() if k not in ["energy","frequencies"]})
	return outdict

def RX_builder(workfolder,data,orig_selec_mode=False):
	'''From the connectivity data parsed in RX_parser(), generate the graph structure and querying information from the corresponding 
	SQL databases to add information (energy, geometry and lists of frequencies) to the network. Relative energies are also computed here,
	taking the lowest-energy intermediate as reference.
	Input:
	- workfolder. String, folder where the RXNet file and the SQL databases will be read from
	- data. List of lists, containing the tags and indices of the node - edge - node pairs in the reaction
	network required to build the corresponding graph -> RX_parser() direct output.
	- orig_selec_mode. Boolean, if True, remove entries starting from PROD structures following original AMK coding
	Output:
	- G. NetworkX.Graph() object including intermediate and TS information embedded in every node and edge, as fetched
	from the corresponding databases. Additional information used on graph creation (e.g. the reference structure) is passed
	as Graph properties.
	'''
	network_info = {}
	# Establish DB connections and fetch all required information to generate a NetworkX graph
	# Prepare DB connections, as read-only
	dbnames = ["file:" + workfolder + "/%s.db?mode=ro" % entity for entity in ["min","ts","prod"]]
	dbconnections = [sqlite3.connect(db,uri=True) for db in dbnames]
	dbmin,dbts,dbprod = [dbcon.cursor() for dbcon in dbconnections]
	dbdict = {"min":dbmin, "ts":dbts, "prod":dbprod}
	
	# We first need to have a energy reference: for consistency with the original implementation,
	# we will be using the middle column element with the minimum index, which due to ordering has the minimum energy
	smin = min([entry[1] for entry in data[1]])
	e_ref = query_energy(dbmin,"WHERE id==%s" % smin,"min")[0]

	# the units of this energy are kcal/mol from MOPAC (LL) and hartree from G09 (HL)
	if ("LL" in workfolder):
		network_info["e_units"] = "kcal/mol"
	else:
		network_info["e_units"] = "hartree"

	# save in output dict
	network_info.update({"ref_struc":"MIN%d" % smin, "ref_energy":e_ref})
	# Build the graph
	nodelist = []
	edgelist = []
	node_build_dict = {}

	for tag,ndx in zip(data[0],data[1]):
		# Original implementation avoids routes starting with PROD: allow to keep them
		# Prepare queries and analyze the right side which can be PROD or MIN, treated differently
		ts_ii,side1_ii,side2_ii = ndx
		qts,qm1,qm2 = ["WHERE id==%s" % ival for ival in ndx]
		#e_ts,geom_ts,freq_ts,fname_ts = [elem[0] for elem in query_all(dbdict["ts"],qts,"ts")]
		data_ts = query_all(dbdict["ts"],qts,"ts")
		e_ts = data_ts["energy"]

		# use the lowercased tags as selectors for the DBs
		sel_m1,sel_m2 = [sel.lower() for sel in tag[1:]]

		if (sel_m2 == "prod" and orig_selec_mode):
			# Skip cases starting with PROD for consistency IN original implementation
			continue
		data_m1 = query_all(dbdict[sel_m1],qm1,sel_m1)
		e_m1 = data_m1["energy"]
		data_m2 = query_all(dbdict[sel_m2],qm2,sel_m2)
		e_m2 = data_m2["energy"]
		
		#e_m1,geom_m1,freq_m1,fname_m1 = [elem[0] for elem in query_all(dbdict[sel_m1],qm1,sel_m1)]
		#e_m2,geom_m2,freq_m2,fname_m2 = [elem[0] for elem in query_all(dbdict[sel_m2],qm2,sel_m2)]

		# check for self-loops and remove them: CONTINUE if we have same index and same selector
		if ((side1_ii == side2_ii) and (sel_m1 == sel_m2)):
			continue

		# Compute relative energies and generate consistent labels: if we have hartree, convert
		relvals = [e - e_ref for e in [e_ts,e_m1,e_m2]]
		if (network_info["e_units"] == "hartree"):
			relvals = [hartree_kcal*value for value in relvals]
		labels = [name + str(ii) for name,ii in zip(tag,ndx)]
		
		# Construct energy dictionary and a full edge constructor with connections, name and energy parameters
		nn1,nn2 = labels[1:3]
		# Dictionary updaters
		if (nn1 not in node_build_dict.keys()):
			nodelist.append((nn1,dict_updater(data_m1,relvals[1],nn1)))
		if (nn2 not in node_build_dict.keys()):
			nodelist.append((nn2,dict_updater(data_m2,relvals[2],nn2)))
		edgelist.append((nn1,nn2,dict_updater(data_ts,relvals[0],labels[0])))

	# Now generate the graph and then add the corresponding node energies
	G = nx.Graph()
	G.add_edges_from(edgelist)
	G.add_nodes_from(nodelist)
	# Add the network_info dict as a Graph.Graph property
	G.graph.update(network_info)
	return G

def vibr_displ_parser(workfolder,G):
	'''Add vibrational displacements to a existing graph, for further visualization. This information is taken from the
	MOLDEN files stored in the normal_modes/ folder as generated by AutoMeKin. PROD structures do not have this vibrational info.
	Input:
	- workfolder. String, folder where the RXNet file and the SQL databases will be read from
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	Output:
	The G object is modified in-place.
	'''
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
	return None

def graph_plotter(G,posx=None):
	'''Basic MPL plotting for reaction networks generated with RX_builder()
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- posx. Layout for NetworkX.Graph. If None, a shell_layout is assumed.
	Output:
	- fig, ax. Matplotlib figure and axis.'''
	fig = plt.figure()
	ax = fig.add_subplot(111)
	# We need position layout & node energies for coloring
	if (not posx):
		posx = nx.shell_layout(G)
	energvals = [entry[1] for entry in nx.get_node_attributes(G,"energy").items()]
	nx.draw(G,posx,with_labels=False,node_color=energvals,alpha=0.5,ax=ax)
	nx.draw_networkx_labels(G,posx,nx.get_node_attributes(G,"name"),font_weight="bold",ax=ax)
	edgelabs = nx.get_edge_attributes(G,"name")
	nx.draw_networkx_edge_labels(G,posx,edgelabs,font_color="red",ax=ax)
	return fig,ax

def graph_plotter_interact(G,figsize=(10,10),posx=None):
	'''Interactive MPL plotting for reaction networks generated with RX_builder(), in which nodes can be clicked to get 
	energy values, in kcal/mol.
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- figsize. Tuple of floats, size of the resulting plot
	- posx. Layout for NetworkX.Graph. If None, a shell_layout is assumed.
	Output:
	- fig, ax. Matplotlib figure and axis.
	'''

	# Set nested event response functions here so they can freely access all params in the function
	
	def annotator(axis,text,loc):
		'''
		Adds an Annotation artist in a requested location, with a given text.
		'''
		note = matplotlib.text.Annotation(text,loc,backgroundcolor="#FF1493",fontsize=14)
		axis.add_artist(note)
		return note
	
	def onpickx(event):
		'''Event picker, responsive to nodes (PathCollection from plt.scatter) or edges (LineCollection from ArrowPaths).
		Adds annotations for the energy of the nodes/edges that are clicked. If the note does already exist, it is deleted.
		Returns the index of the selection.'''
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
			current_note.remove()
			del note_track[idname]
		ax.figure.canvas.draw_idle()
		return ind
	
	# Unit checking	
	unit_change = G.graph["e_units"] == "hartree"
	note_track = {}
	# Figure generation
	fig = plt.figure(figsize=figsize)
	ax = fig.add_subplot(111)
	ax.set_frame_on(False)
	if (not posx):
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

# Basic profile plotting
def node_position(edge,pos_dict,x0,cntr):
	'''Assignment of node positions for profile generation, inspired by the gnuplot scripts built in AutoMeKin. Sign assignment allows to set 
	nodes at both sides of the plot by starting by a single node, having a more compact plot.
	Input:
	- edge. Tuple of strings defining an edge whose nodes we are plotting.
	- pos_dict. Dictionary mapping node names to the position of the 1st point on the segment, to position structures that 
	  were already defined
	- x0. Float, default initial position for the first node in the edge.
	- cntr. Integer, a counter to handle signs.
	Output:
	- [n1_x0,n1_x1]. List of floats, x positions for the segment of the first node.
	- [n2_x0,n2_x1]. List of floats, x positions for the segment of the second node.
	- sign. Integer, +-1, current value of the sign.
	-'''
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

def line_iterator(G,node,pos_dict,used_ts,x0=0):
	'''For a given node in the graph, detects its inmediate neighbors and prepares the corresponding lines in the graph
	(node A, edge AB and node B, as a XY array to be then passed to matplotlib), for any number of neighbors. Labels
	with names are also generated.
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- node. String, name of the node whose neighbors are being currently detected to draw.
	- pos_dict. Dictionary mapping node names to the position of the 1st point on the segment, to position structures that 
	  were already defined
	- used_ts. List of strings, containing the TSx names for the edges that have already been drawn (as TSs must not be repeated)
	- x0. Float, default initial position for the first node.
	Output:
	- list_arrs. List of XY arrays for each node A/edge AB/node B combination found from 1st neighbor search from /node/.
	- list_labels. List of labels for each node A/edge AB/node B combination.
	- neighborhood. List of strings with the neighbors of the input node.
	'''
	list_arrs = []
	list_labels = []
	# Check positioning
	neighborhood = list(G.neighbors(node))
	neigh_edges = [(node,neigh) for neigh in neighborhood]
	cnt = 0
	for ie,ed in enumerate(neigh_edges):     
		tsname = G.edges[ed]["name"]
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
		n1_e = G.nodes[ed[0]]["energy"]
		n2_e = G.nodes[ed[1]]["energy"]
		ts_e = G.edges[ed]["energy"]
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

def simple_profile_generator(G,startnode):
	'''Applies line_iterator() through the whole network by iterating along the neighbors until no
	connections remain to be read. The graph must be fully connected, as it should in this kind of setup.
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- startnode. String, name of the first node to be analyzed in the recursive neighbor search protocol.
	Output:
	- list_profiles. Full list of XY arrays for each node A/edge AB/node B combination 
	- list_labels. Full list of strings with labels for each node A/edge AB/node B combination 
	 '''

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
	return list_profiles,list_labels

def simple_profile_plotter(arrlist,lablist,put_energy=False,figsize=(10,10)):
	'''Matplotlib plotting of the XY arrays generated by simple_profile_generator()
	Input:
	- arrlist. Full list of XY arrays for each node A/edge AB/node B combination, from simple_profile_generator()
	- lablist. Full list of strings with labels for each node A/edge AB/node B combination, from full_profile_generator().
	- put_energy. Boolean, if True add additional annotations with the energy of each node and edge.
	- figsize. Tuple of floats, figure sizing.
	Output:
	- fig,ax. Matplotlib figure and axis
	'''
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
	'''From a node or edge, extract the "geometry" fiels in a xyz file, indicating name & energy in comment line
	Input:
	- item_select. Node or edge selected from the nx.Graph object, taken via G.nodes[NODE] or G.edges[(NODE_A,NODE_B)]
	Output:
	- xyzblock. String, newline-separated XYZ-compliant block.
	'''
	# get no. of atoms
	geoblock = item_select["geometry"]
	Nat = geoblock.count("\n") + 1
	xyzblock = "%d\n" % Nat
	xyzblock += "%s    energy=%.2f \n" % (item_select["name"],item_select["energy"])
	xyzblock += item_select["geometry"]
	return xyzblock

def profile_geometrist(G,label_list,out_folder="xyz_profiles"):
	'''From the list of labels generated by full_profile(), fetch all geometries and join
	them in XYZ files for the 3-structure reaction stage.
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- lablist. Full list of strings with labels for each node A/edge AB/node B combination, from full:profile().
	- out_folder. String, folder to write XYZ files to.
	Output:
	Generates XYZ files in the requested folder
	'''
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


# Alternative profile generation: generate complete profiles instead of isolated INT/TS/INT triplets
def neighborizer(G,nodein,edgelist):
	'''For a given node nodein in a graph G, obtain the neighbors and the triples
	defining the connection (node - edge - node), without repetition of already
	traversed edges
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- nodein. String, name of the first node to be analyzed in the recursive neighbor search protocol.
	- edgelist. List, containing strings for the edge names (of the form TSn) that have been already traversed
	in previous steps of the iterative search.
	Output:
	edgelist is updated with the edges between nodein and its neighbors. Only entries whose edge is not
	already in edgelist are processed
	- clean_neighbors. List of strings, containing node names for the valid neighbors of nodein 
	- triplets. List of lists, containing NODE1,EDGE,NODE2 strings for the valid neighbors of nodein
	'''
	neighbors = [nd for nd in G.neighbors(nodein)]
	edges = [(nodein,nd) for nd in neighbors]
	edgenames = [G.edges[edg]['name'] for edg in edges]
	# but we only want those neighbors coming from new edges
	clean_ndx = [ie for ie,edg in enumerate(edgenames) if edg not in edgelist]
	clean_neighbors = [neighbors[ie] for ie in clean_ndx]
	clean_edgenames = [edgenames[ie] for ie in clean_ndx]
	edgelist.extend(clean_edgenames)
	triplets = [[nodein,edg,nd] for nd,edg in zip(clean_neighbors,clean_edgenames)]
	return clean_neighbors,triplets

def theor_profile_finder(G,work_nodes,count,track_dict,edge_bookkeep):
	'''
	For a list of lists containing nodes in a graph G, obtain the neighbors via neighborizer()
	and build the profiles, keeping track of the NODE1,EDGE,NODE2 triples generated in previous
	iterations. Shall be done in a loop.
	Ìnput:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- work_nodes. List of lists, containing the nodes to be explored in the current iteration.
	- count. Integer, counter used to keep track of iterations.
	- track_dict. Dictionary, using count as key to store all NODE1,EDGE,NODE2 triplets generated by neighborizer()
	in a pass of the function.
	- edge_bookkeep. List of strings containing TSn names for the previously explored edges in G.
	Output:
	track_dict and edge_bookkeep are updated during the run.
	- next_neighbors. List of lists, containing node names for the valid neighbors of every node in work_nodes.
	- final_profiles. List of lists, containing profiles for which no further exploration is possible and are already finished.
	'''
	final_profiles = []
	next_neighbors = []
	track_dict[count] = []
	for ii,ndlist in enumerate(work_nodes):
		for jj,nd in enumerate(ndlist):
			current_neigh,current_triple = neighborizer(G,nd,edge_bookkeep)
			next_neighbors.append(current_neigh)
			# Keep track of the profiles
			if (count > 0):
				prev = track_dict[count - 1][ii][jj]
				if (current_triple):
					current_prof = [prev + trip[1:] for trip in current_triple]
				else:
					current_prof = prev
					final_profiles.append(current_prof)
				track_dict[count].append(current_prof)
			else:
				track_dict[count].append(current_triple)

	# update orig. list
	work_nodes_flat = [node for sublist in next_neighbors for node in sublist]

	return next_neighbors,final_profiles

def theor_profile_builder(G,start_node,print_flag=False):
	'''Iterative generation of profiles via theor_profile_finder(), starting
	from a given node.
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- start_node. String, name of the node to start profiles from.
	- print_flag. Boolean, if True print state of the iterator.
	Output:
	- final_profiles. List of lists, containing profiles for which no further exploration is possible'''

	# Definition of lists, dicts and counters
	Nedges = len(G.edges)
	edge_bookkeep = []
	ct = 0
	track_dict = {}
	final_profiles = []
	work_nodes = [[start_node]]
	# Keep iterating until no new connections are found
	iterate_flag = True
	while (iterate_flag):
		next_neighbors,ready_profiles = theor_profile_finder(G,work_nodes,ct,track_dict,edge_bookkeep)
		final_profiles.extend(ready_profiles)
		work_nodes = next_neighbors
		if (print_flag):
			print(next_neighbors)
			print("------iteration %d------" % ct) 
		ct += 1
		# Finish checks: no new connections, too many iterations OR all edges explored
		if (not next_neighbors):
			print("No new connections found")
			iterate_flag = False
		if (ct >= 50):
			print("Iterations exceeded: graph is not well defined")
			iterate_flag = False
		elif (len(edge_bookkeep) == Nedges):
			print("All edges have been traversed")
			iterate_flag = False

	# Do an additional iteration to dump the final profile
	next_neighbors,ready_profiles = theor_profile_finder(G,work_nodes,ct,track_dict,edge_bookkeep)
	final_profiles.extend(ready_profiles)
	return final_profiles

def theor_cycle_builder(G,start_node=None):
	'''Generation of cyclic profiles for a given graph: optionally use start_node as starting point of all cycles
	in which it is contained.
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- start_node. String, name of the node which will be considered as cycle start.
	Output.
	- traversed_edges. List, containing names of the transition states that have been already considered in cycle determination.
	- cyclic_profiles. List of lists, containing profiles that comprise the cycles
	'''
	# Determine the cycle basis and reorder it to use start_node as reference where possible
	cycle_basis = nx.cycle_basis(G)
	for ic,cycle in enumerate(cycle_basis):
		if (start_node in cycle):
			ndx = cycle.index(start_node)
			nw_cycle = cycle[ndx:] + cycle[:ndx+1]
			cycle_basis[ic] = nw_cycle
		else:
			cycle += [cycle[0]]
			cycle_basis[ic] = cycle
	
	# Now go along the cycles to check the edges that have been used and to properly define cyclic profiles
	cyclic_profiles = []
	traversed_edges = []
	for cycle in cycle_basis:
		edge_list = [(cycle[ii],cycle[ii+1]) for ii,item in enumerate(cycle) if ii < len(cycle) - 1]
		edge_names = [G.edges[ed]["name"] for ed in edge_list]
		# Generate the profiles from this list and the list of minima (cycle)
		profile_base = [[cycle[ii],edge_names[ii]] for ii in range(len(edge_list))]
		profile = [item for sublist in profile_base for item in sublist] + [cycle[-1]]
		cyclic_profiles.append(profile)
		# Check the edges that have been traversed
		new_edge_names = [ename for ename in edge_names if ename not in traversed_edges]
		traversed_edges.extend(new_edge_names)

	return traversed_edges,cyclic_profiles

def theor_cycle_branch_builder(G,start_node=None):
	'''Generation of cyclic profiles for a given graph via theor_cycle_builder() AND
	addition of branches via theor_profile_builder()
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- start_node. String, name of the node which will be considered as cycle start.
	Output.
	- traversed_edges. List, containing names of the transition states that have been already considered in cycle determination.
	- cyclic_profiles. List of lists, containing profiles that comprise the cycles
	'''
	trav_edges,cycles = theor_cycle_builder(G,start_node)
	# Check whether all edges have already been traversed
	if (len(trav_edges) == len(G.edges)):
		print("All edges traversed: no branching")
		return cycles
	else:
		# Fetch all names
		all_edges = [item for item in G.edges(data="name")]
		new_edges = [item for ndx,item in enumerate(all_edges) if item[2] not in trav_edges]
		# By now, only add these branches as triplets
		branches = [[item[0],item[2],item[1]] for item in new_edges]
		return cycles + branches

def theor_profile_plotter(G,profile_list,figsize=(10,10),cmap="tab10"):
	'''For a list of profiles obtained through theor_profile_builder() or theor_cycle_branch_builder(),
	plot, fetch energy values from the graph G and plot the profile through Matplotlib
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- profile_list. List of lists of strings containing the labels defining the profiles generated by
	theor_profile_builder() or theor_cycle_branch_builder().
	- figsize. Tuple of floats, figure sizing.
	- cmap. Matplotlib colormap to be used.
	Output:
	- fig,ax. Matplotlib figure and axis.
	'''
	# Define colormap: fetch it if it is a string, use as-is elsewhere
	if (isinstance(cmap,str)):
		cmap = plt.get_cmap(cmap)

	# Inner plotting parameters
	xshift,yshift = [0.5,2]
	# Figure definition
	fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(10,10))
	ax.set_frame_on(False)
	ax.axes.get_xaxis().set_visible(False)
	ax.get_yaxis().tick_left()
	ax.axes.grid(b=True,which='both',linestyle='--',c="silver")		

	# System to keep track of labels that have been already been used in a certain position
	# Save xpos and name in a dictionary and check before using the label
	label_tracking = {}

	for ii,prof in enumerate(profile_list):
		# Recover energies for nodes (as strings) and edges (as tuples of node names)
		nodes = prof[::2]
		edges = [(nodes[jj],nodes[jj+1]) for jj in range(len(nodes)-1)]
		Enodes = [simple_prop_fetch(G,nd,"energy") for nd in nodes]
		Eedges = [simple_prop_fetch(G,ed,"energy") for ed in edges]
		# transform back to profile
		prof_energies_raw = [(End,Eed) for End,Eed in zip(Enodes[:-1],Eedges)]
		prof_energies = [item for sublist in prof_energies_raw for item in sublist] + [Enodes[-1]]
		# build XY array from the energies
		xvals = np.arange(0,2*len(prof))
		yvals = [item for yval in prof_energies for item in (yval,yval)]
		xyarr = np.column_stack([xvals,yvals])
		
		# First plot an skeleton for the profile, then add wider lines for all species and the corresponding labels
		ax.plot(xyarr[:,0],xyarr[:,1],color=cmap(ii))
		for jj in range(len(prof)):
			single_line = xyarr[2*jj:2*jj+2]
			ax.plot(single_line[:,0],single_line[:,1],lw=4,c=cmap(ii))
		# prepare labels
		for loc,text in zip(xyarr[::2,:],prof):
			loc[0] += xshift
			# Check whether this x-value was already used to define the label
			prev_xloc = label_tracking.get(text)
			if (prev_xloc == None):
				label_tracking[text] = []
			elif (loc[0] in prev_xloc):
				continue
			# Add to dict and add labels to plot
			label_tracking[text].append(loc[0])	
			ax.annotate(text,loc,horizontalalignment='center',verticalalignment='bottom').draggable()
			evalue = "%.2f" % loc[1]
			ax.annotate(evalue,[loc[0],loc[1]-yshift],horizontalalignment='center',verticalalignment='top').draggable()
	return fig,ax	