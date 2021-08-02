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
import itertools
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

def query_energy(dbcursor,filterstruc,tablename,add_zpe=True):
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

def query_all(dbcursor,filterstruc,tablename,add_zpe=True,e_units="kcal/mol",multi_query=False):
	'''Get energy and ZPEs from a SQL table and sum them
	Input:
	- dbcursor. SQLite3 cursor for a DB connection.
	- filterstruc. String, SQL statement of the type WHERE to filter a given entry
	- tablename. String, name of the table in the SQL file.
	- add_zpe. Boolean, if True, sum zero-point vibrational energy to the energy results.
	- e_units. String, units of the "energy" field, either "kcal/mol" or "hartree".
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
	# extract energies and handle units here
	e_elec = np.array([match[0] for match in matches])
	if (e_units == "hartree"):
		e_elec *= hartree_kcal
	if (add_zpe):
		# ZPE is always in kcal/mol. If energies are in hartree, convert them.
		zpe = np.array([match[1] for match in matches])
		energies = list(e_elec + zpe)
	else:
		energies = list(e_elec)

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
def RX_parser(workfolder,rxnfile="RXNet",check_connection=False):
	'''Direct parsing of RXNet files generated by AutoMeKin, preparing the general structure (connectivity) to later
	instantiate the network
	Input:
	- workfolder. String, folder where the RXNet file and the SQL databases will be read from
	- rxnfile. String, specific name of the RXNet file to be used
	- check_connection. Boolean, if True, remove entries labeled as DISCONN in the RXNet file (if present)
	Output:
	- data. List of lists, containing the tags and indices of the node - edge - node pairs in the reaction
	network required to build the corresponding graph'''
	froute = "%s/%s" % (workfolder,rxnfile)
	# Flag barrierless networks so we do not wrongly assign a TS tag
	# Use negative integers for the index and TSb as label instead of TS
	barrierless = "barrless" in rxnfile
	with open(froute,"r") as frxn:
		tstag = "TS"    
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
		
		for ii,line in enumerate(dump_proc):
			# Evaluate if the line is connected
			if (check_connection):
				conn_status = dump[ii].split()[-1]
				if (conn_status not in ["CONN","DISCONN"]):
					print("%s does not contain information about connectivity of the network" % rxnfile)
					check_connection = False
				elif (conn_status == "DISCONN"):
					continue
				else:
					print(ii,conn_status)

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
			if (barrierless):
				tstag = "TSb"
				ts = -ts
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

def check_mininfo(workfolder):
	'''Take the MINinfo file in workfolder and read it until the species with energy ZERO is reached, to fetch the
	energy reference in AMK.
	Input:
	- workfolder. String, folder to read the MINinfo file from
	Output:
	- min. Integer, index of the reference minimum'''
	# Use a try/except block in case there is no MINinfo file
	try:
		with open("%s/MINinfo" % workfolder,"r") as fmins:
			next(fmins)			# skip 1st line
			for line_raw in fmins:
				line = line_raw.strip().split()
				e = float(line[1])
				if (abs(e) < 1E-8):
					minstring = line[0].strip()
					smin = float(minstring.replace("MIN",""))
					return smin
	except:
		return None

def RX_builder(workfolder,data,orig_selec_mode=False,add_zpe=True):
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

	# the units for the energy are kcal/mol from MOPAC (LL) and hartree from G09 (HL)
	if ("LL" in workfolder):
		e_units = "kcal/mol"
	else:
		e_units = "hartree"

	# as query_all() handles all energy conversions, units in the graph will always be kcal/mol
	network_info["e_units"] = "kcal/mol"
	# We first need to have a energy reference: for consistency with the original implementation,
	# Check MINinfo and get the element whose energy is ZERO (energy reference consistent /w AMK)
	# If MINinfo is not available, select the MIN with lowest index as fallback
	smin = check_mininfo(workfolder)
	if (not smin):
		smin = min([entry[1] for entry in data[1]])
	e_ref = query_all(dbmin,"WHERE id==%s" % smin,"min",add_zpe=add_zpe,e_units=e_units)["energy"]

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
		# In the queries, check for POSITIVE indices, else we will be in a barrierless case
		qts,qm1,qm2 = ["WHERE id==%s" % ival if ival >= 0 else None for ival in ndx]
		if (qts):
			data_ts = query_all(dbdict["ts"],qts,"ts",add_zpe=add_zpe,e_units=e_units)
			e_ts = data_ts["energy"]
			barrierless = False
		else:
			# Check barrierless transition states so we do not assign anything to this edge apart from the name
			barrierless = True
			data_ts = {"name":"TSb_%d" % abs(ndx[0])}
		data_ts["barrless"] = barrierless
		# use the lowercased tags as selectors for the DBs
		sel_m1,sel_m2 = [sel.lower() for sel in tag[1:]]

		if (sel_m2 == "prod" and orig_selec_mode):
			# Skip cases starting with PROD for consistency IN original implementation
			continue

		data_m1 = query_all(dbdict[sel_m1],qm1,sel_m1,add_zpe=add_zpe,e_units=e_units)
		e_m1 = data_m1["energy"]
		data_m2 = query_all(dbdict[sel_m2],qm2,sel_m2,add_zpe=add_zpe,e_units=e_units)
		e_m2 = data_m2["energy"]

		# assign m2 energy to barrierless TSs
		if (barrierless):
			e_ts = max(e_m1,e_m2)
			data_ts["energy"] = max(e_m1,e_m2)
			data_ts["frequencies"] = None
			data_ts["geometry"] = None
		# check for self-loops and remove them: CONTINUE if we have same index and same selector
		if ((side1_ii == side2_ii) and (sel_m1 == sel_m2)):
			continue

		# Compute relative energies and generate consistent labels: unit handling is done at query_all()
		relvals = [e - e_ref for e in [e_ts,e_m1,e_m2]]
		labels = [name + str(ii) for name,ii in zip(tag,ndx)]
		
		# Construct energy dictionary and a full edge constructor with connections, name and energy parameters
		nn1,nn2 = labels[1:3]
		# Dictionary updaters
		if (nn1 not in node_build_dict.keys()):
			nodelist.append((nn1,dict_updater(data_m1,relvals[1],nn1)))
		if (nn2 not in node_build_dict.keys()):
			nodelist.append((nn2,dict_updater(data_m2,relvals[2],nn2)))
		# For edges, consider possible barrierless cases
		if (barrierless):
			data_ts["energy"] = relvals[0]
			edgelist.append((nn1,nn2,data_ts))
		else:
			edgelist.append((nn1,nn2,dict_updater(data_ts,relvals[0],labels[0])))

	# Now generate the graph and then add the corresponding node energies
	G = nx.Graph()
	G.add_edges_from(edgelist)
	G.add_nodes_from(nodelist)
	# Add the network_info dict as a Graph.Graph property
	G.graph.update(network_info)
	return G

def vibr_displ_parser(workfolder,G,Nmodes=-1):
	'''Add vibrational displacements to a existing graph, for further visualization. This information is taken from the
	MOLDEN files stored in the normal_modes/ folder as generated by AutoMeKin. PROD structures do not have this vibrational info.
	Input:
	- workfolder. String, folder where the RXNet file and the SQL databases will be read from
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- Nmodes. Integer, no. of vibr. modes to extract. If -1, ALL modes are fetched.
	Output:
	The G object is modified in-place.
	'''
	nm_folder = workfolder + "/normal_modes"

	# Process Nmodes variable for slicing: if None, lists will be sliced until the end
	if (Nmodes == -1 or not isinstance(Nmodes,int)): 
		Nmodes = None

	# Nodes first, then edges, remembering that there is no info for PROD species
	for nd in G.nodes(data=True):
		if ("MIN" in nd[0]):
			ndid = int(nd[0].replace("MIN",""))
			nm_file = nm_folder + "/MIN%04d.molden" % ndid
			frq,coords,displ = molden_vibration_parser(nm_file)
			frq_block = [str(frqval) for frqval in frq[0:Nmodes]]
			nd[1]["frequencies"] = ";".join(frq_block)
			nd[1]["vibr_displace"] = displ[0:Nmodes]

	for ed in G.edges(data=True):
		# Only proceed when there is data: skip edges without information such as barrierless TSs
		if (ed[2]["geometry"]):
			tsid = int(ed[2]["name"].replace("TS",""))
			nm_file = nm_folder + "/TS%04d.molden" % tsid
			frq,coords,displ = molden_vibration_parser(nm_file)
			frq_block = [str(frqval) for frqval in frq[0:Nmodes]]
			ed[2]["frequencies"] = ";".join(frq_block)
			ed[2]["vibr_displace"] = displ[0:Nmodes]
		else:
			ed[2]["vibr_displace"] = None
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

# Location of paths

def formula_dict_constructor(G):
	'''Generate a mapping between possible formulas and product(s) in the graph, and store as graph property
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	Output:
	- formula_map. Dictionary, mapping a formula-defining string to the corresponding list of product names
	- None, adds the "formula_map" property to the input graph
	'''
	# Select products and get the formulas
	prod_list = [nd for nd in G.nodes(data=True) if "PROD" in nd[0]]
	formulae = [nd[1]["formula"].replace(" ","") for nd in prod_list]
	prod_names = [nd[0] for nd in prod_list]
	formula_map = {}
	# Take care: a formula may map to several products
	for name,form in zip(prod_names,formulae):
		if (form not in formula_map.keys()):
			formula_map[form] = [name]
		else:
			formula_map[form].append(name)
	G.graph["formula_map"] = formula_map
	return formula_map


def formula_locator(G,formula):
	'''Helper function to fetch the PRODXX code(s) corresponding to a given molecular formula
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- formula. String, defining the A+B+C... molecular formula to be fetched
	'''
	# We should have a dictionary mapping formulas to names as a graph property. Fetch it, and generate if not present
	try:
		formula_map = G.graph["formula_map"]
	except:
		formula_map = formula_dict_constructor(G)
	formula_format = formula.replace(" ","")
	prod_list = formula_map[formula_format]
	return prod_list

def product_collapser(G):
	'''Contract all equivalent products in the graph to a single node. Check energies to select the most stable product.
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	Output:
	- Gwork. NetworkX.Graph() where all equivalent products have been collapsed'''
	
	Gwork = G.copy()

	# Regenerate formula mapping to avoid problems
	formula_map = formula_dict_constructor(Gwork)

	collapsed_dict = {}
	for formula,prod_list in formula_map.items():
		if (len(prod_list) == 1):
			continue
		# Get energies and sort the list to select the lowest-energy product complex as parent
		prod_energies = [G.nodes[prod]["energy"] for prod in prod_list]
		prod_list_sort = [prod_list[ndx] for ndx in np.argsort(prod_energies)]
		parent = prod_list_sort[0]
		[nx.contracted_nodes(Gwork,parent,child,copy=False) for child in prod_list_sort[1:]]
		# and relabel using the formula!
		label = "PROD_" + formula
		nx.relabel_nodes(Gwork,mapping={parent:label},copy=False)
		# prepare the dict mapping the nodes to the new label
		collapsed_dict.update({prod:label for prod in prod_list_sort})
	# If paths are present, modify them accordingly
	if (Gwork.graph["path_list"]):
		collapsed_paths = []
		for path in Gwork.graph["path_list"]:
			nw_path = [collapsed_dict[item] if item in collapsed_dict.keys() else item for item in path]
			collapsed_paths.append(nw_path)
		Gwork.graph["path_list"] = collapsed_paths
	return Gwork

def node_synonym_search(G,nodelist):
	'''Standardize a given list of node elements, containing: i) single name string, ii) list of name strings or iii) a product molecular formula,
	to obtain a list of lists of name strings. Product formulas are fetched in the graph to obtain all the PRODxxx strings matching the formula
	Input
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- nodelist. List containing i) single name strings, ii) list of name strings or iii) product molecular formulas.
	Output:
	- standard_nodes. List of lists of strings corresponding to the input nodes.'''
	standard_nodes = []
	for item in nodelist:
		if (isinstance(item,list)):
			nodes = [formula_locator(G,entry) if "+" in entry else entry for entry in item]
			# and we want to flatten
			flat_nodes = [element for entry in nodes for element in entry]
			nodes = flat_nodes
		else:
			if ("+" in item):
				nodes = formula_locator(G,item)
			else:
				nodes = [item]
		standard_nodes.append(nodes)
	return standard_nodes

def single_path_finder(G,source,target,cutoff,skip_int_frags=True):
	''' Finds all simple paths connecting two nodes, allowing to skip all paths that have fragmented
	PROD structures in the middle  .
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- source. String, name of the node where paths will start.
	- target. String, name of the node where paths will end.
	- cutoff. Integer, maximum depth for the search.
	- skip_int_frags. Boolean, if True, discard paths that have fragmented structures between source and target.
	Output:
	- valid_path_list. List of lists, containing the fetched paths as lists of node names.
	'''
	valid_path_list = []
	for path in nx.all_simple_paths(G,source=source,target=target,cutoff=cutoff):
		# Exclude paths that go through products: we only want these that only pass through minima
		interm_elements = path[1:-1]
		prod_loc = ["PROD" in elem for elem in interm_elements]
		if (np.any(prod_loc) and skip_int_frags):
			continue
		valid_path_list.append(path)
	return valid_path_list

def poly_path_finder(G,source_collection,target_collection,cutoff,skip_int_frags=True):
	'''Finds all simple paths beginning in a set of source nodes and ending in a set of target nodes. Allows
	to skip paths that have fragmented PROD structures in the middle.
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- source. String, name of the node where paths will start.
	- target. String, name of the node where paths will end.
	- cutoff. Integer, maximum depth for the search
	- skip_int_frags. Boolean, if True, discard paths that have fragmented structures between source and target.
	Output:
	- all_paths. List of lists, containing the fetched paths as lists of node names.
	'''
	all_paths = []
	# Get all combinations, via the Cartesian product implemented in itertools
	combination_list = list(itertools.product(source_collection,target_collection))
	all_paths = []
	for comb in combination_list:
		current_paths = single_path_finder(G,source=comb[0],target=comb[1],cutoff=cutoff,skip_int_frags=skip_int_frags)
		all_paths.extend(current_paths)
	return all_paths

def path_reformatter(G,path_list):
	'''Add transition states to a path specified as a list of nodes, for compatibility with profile representation
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- path_list. List of lists, containing the fetched paths from single_path_finder() or poly_path_finder()
	Output:
	- format_path_list. List of lists, containing paths with transition state names in between nodes'''
	format_path_list = []
	for path in path_list:
		edges = [(path[ii],path[ii+1]) for ii in range(len(path) - 1)]
		edge_labels = [G.edges[ed]["name"] for ed in edges]
		format_path_raw = [[path[jj],edge_labels[jj]] for jj in range(len(edges))] 
		format_path = [item for sublist in format_path_raw for item in sublist] + [path[-1]]
		format_path_list.append(format_path)
	return format_path_list

def path_filter(G,path_list,threshold):
	'''For a given list of paths, filter out all paths where any TS is above a given energy threshold
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- path_list. List of lists, containing the fetched paths from single_path_finder() or poly_path_finder()
	- threshold. Float, maximum energy to be accepted.
	Output:
	- valid_paths. List of lists, containing paths that fulfill the input condition'''
	valid_paths = []
	for path in path_list:
		edges = [(path[ii],path[ii+1]) for ii in range(len(path) - 1)]
		edge_energies = np.array([G.edges[ed]["energy"] for ed in edges])
		discard_path = np.any(edge_energies >= threshold)
		if (not discard_path):
			valid_paths.append(path)
	return valid_paths

def graph_path_selector(G,path_list):
	'''Generates a subgraph containing the nodes pertaining to the paths specified in path_list.
	Input:
	- G. NetworkX.Graph() as generated by RX_builder() -> connectivity, energy, geometry & frequencies.
	- path_list. List of lists, containing the fetched paths from single_path_finder() or poly_path_finder()
	Output:
	- Gsub. Subgraph limited to the nodes in the path input.'''
	# Get all unique nodes belonging to the paths and generate the subgraph
	selected_nodes = [node for path in path_list for node in path]
	unique_nodes = list(set(selected_nodes))
	Gsub = G.subgraph(unique_nodes).copy()
	fmap = formula_dict_constructor(Gsub)
	# Set the list of paths as a graph property, including the transition states
	#Gsub.graph["path_list"] = path_reformatter(G,path_list)
	Gsub.graph["path_list"] = path_list
	return Gsub

# Profile generation: generate complete profiles instead of isolated INT/TS/INT triplets
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
	- neigh_connections. List of lists, containing NODE1,NODE2 strings for the valid neighbors of nodein
	'''
	neighbors = [nd for nd in G.neighbors(nodein)]
	edges = [(nodein,nd) for nd in neighbors]
	edgenames = [G.edges[edg]['name'] for edg in edges]
	# but we only want those neighbors coming from NEW edges
	clean_ndx = [ie for ie,edg in enumerate(edgenames) if edg not in edgelist]
	clean_neighbors = [neighbors[ie] for ie in clean_ndx]
	neigh_connections = [[nodein,nodeout] for nodeout in clean_neighbors]
	return neigh_connections

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
			current_connection = neighborizer(G,nd,edge_bookkeep)
			current_neigh = [item[1] for item in current_connection]
			next_neighbors.append(current_neigh)
			# Keep track of the profiles
			if (count > 0):
				prev = track_dict[count - 1][ii][jj]
				if (current_connection):
					current_prof = [prev + item[1] for item in current_connection]
				else:
					current_prof = prev
					final_profiles.append(current_prof)
				track_dict[count].append(current_prof)
			else:
				track_dict[count].append(current_connection)

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
		cyclic_profiles.append(cycle)
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
		# By now, only add these branches as NODE1,NODE2 pairs
		branches = [[item[0],item[1]] for item in new_edges]
		return cycles + branches

def theor_profile_plotter(G,profile_list,figsize=(10,10),cmap="tab10"):
	'''For a list of profiles obtained through theor_profile_builder(), theor_cycle_branch_builder(), or
	poly_path_finder(), fetch energy values from the graph G and plot the profile through Matplotlib
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

	for ii,nodes in enumerate(profile_list):
		# Recover energies for nodes (as strings) and edges (as tuples of node names)
		edges = [(nodes[jj],nodes[jj+1]) for jj in range(len(nodes) - 1)]
		Enodes = [simple_prop_fetch(G,nd,"energy") for nd in nodes]
		Eedges = [simple_prop_fetch(G,ed,"energy") for ed in edges]
		# transform back to profile
		prof_energies_raw = [(End,Eed) for End,Eed in zip(Enodes[:-1],Eedges)]
		prof_energies = [item for sublist in prof_energies_raw for item in sublist] + [Enodes[-1]]
		# and insert the NODE1 - EDGE - NODE2... structure 
		edgenames = [G.edges[ed]["name"] for ed in edges]
		prof_raw = [[nodes[jj],edgenames[jj]] for jj in range(len(edgenames))]
		prof = [item for sublist in prof_raw for item in sublist] + [nodes[-1]]

		# Barrierless TSs (TSb) shall be removed from prof and prof_energies
		kept_ndx = [ii for ii,entry in enumerate(prof) if "TSb" not in entry]
		prof_clean = [prof[ii] for ii in kept_ndx]
		prof_energies_clean = [prof_energies[ii] for ii in kept_ndx]
		prof = prof_clean
		prof_energies = prof_energies_clean

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

