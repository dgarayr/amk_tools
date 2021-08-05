import bokeh.plotting
import bokeh.palettes
import bokeh.models as bkm
import bokeh.transform
import bokeh.layouts
import argparse
import sys
import RXReader as arx
from jsmol_bokeh_extension import JSMol
import networkx as nx
import numpy as np

# Basic management functions

def name_prop_fetch(G,key,prop):
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

def generate_inline_mol(xyzstring):
	'''Create a JSMol visualization from a newline-separated XYZ-format block
	Input:
	- xyzstring. String, XYZ formatted block of text, separated by newlines
	Output:
	- mol_script_string. String, containing a load data statement to load the molecule in JSMol'''
	xyzsubst = xyzstring.replace("\n","|")
	mol_script_string = r"data 'model X'|" + xyzsubst + r"|end 'model X';show data" 
	return mol_script_string 

def generate_inline_mol_freqs(data_dict):
	'''Create a JSMol visualization including vibrations: processes a whole ColumnDataSource object from
	bokeh to handle all required data: frequency values, mol. coordinates and displacements
	Input:
	- data_dict. Dictionary containing geometry, name, frequencies and vibr_displace fields (to make it compatible
	with both ColumnDataSource objects or with node/edge views)
	- bkdata. ColumnDataSource object containing geometry, name, frequencies and vibr_displace fields.
	- ndx. Integer, index of the requested structure
	Output:
	- model_freq_block. String, containing a full LOAD DATA statement for JSMol, with separate models for each vibration'''
	# The try-except block allows to skip structures that lack required fields
	try:
		geo = data_dict["geometry"]
		name = data_dict["name"]
		frq = data_dict["frequencies"]
		displ = data_dict["vibr_displace"]
	except:
		return None
	# Handle cases where the edge has no info
	if (not (geo and frq and displ)):
		return None
	# we need to merge each displacement with the coordinates (and with a charge placeholder), and add the freq value as a comment
	freqvals = frq.split(";")
	coords = geo.split("\n")
	Nat = len(coords)
	displ_strings = [[" ".join(line) for line in dx] for dx in displ]
	model_freq_block = "load data 'model FRQS'\n"
	for jj,dx in enumerate(displ_strings):
		out_list = [line_crd + " 0 " + line_displ for line_crd,line_displ in zip(coords,dx)]	
		model_freq_block += "%d\n" % Nat
		model_freq_block += "frq = %.2f cm-1 \n" % float(freqvals[jj])
		model_freq_block += ("\n".join(out_list) + "\n")

	model_freq_block += "end 'model FRQS'"
	return model_freq_block	

def xyz_from_atoms(atomblock,comment=""):
	'''Generate a XYZ block from a block of text containing atoms and positions separated
	by newlines, as stored in the "geometry" property in RXReader
	Input:
	- atomblock. String, newline-separated geometry specification.
	- comment. String, text to be added in the comment line of the XYZ file
	Output:
	- block. String, newline-separated block in XYZ format (including no. of atoms and an optional comment)'''
	Nat = len(atomblock.split("\n"))
	block = "%d\n%s\n" % (Nat,comment)
	block += atomblock
	return block

def add_models(G):
	'''Add JSMol-compatible models from generate_inline_mol() for all nodes and edges in a graph
	Input:
	- G. nx.Graph object as generated from RXReader.
	Output:
	- Adds the "model" property to every node and edge
	'''
	for nd in G.nodes(data=True):
		geo = nd[1]["geometry"]
		xyz = xyz_from_atoms(geo,comment=nd[0])
		nd[1]["model"] = generate_inline_mol(xyz)

	for ed in G.edges(data=True):
		geo = ed[2]["geometry"]
		# Do not fail when there are empty blocks (e.g. for barrierless transition states)
		if (geo):
			xyz = xyz_from_atoms(geo,comment=ed[2]["name"])
			ed[2]["model"] = generate_inline_mol(xyz)
		else:
			ed[2]["model"] = None
	return None

def add_vibr_models(G):
	'''Add JSMol-compatible models from generate_inline_mol() for the vibrational displacements
	for all nodes and edges in a graph -> which must have been initialized through vibr_displ_parser()
	the property "vibr_displace" is available
	Input:
	- G. nx.Graph object as generated from RXReader -> with vibr_displ_parser().
	Output:
	- Adds the "vibr_models" property to every node and edge
	'''
	keys = ["geometry","name","frequencies","vibr_displace"]
	for nd in G.nodes(data=True):
		model_freq_block = generate_inline_mol_freqs(nd[1])
		nd[1]["vibr_models"] = model_freq_block

	for ed in G.edges(data=True):
		# Do not fail when there are empty blocks (e.g. for barrierless transition states)
		model_freq_block = generate_inline_mol_freqs(ed[2])
		ed[2]["vibr_models"] = model_freq_block

	return None

def generate_applet(local=False,local_route=None):
	'''JSMol applet generation for a molecule with a given XYZ block, borrowed
	from example in original repo of jsmol_bokeh_extension. Relies on three main elements
	Input
	- local. Boolean, if True, use a local instance of JSMol, located in local_route
	- local_route. String, path to the local JSMol instance
	Output:
	- info_dict. Dictionary with basic information to run JSMol.
	- applet. JSMol applet.
	- script_source. bokeh.ColumnDataSource object, allowing interactivity
	'''
	if (local and local_route):
		serverURL = "file://%s/php/jsmol.php" % local_route
		j2spath = "file://%s/j2s" % local_route
	else:
		serverURL = "https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php"
		j2spath = "https://chemapps.stolaf.edu/jmol/jsmol/j2s"
	
	info_dict = dict(
		height="100%",
		width="100%",
		serverURL=serverURL,
		use="HTML5",
		j2sPath=j2spath,
		script="background white;" 
	)
	# Instantiate a new source
	script_source = bkm.ColumnDataSource({"script":[]})
	# Build applet
	applet = JSMol(
		width=600,
		height=600,
		script_source=script_source,
		info=info_dict,
	)
	return info_dict,applet,script_source

def bokeh_network_view(G,positions=None,width=800,height=600,graph_title="Reaction network representation"):
	'''
	Generation of an interactive visualization of a RXReader-processed reaction network, via the
	from_networkx() method in Bokeh. 
	Input:
	- G. nx.Graph object as generated from RXReader.
	- positions. Dictionary, defining a layout for the graph nodes, as generated by any NetworkX layout generator.
	- width, height. Floats, dimensions of the figure
	- graph_title. String, label for the graph.
	Output:
	- bfig. bokeh.plotting.Figure() object, containing the graph, the corresponding node/edge labels & the basic interactive tools (pan, zoom,
	  reset, hovering...)
	- Gbok. Bokeh plot generated via from_networkx(),
	'''

	# Instantiate figure and graph, setting out basic parameters
	bfig = bokeh.plotting.Figure(title=graph_title,width=width,height=height,tools="pan,wheel_zoom,box_zoom,reset",
								x_range=bkm.Range1d(-1.2,1.2),y_range=bkm.Range1d(-1.2,1.2))
	bfig.axis.visible = False
	bfig.xgrid.grid_line_color = None
	bfig.ygrid.grid_line_color = None

	# Check if positions were passed, and create inside if necessary
	if (not positions):
		print("Generating network layout")
		positions = nx.spring_layout(G)

	Gbok = bokeh.plotting.from_networkx(G,layout_function=positions,scale=1,center=(0,0))
	bfig.renderers.append(Gbok)

	# Modify how nodes are rendered: access the ColumnDataSource for energies and build the transformation
	nodesource = Gbok.node_renderer.data_source
	evalues = nodesource.data["energy"]
	colormap = bokeh.transform.linear_cmap(field_name="energy",palette="Viridis256",
									  low=min(evalues),high=max(evalues))

	Gbok.node_renderer.glyph = bkm.Circle(size=30,fill_color=colormap,fill_alpha=0.5)
	Gbok.node_renderer.selection_glyph = bkm.Circle(size=30,fill_color=colormap,fill_alpha=0.8,line_width=2)
	
	# For edges, only modify selection glyph by now
	edgesource = Gbok.edge_renderer.data_source
	Gbok.edge_renderer.selection_glyph = bkm.MultiLine(line_width=3,line_color="orange")

	# Add labels according to the positioning object, first for nodes and then for edges
	nodenames,coords = zip(*positions.items())
	xn,yn = zip(*coords)
	posnode_dict = {'xn':xn,'yn':yn,'nnames':nodenames}
	posnode = bkm.ColumnDataSource(posnode_dict)
	labels_node = bkm.LabelSet(x='xn', y='yn', text='nnames', text_font_style='bold',
				  x_offset=0, y_offset=5, source=posnode, render_mode='canvas')
	bfig.add_layout(labels_node)

	# TS labels: fetch central position and name for every edge
	posedge_dict = {'xe':[],'ye':[],'enames':[],'etuples':[]}
	for ed in G.edges(data=True):
		xy1,xy2 = [positions[nn] for nn in ed[0:2]]
		mid = 0.5*(xy1+xy2)
		output = [mid[0],mid[1],ed[2]["name"],(ed[0],ed[1])]
		for ik,key in enumerate(['xe','ye','enames','etuples']):
			posedge_dict[key].append(output[ik])
	posedge = bkm.ColumnDataSource(posedge_dict)
	labels_edge = bkm.LabelSet(x='xe', y='ye', text='enames', text_font_style='bold',text_color='red',
				  x_offset=0, y_offset=0, source=posedge, render_mode='canvas')
	bfig.add_layout(labels_edge)

	# Adding tools: hover & tap. To be able to see info about both nodes and edges, we must change the inspection policy
	# and add the two renderers (nodes & edges) explicitly to the constructor. 
	Gbok.inspection_policy = bkm.EdgesAndLinkedNodes()
	# Select the fields to be shown (with @) and optionally add formatting between curly braces

	# Custom JS selector for node hovering to have the formula field in the products
	hoverJS = '''
	var nrend = graph.node_renderer.data_source
	if (cb_data.index.indices.length > 0) {
		var ndx = cb_data.index.indices[0]
		var formula_list = nrend.data["formula"]
		var formula = formula_list[ndx]
		if (!formula){
		hover.tooltips = [["tag","@name"],["E","@energy{%.2f}"]]
 		} else {
		hover.tooltips = [["tag","@name"],["E","@energy{%.2f}"],["formula","@formula"]]
		}
	}
	'''

	hover_node = bkm.HoverTool(description="Node hover",renderers=[Gbok.node_renderer],formatters={"@energy":"printf"})
	hover_node.callback = bkm.CustomJS(args=dict(hover=hover_node,graph=Gbok),code=hoverJS)
	hover_edge = bkm.HoverTool(description="Edge hover",tooltips=[("tag","@name"),("E","@energy{%.2f}")],
							   formatters={"@energy":"printf"},renderers=[Gbok.edge_renderer])

	bfig.add_tools(hover_node)
	bfig.add_tools(hover_edge)

	Gbok.selection_policy = bkm.EdgesAndLinkedNodes()
	tap = bkm.TapTool(renderers=[Gbok.node_renderer,Gbok.edge_renderer])
	bfig.add_tools(tap)

	return bfig,Gbok

def profile_datasourcer(G,profile_list):
	'''Convert a list of profiles (obtained through arx.theor_profile_builder(), arx.theor_cycle_branch_builder(), or
	arx.poly_path_finder()) to a list of ColumnDataSource objects that can be used by Bokeh, with the structure
	ii | jj | energy | label, where ii is the index of the profile and jj the order of the species
	Input:
	- G. nx.Graph object as generated from RXReader.
	- profile_list. List of lists of strings containing labels defining energy profiles, including TSs (via arx.path_reformatter()).
	Output:
	- cds_list. List of ColumnDataSource objects defining all required energy profiles'''
	cds_list = []
	for ii,prof in enumerate(profile_list):
		# Filter out barrierless TS 
		kept_ndx = [ii for ii,entry in enumerate(prof) if "TSb" not in entry]
		cleaned_prof = [prof[ii] for ii in kept_ndx]
		# fix TS names, taking the corresponding nodes as reference
		working_prof = [(prof[jj-1],prof[jj+1]) if "TS" in entry else entry for jj,entry in enumerate(cleaned_prof)]
		# Prepare the lines, as in arx.theor_profile_plotter(), duplicating entries
		xvals = np.arange(0,2*len(working_prof))
		energies = [arx.simple_prop_fetch(G,entry,"energy") for entry in working_prof]
		yvals = [e for evalue in energies for e in (evalue,evalue)]
		labels = [lab for labvalue in cleaned_prof for lab in (labvalue,labvalue)]
		# save as ColumnDataSource
		prof_cds = bkm.ColumnDataSource(data={"x":xvals,"y":yvals,"lab":labels})
		cds_list.append(prof_cds)
	return cds_list

def profile_bokeh_plot(G,profile_list,condition=[]):
	'''Generate a Bokeh figure for energy profiles, using profile_datasourcer() to transform the profile list to a 
	ColumnDataSource list.
	Input:
	- G. nx.Graph object as generated from RXReader.
	- profile_list. List of lists of strings containing labels defining energy profiles, including TSs (via arx.path_reformatter()).
	Output:
	- bfig. Bokeh figure containing line plots & labels for every profile in profile_list
	'''
	# Initialization: instantiate figure, remove X-axis and add palette
	bfig = bokeh.plotting.Figure(width=600,height=600,tools="pan,wheel_zoom,reset",name="PROFILE")
	bfig.xaxis.visible = False
	bfig.yaxis.axis_label = "E/kcal mol-1"
	palette = bokeh.palettes.d3['Category10'][10]

	# Generate the list of ColumnDataSources
	cds_paths = profile_datasourcer(G,profile_list)

	# Check condition
	if (len(condition) != len(profile_list)):
		condition = [True for prof in profile_list]

	# Iterate & plot
	for ii,cdspath in enumerate(cds_paths):
		Nentries = len(cdspath.data["lab"])
		cndx = (ii % 10)
		rx_line = bfig.line(x="x",y="y",color=palette[cndx],source=cdspath)
		# Prepare labels, slicing the original CDS to avoid repetitions
		cds_lab = bkm.ColumnDataSource({k:cdspath.data[k][::2] for k in ["x","y","lab"]})
		cds_lab.data["x"] = [(float(item) + 0.5) for item in cds_lab.data["x"]]
		rect = bkm.Rect(x="x",y="y",fill_color=palette[cndx],fill_alpha=0.9,width=1,height=2,line_width=0,name="RECTANGLE")
		rx_rect = bfig.add_glyph(cds_lab,rect)
		rx_label = bkm.LabelSet(x="x",y="y",text="lab",source=cds_lab,x_offset=-1,y_offset=1,name="THELABELS")
		bfig.add_layout(rx_label)
		### Add filtering: hide when filter condition is not fulfilled
		rx_line.visible = condition[ii]
		rx_rect.visible = condition[ii]
		rx_label.visible = condition[ii]

	# And add the hover
	hover_prof = bkm.HoverTool(description="Profile hover",tooltips=[("tag","@lab"),("E","@y{%.2f}")],
							   formatters={"@y":"printf"},renderers=bfig.renderers)
	bfig.add_tools(hover_prof)
	bfig.js_on_event('reset',bkm.CustomJS(args = {"prof":bfig}, code = js_callback_dict["resetProfile"]))
	return bfig
		
js_callback_dict = {
	# Dictionary storing JS code for bokeh.models.CustomJS definitions
	"loadMolecule":"""
		// from graph, we fetch nrend - node renderer, erend - edgerenderer
		//source - source object for JSMol
		var nrend = graph.node_renderer.data_source
		var erend = graph.edge_renderer.data_source
		if (cb_obj.indices.length) {
			var ninds = nrend.selected.indices
			var einds = erend.selected.indices
			if (ninds.length) {
				var rend = nrend
			} else {
				var rend = erend
			}
			var ndx = rend.selected.indices[0]
			var model = rend.data["model"][ndx]
			var molname = rend.data["name"][ndx]
			textField.text = "<font size=+1><b>" + molname + "</b></font>"
			source.data["script"] = [model]
			source.properties.data.change.emit()
		}
		""",

	"loadVibrations":"""
		// from graph, we fetch nrend - node renderer, erend - edgerenderer, 
		// source - source object for JSMol, menu - dropdown menu to change vibrations
		var nrend = graph.node_renderer.data_source
		var erend = graph.edge_renderer.data_source
		
		var ninds = nrend.selected.indices
		var einds = erend.selected.indices
		if (ninds.length) {
			var rend = nrend
		} else {
			var rend = erend
		}
		var ndx = rend.selected.indices[0]
		var vibr_models = rend.data["vibr_models"][ndx]
		if (vibr_models.length){
			source.data["current_model"] = [vibr_models]
			source.data["script"] = [vibr_models + "; vibration on ; vibration scale 0.25"]
			source.data["current_mol"] = rend.data["name"][ndx]
			// and also populate the normal mode dropdown: get the frequencies and assign them to menu
			// we must pass indices as strings for it to work, add 1 to match model numbering
			menu.disabled = false
			var freqlist = rend.data["frequencies"][ndx].split(";")
			var menulist = freqlist.map(function (frq,ii) {return [frq+" cm-1",(ii+1).toString()]})
			menu.menu = menulist
			var nw_label = "Normal modes (" + source.data["current_mol"] + ")"
			menu.label = nw_label
			// modify text element using the BASE frequency
			var molname = source.data["current_mol"]
			var frqval = freqlist[0].toString().trim() + " cm-1"
			textField.text = "<font size=+1><b>" + molname + " (" + frqval + ")" + "</b></font>"
		}
	
		source.properties.data.change.emit()
		""",

	"molToClipboard":"""
		// from graph, we fetch nrend - node renderer, erend - edgerenderer
		// source - source object for JSMol, button - button object
		var nrend = graph.node_renderer.data_source
		var erend = graph.edge_renderer.data_source

		var ninds = nrend.selected.indices
		var einds = erend.selected.indices
		if (ninds.length) {
			var rend = nrend
		} else {
			var rend = erend
		}
		var ndx = rend.selected.indices[0]
		var geo = rend.data["geometry"][ndx]
		var name = rend.data["name"][ndx]
		var orig_label = button.label
		// use a function to apply setTimeout and recover original label
		function recover_label(button,orig_label){
			button.label = orig_label
		}	

		if (geo.length){
			navigator.clipboard.writeText([geo])
			button.label = "Copied " + name + "!"
			setTimeout(recover_label,2000,button,orig_label)
		}
		""",

		"locateMolecule":"""
		// source - source object for JSMol
		// pass graph and fetch node and edge renderers
		// from fig, we modify x_range and y_range. Default plot starts from -1.2 to 1.2,
		var nrend = graph.node_renderer.data_source
		var erend = graph.edge_renderer.data_source
		var layout = graph.layout_provider.graph_layout
		// fetch the query in the data sources, choosing the appropiate renderer depending on the query
		var mol_query = text_input.value
		if (mol_query.includes("TS")) {
			var renderer = erend
		} else {
			var renderer = nrend
		}
		var pool_names = renderer.data["name"]
		var ndx = pool_names.indexOf(mol_query)
		// locate positions of the node or of the nodes defining an edge
		if (mol_query.includes("TS")) {
			var n1 = renderer.data["start"][ndx]
			var n2 = renderer.data["end"][ndx]
			var pos1 = layout[n1]
			var pos2 = layout[n2]
			var positions = new Array(2)
			positions[0] = 0.5*(pos1[0]+pos2[0])
			positions[1] = 0.5*(pos1[1]+pos2[1])
		} else {
			var positions = layout[mol_query]
		}
		// returns -1 if element is not present
		if (ndx >= 0) {
			renderer.selected.indices = [ndx]
			fig.x_range.start = positions[0] - 0.5
			fig.x_range.end = positions[0] + 0.5
			fig.y_range.start = positions[1] - 0.5
			fig.y_range.end = positions[1] + 0.5
		}
		""",

		"chooseVibrationMenu":"""
		// source - source object for JSMol ; menu - dropdown menu object
		source.data["script"] = [source.data["current_model"] + "; vibration on ; vibration scale 0.25 ; model " + this.item]
		var ndx = parseInt(this.item) - 1
		var molname = source.data["current_mol"]
		var frqval = menu.menu[ndx][0]
		source.data["current_vibration"] = frqval
		textField.text = "<font size=+1><b>" + molname + " (" + frqval.trim() + ")" + "</b></font>"
		source.properties.data.change.emit()
		""",

		"replacePlot":"""
		//layout - full layout, jsmol - jsmol app (where the profile window will appear), prof - profile

		var current_fig = layout.children[1][0].children[1]
		var current_plot = current_fig.children[0]
		var controls = current_fig.children[1]
		if (current_plot.properties.name.spec.value == "PROFILE"){
			var side_array = [jsmol,controls]
		} else {
			var side_array = [prof,controls]
		}
		// we need to modify the children list, does not work for single substitutions
		layout.children[1][0].children[1].children = side_array
		""",

		"selectProfileByMol":"""
		// graph - bokeh graph for the network view, use to fetch node and edge renderers
		// prof - bokeh figure for the profiles
		var nrend = graph.node_renderer.data_source
		var erend = graph.edge_renderer.data_source
		var ninds = nrend.selected.indices
		var einds = erend.selected.indices
		if (ninds.length) {
			var rend = nrend
		} else if (einds.length) {
			var rend = erend
		} else {
			var rend = null
			// make everything visible
			for (let [index,line] of prof.renderers.entries()){
				prof.renderers[index].visible = true
				prof.center[index+2].visible = true
			} 
			return true
		}
		var ndx = rend.selected.indices[0]
		var sel_spc = rend.data["name"][ndx]
		// labels are stored in "center" property of the figure
		// and we have TWO renderers per entry: skeleton and rectangle
		var Nlines = prof.renderers.length/2
		for (let i = 0 ; i < Nlines ; i++){
			const index = 2*i
			var current_data = prof.renderers[index].data_source
			var labels = current_data.data["lab"]
			// use ! to negate inclusion
			if (!labels.includes(sel_spc)){
				prof.renderers[index].visible = false
				prof.renderers[index+1].visible = false
				// use i to map the labels, which are not duplicated
				prof.center[i+2].visible = false
			} 
		prof.properties.renderers.change.emit()
		}
		""",

		"selectProfileByEnergy":"""
		// graph - bokeh graph for the network view, use to fetch node and edge renderers
		// prof - bokeh figure for the profiles, thrbox - textbox for threshold energy
		var thr = parseFloat(thrbox.value)
		var range = {min: thr, max: thr}
		function overThreshold(element){
			return element > thr
		}
		var Nlines = prof.renderers.length/2
		for (let i = 0 ; i < Nlines ; i++){
			const index = 2*i
			var current_data = prof.renderers[index].data_source
			var energies = current_data.data["y"]
			//hide if anything is above the threshold
			if (energies.some(overThreshold)){
				prof.renderers[index].visible = false
				prof.renderers[index+1].visible = false
				// use i to map the labels, which are not duplicated
				prof.center[i+2].visible = false
			} 
		}
		""",

		"resetProfile":"""
		// prof - bokeh figure for the profiles
		// make everything visible
		var Nlines = prof.renderers.length/2
		for (let i = 0 ; i < Nlines ; i++){
			const index = 2*i
			prof.renderers[index].visible = true
			prof.renderers[index+1].visible = true
			prof.center[i+2].visible = true
		}
		""" 
}

def full_view_layout(bokeh_figure,bokeh_graph,G=None,local_jsmol=False,local_jsmol_route=None):
	'''
	Combination of the interactive graph visualization generated by bokeh_network_view with a JSMol instance and the
	corresponding interaction buttons. Molecules are shown by clicking in each edge or node.
	All callbacks use JavaScript code, not Python.
	Nested functions are defined here to allow data passing.
	Input:
	- bokeh_figure, bokeh_graph. Output of bokeh_network_view()
	- G. NetworkX.Graph(). Used for profile support.
	- local_jsmol. Boolean, if True, use a local instance of JSMol, located in local_route
	- local_jsmol_route. String, path to the local JSMol instance
	'''
	# Add access to data sources for nodes and edges
	nodesource = bokeh_graph.node_renderer.data_source
	edgesource = bokeh_graph.edge_renderer.data_source
	
	# Instantiate the applet
	info,app,script_source = generate_applet(local_jsmol,local_jsmol_route)

	# Instantiate the required widgets: buttons, text inputs, menus...
	b1 = bkm.Button(label="Load vibrations",width=100) 
	b2 = bkm.Button(label="Geometry to clipboard",width=100)
	text_input = bkm.TextInput(value=nodesource.data["name"][0],width_policy="fit")
	b3 = bkm.Button(label="Locate molecule",width_policy="fit")
	menu = bkm.Dropdown(label="Normal modes",menu=[("No vibr. loaded","None"),None],disabled=True,width=200)
	spc1 = bkm.Spacer(width=650)
	text = bkm.Div(text="",width=300)
	# For profile-based elements
	b4 = bkm.Button(label="Molec. filter",width=100)
	cbox = bkm.CheckboxGroup(labels=["Show profile"])
	b5 = bkm.Button(label="Energy filter",width=100)
	thrbox = bkm.TextInput(value="%.2f" % max(edgesource.data["energy"]),width=100)

	# Write the JavaScript callback to allow to avoid the Bokeh server: all is ported to JavaScript
	# For this to work, we need to pre-load all models in the graph
	js_load_mol = bkm.CustomJS(args = {"graph":bokeh_graph,"source":script_source,"textField":text}, 
							   code = js_callback_dict["loadMolecule"])

	js_load_vibrations = bkm.CustomJS(args = {"graph":bokeh_graph,"source":script_source,"menu":menu,"textField":text}, 
									  code = js_callback_dict["loadVibrations"])

	js_geo_clipboard = bkm.CustomJS(args = {"graph":bokeh_graph,"source":script_source,'button':b2}, 
									code = js_callback_dict["molToClipboard"])

	js_mol_locator = bkm.CustomJS(args = {"graph":bokeh_graph,"fig":bokeh_figure,"text_input":text_input},
								  code = js_callback_dict["locateMolecule"])
	
	js_menu_selector = bkm.CustomJS(args = {"source":script_source,"menu":menu,"textField":text},
					 			    code = js_callback_dict["chooseVibrationMenu"])

	# Set up the JS callbacks for nodes and edges
	bokeh_graph.node_renderer.data_source.selected.js_on_change('indices',js_load_mol)
	bokeh_graph.edge_renderer.data_source.selected.js_on_change('indices',js_load_mol)

	# Widget callbacks
	# Button 1: load vibrational models
	b1.js_on_click(js_load_vibrations)
	# Button 2: pass current geometry to clipboard
	b2.js_on_click(js_geo_clipboard)
	# Allow to select by name
	b3.js_on_click(js_mol_locator)
	# Vibration selector
	menu.js_on_event("menu_item_click",js_menu_selector)

	# Layout: 
	# Upper row: controls & info. 
	# Lower row: main visualizations
	# | Load vibr.	|Select mode | To clipboard |  | Mol. info 		|
	# | 	Network visualization	|		JSMol instance			|
	# | Location text | Loc. button |	Profile options & filters	|
	row1 = bkm.Row(b1,menu,b2,spc1,text)
	row2 = bkm.Row(bkm.Column(bokeh_figure,bkm.Row(text_input,b3)),bkm.Column(app))
	layout = bokeh.layouts.grid([row1,row2])

	# Add the options to see profiles if G was passed to the function and contains a pathList property
	if (G and "pathList" in G.graph):
		# Define elements (checkbox & button) and append them to the column containing the JSMol widget
		# Then instantiate figure and prepare callbacks
		layout.children[1][0].children[1].children.append(bkm.Row(cbox,b4,thrbox,b5))
		fig_prof = profile_bokeh_plot(G,G.graph["pathList"])
		js_plot_modif = bkm.CustomJS(args = {"layout":layout,"prof":fig_prof,"jsmol":app}, 
									 code = js_callback_dict["replacePlot"])
		cbox.js_on_click(js_plot_modif)
		js_select_profile_mol = bkm.CustomJS(args = {"graph":bokeh_graph,"prof":fig_prof}, 
											 code = js_callback_dict["selectProfileByMol"])
		b4.js_on_click(js_select_profile_mol)
		js_select_profile_e = bkm.CustomJS(args = {"prof":fig_prof,"thrbox":thrbox}, 
										   code = js_callback_dict["selectProfileByEnergy"])
		b5.js_on_click(js_select_profile_e)
	return layout

def generate_visualization(G,title,outfile,workfolder,Nvibrations=-1,with_profiles=False):
	'''Wrapper function to generate HTML visualizations for a given network
	Input:
	- G. nx.Graph object as generated from RXReader. For profile support, it should contain a graph["pathList"] property.
	- title. String, title for the visualization.
	- outfile. String, name of the output HTML file.
	- workfolder. String, name of the folder to fetch calculations from (FINAL_LL/FINAL_HL)
	- Nvibrations. Integer, no. of vibrations to be loaded in the visualization. 0: skip vibration loading, -1, load all
	- with_profiles. Boolean, if True, add profile visualization.
	Output:
	- lay. Bokeh layout as generated by full_view_layout()
	'''
	posx = nx.spring_layout(G)
	# Add model field to all nodes and edges & also vibrations
	add_models(G)
	if (Nvibrations != 0):
		arx.vibr_displ_parser(workfolder,G,Nvibrations)
		add_vibr_models(G)
	
	# Bokeh-powered visualization via RXVisualizer
	bk_fig,bk_graph = bokeh_network_view(G,positions=posx,graph_title=title)
	if (with_profiles):
		# Call full_view_layout() passing the graph so it loads profiles
		lay = full_view_layout(bk_fig,bk_graph,G)
	else:
		lay = full_view_layout(bk_fig,bk_graph)

	bokeh.plotting.output_file(outfile,title=title)
	bokeh.plotting.save(lay)
	
	return lay

def pass_args_cmd():
	'''Use argparse to pass command-line arguments to generate a visualization
	Input:
	- Args from STDIN
	Output:
	- args. argparse.ArgumentParser() object.
	'''
	argparser = argparse.ArgumentParser()
	argparser.add_argument("workdir",help="Directory with AutoMeKin FINAL calculations",type=str)
	argparser.add_argument("rxnfile",help="Name of the RXNet file to use. Options: RXNet, RXNet.cg, RXNet.rel",type=str)
	argparser.add_argument("--barrierless",help="Include barrierless routes from RXNet.barrless",action='store_true')
	argparser.add_argument("--vibrations",help="Number of normal modes to add to visualization: use -1 for all",type=int,default=-1,metavar="Nvibr")
	argparser.add_argument("--paths",help="Include paths from SOURCE to TARGET",type=str,nargs="*",metavar=("SOURCE","TARGET"))
	argparser.add_argument("--cutoff_path",help="Set cutoff for the paths",type=int,default=4,metavar="CUTOFF")
	argparser.add_argument("--unreduced",help="Use full graph, without selecting nodes in paths",action='store_true')
	argparser.add_argument("--outfile",help="Name for the HTML output file",type=str,default="network.html")
	argparser.add_argument("--title",help="Title for the HTML visualization",type=str,default="Reaction network visualization")
	try:
		args = argparser.parse_args()
	except:
		print("workdir and rxnfile must be passed")
		sys.exit()
	return args

def gen_view_cmd(args):
	'''Generate a visualization via generate_visualization() through the commandline arguments handled by pass_args_cmd()
	Input:
	- args. Filled argparse.ArgumentParser() object.
	Output:
	- Gwork. NetworkX.Graph object generated by RX_builder(), possibly including paths in pathList property.
	- view. Bokeh visualization from full_view_layout().
	'''
	# Read the graph
	data = arx.RX_parser(workfolder=args.workdir,rxnfile=args.rxnfile,check_connection=True)
	if (args.barrierless):
		data_barrless = arx.RX_parser(workfolder=args.workdir,rxnfile="RXNet.barrless")
		joined_data = [data[ii]+data_barrless[ii] for ii in range(len(data))]
		data = joined_data
	G = arx.RX_builder(workfolder=args.workdir,data=data)
	# Path handling: several situations are possible: i) SOURCE and TARGET, ii) only SOURCE, iii) not SOURCE nor TARGET, iv) --paths not passed
	# For i), ii) and iii), paths shall be added, and args.paths will be a LIST
	paths_passed = isinstance(args.paths,list)
	if (paths_passed):
		# Fully specified path
		if (len(args.paths) == 2):
			source,target = args.paths
			limits = arx.node_synonym_search(G,[source,target])
		elif (len(args.paths) == 1):
			source = args.paths
			target = None
			limits = arx.node_synonym_search(G,[source])
			limits.append([])
		else:
			source = None
			target = None
			limits = [[],[]]
		#args.title = " to ".join([source,target])
		paths = arx.add_paths(G,limits[0],limits[1],cutoff=args.cutoff_path)
		if (not args.unreduced):
			Gwork = arx.graph_path_selector(G,paths)
		else:
			Gwork = G
	else:
		Gwork = G
	# Generate the visualization
	view = generate_visualization(Gwork,title=args.title,outfile=args.outfile,workfolder=args.workdir,
								  Nvibrations=args.vibrations,with_profiles=paths_passed)

	return Gwork,view