import bokeh.plotting
import bokeh.models as bkm
import bokeh.transform
import bokeh.layouts
from jsmol_bokeh_extension import JSMol
import networkx as nx
import numpy as np

# Basic management functions

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
		xyz = xyz_from_atoms(geo,comment=ed[2]["name"])
		ed[2]["model"] = generate_inline_mol(xyz)

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
		model_freq_block = generate_inline_mol_freqs(ed[2])
		ed[2]["vibr_models"] = model_freq_block

	return None

def generate_applet():
	'''JSMol applet generation for a molecule with a given XYZ block, borrowed
	from example in original repo of jsmol_bokeh_extension. Relies on three main elements
	- Dictionary with basic information to run JSMol
	- bokeh ColumnDataSource object, allowing interactivity
	- JSMol applet	
	'''
	info_dict = dict(
		height="100%",
		width="100%",
		serverURL="https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php",
		use="HTML5",
		j2sPath="https://chemapps.stolaf.edu/jmol/jsmol/j2s",
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

def bokeh_network_view(G,positions=None,width=800,height=600):
	'''
	Generation of an interactive visualization of a RXReader-processed reaction network, via the
	from_networkx() method in Bokeh. 
	Input:
	- G. nx.Graph object as generated from RXReader.
	- positions. Dictionary, defining a layout for the graph nodes, as generated by any NetworkX layout generator.
	- width, height. Floats, dimensions of the figure
	Output:
	- bfig. bokeh.plotting.Figure() object, containing the graph, the corresponding node/edge labels & the basic interactive tools (pan, zoom,
	  reset, hovering...)
	- Gbok. Bokeh plot generated via from_networkx(),
	'''

	# Instantiate figure and graph, setting out basic parameters
	bfig = bokeh.plotting.Figure(title="Reaction network representation",width=width,height=height,tools="pan,wheel_zoom,box_zoom,reset",
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
	tooltips = [("tag","@name"),("E","@energy{%.2f}")]
	hovering = bkm.HoverTool(tooltips=tooltips,formatters={"@energy":"printf"},renderers=[Gbok.node_renderer,Gbok.edge_renderer])
	bfig.add_tools(hovering)

	Gbok.selection_policy = bkm.EdgesAndLinkedNodes()
	tap = bkm.TapTool(renderers=[Gbok.node_renderer,Gbok.edge_renderer])
	bfig.add_tools(tap)

	return bfig,Gbok

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
			var freqlist = rend.data["frequencies"][0].split(";")
			var menulist = freqlist.map(function (frq,ii) {return [frq+" cm-1",(ii+1).toString()]})
			menu.menu = menulist
			var nw_label = "Normal modes (" + source.data["current_mol"] + ")"
			menu.label = nw_label
			// modify text element using the BASE frequency
			var molname = source.data["current_mol"]
			var frqval = freqlist[1].toString() + " cm-1"
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
		// from graph, we fetch nrend (node renderer), erend and layout
		// from fig, we modify x_range and y_range. Default plot starts from -1.2 to 1.2,
		var nrend = graph.node_renderer.data_source
		var layout = graph.layout_provider.graph_layout
		// fetch the query in the data sources
		var mol_query = text_input.value
		var all_names = nrend.data["name"]
		var ndx = all_names.indexOf(mol_query)
		// returns -1 if element is not present
		if (ndx >= 0) {
			nrend.selected.indices = [ndx]
			var positions = layout[mol_query]
			fig.x_range.start = positions[0] - 1.0
			fig.x_range.end = positions[0] + 1.0
			fig.y_range.start = positions[1] - 1.0
			fig.y_range.end = positions[1] + 1.0
		}
		""",

		"chooseVibrationMenu":"""
		// source - source object for JSMol ; menu - dropdown menu object
		source.data["script"] = [source.data["current_model"] + "; vibration on ; vibration scale 0.25 ; model " + this.item]
		var ndx = parseInt(this.item) - 1
		var molname = source.data["current_mol"]
		var frqval = menu.menu[ndx][0]
		source.data["current_vibration"] = frqval
		textField.text = "<font size=+1><b>" + molname + " (" + frqval + ")" + "</b></font>"
		source.properties.data.change.emit()
		"""
}

def full_view_layout(bokeh_figure,bokeh_graph,py_callbacks=True):
	'''
	Combination of the interactive graph visualization generated by bokeh_network_view with a JSMol instance and the
	corresponding interaction buttons. Molecules are shown by clicking in each edge or node.
	Nested functions are defined here to allow data passing.
	Input:
	- bokeh_figure, bokeh_graph. Output of bokeh_network_view()
	- py_callbacks. Boolean, if True, use Python callbacks. Else, use JavaScript callbacks.
	'''
	# Add access to data sources for nodes and edges
	nodesource = bokeh_graph.node_renderer.data_source
	edgesource = bokeh_graph.edge_renderer.data_source
	def load_mol():
		'''Generation of a JSMol model for a given species in the network and loading in the widget'''
		# check node and edge sources
		is_node,is_edge = [bool(source.selected.indices) for source in [nodesource,edgesource]]
		# select the source
		if (is_node and not is_edge):
			sel_source = nodesource
		elif (is_edge and not is_node):
			sel_source = edgesource
		else:
			return None
		ndx = sel_source.selected.indices[0]
		# Fetch geometry and name, then convert to XYZ and pass to JSMol
		geo = sel_source.data["geometry"][ndx]
		name = sel_source.data["name"][ndx]
		xyz = xyz_from_atoms(geo,comment=name)
		mol = generate_inline_mol(xyz)
		script_source.data['script'] = [mol]

	def load_vibration():
		'''Generation of a JSMol model collection (including vibrations) for a given species in the network and loading in the widget'''
		# check node and edge sources
		is_node,is_edge = [bool(source.selected.indices) for source in [nodesource,edgesource]]
		# select the source
		if (is_node and not is_edge):
			sel_source = nodesource
		elif (is_edge and not is_node):
			sel_source = edgesource
		else:
			return None
		ndx = sel_source.selected.indices[0]
		# Add a try-except block for cases where frequencies have not been read or for structures that did not include freqs (such as PROD)
		# And also generate a dictionary with the required parameters
		keys = ["geometry","name","frequencies","vibr_displace"]
		try:
			mol_dict = {k:sel_source.data[k][ndx] for k in keys}
			model_freq_block = generate_inline_mol_freqs(mol_dict)
		except:
			return None
		script_source.data['script'] = [model_freq_block + "; vibration on ; vibration scale 0.5"]

	def click_callback(attr,old,new):
		'''Wrapper to use load_mol() as a callback upon selection changes: molecule loading upon node/edge clicking'''
		load_mol()

	# Instantiate the applet
	info,app,script_source = generate_applet()

	# Instantiate the required widgets: buttons, text inputs, menus...
	b1 = bkm.Button(label="Load vibrations",width=100) 
	b2 = bkm.Button(label="Geometry to clipboard",width=100)
	text_input = bkm.TextInput(value=nodesource.data["name"][0],width_policy="fit")
	b3 = bkm.Button(label="Locate molecule",width_policy="fit")
	menu = bkm.Dropdown(label="Normal modes",menu=[("No vibr. loaded","None"),None],disabled=True,width=200)
	spc1 = bkm.Spacer(width=650)
	text = bkm.Div(text="",width=300)

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

	# Set up the callbacks for nodes and edges, either from Python or from JS
	if (py_callbacks):
		for source in [bokeh_graph.node_renderer,bokeh_graph.edge_renderer]:
			source.data_source.selected.on_change('indices',click_callback) 
	else:
		bokeh_graph.node_renderer.data_source.selected.js_on_change('indices',js_load_mol)
		bokeh_graph.edge_renderer.data_source.selected.js_on_change('indices',js_load_mol)

	# Widget callbacks
	# Button 1: load vibrational models
	if (py_callbacks):
		b1.on_click(load_vibration)
	else:
		b1.js_on_click(js_load_vibrations)
	
	# Button 2: pass current geometry to clipboard
	b2.js_on_click(js_geo_clipboard)
	
	# Allow to select by name
	b3.js_on_click(js_mol_locator)

	# Vibration selector
	menu.js_on_event("menu_item_click",js_menu_selector)
	
	# Layout: left column with network and locator, right row with buttons for JSMol control
	# Row-based layout???
	row1 = bkm.Row(b1,menu,b2,spc1,text)
	row2 = bkm.Row(bkm.Column(bokeh_figure,bkm.Row(text_input,b3)),app)
	layout = bokeh.layouts.grid([row1,row2])
	return layout