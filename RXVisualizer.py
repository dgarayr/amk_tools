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
	freqvals = frq.split("\n")
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
								x_range=bkm.Range1d(-1.1,1.1),y_range=bkm.Range1d(-1.1,1.1))
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

	# Write the JavaScript callback to allow to avoid the Bokeh server: all is ported to JavaScript
	# For this to work, we need to pre-load all models in the graph
	js_load_mol = bkm.CustomJS(args = {'nrend':nodesource,'erend':edgesource,"source":script_source}, code = """
		if (cb_obj.indices.length) {
			var ninds = nrend.selected.indices
			var einds = erend.selected.indices
			if (ninds.length) {
				var rend = nrend
			} else {
				var rend = erend
			}
			var ndx = rend.selected.indices[0]
			var message = "selection callback"
			var model = rend.data["model"][ndx]
			source.data["script"] = [model]
			source.properties.data.change.emit()
		}
		""")

	js_load_vibrations = bkm.CustomJS(args = {'nrend':nodesource,'erend':edgesource,"source":script_source}, code = """
		var ninds = nrend.selected.indices
		var einds = erend.selected.indices
		if (ninds.length) {
			var rend = nrend
		} else {
			var rend = erend
		}
		var ndx = rend.selected.indices[0]
		var message = "selection callback"
		var vibr_models = rend.data["vibr_models"][ndx]
		console.log(vibr_models)
		if (vibr_models.length){
			source.data["script"] = [vibr_models + "; vibration on ; vibration scale 0.5"]
		}
		source.properties.data.change.emit()
		""")

	# Set up the callbacks for nodes and edges, either from Python or from JS
	if (py_callbacks):
		for source in [bokeh_graph.node_renderer,bokeh_graph.edge_renderer]:
			source.data_source.selected.on_change('indices',click_callback) 
	else:
		bokeh_graph.node_renderer.data_source.selected.js_on_change('indices',js_load_mol)
		bokeh_graph.edge_renderer.data_source.selected.js_on_change('indices',js_load_mol)
	# Prepare the button to load vibrations on a selected structure
	b1 = bkm.Button(label="Load vibration")
	if (py_callbacks):
		b1.on_click(load_vibration)
	else:
		b1.js_on_click(js_load_vibrations)

	# Dispose the layout: graph at left, column with JSMol and buttons right
	col2 = bkm.Column(app,b1)
	layout = bokeh.layouts.gridplot([bokeh_figure,col2],ncols=2)

	return layout