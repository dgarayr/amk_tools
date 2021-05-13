import bokeh.plotting
import bokeh.models as bkm
import bokeh.transform
import bokeh.layouts
from jsmol_bokeh_extension import JSMol
import networkx as nx
import numpy as np

# Basic management functions

def generate_inline_mol(xyzstring):
	'''Create a JMol visualization for a valid XYZ block'''
	xyzsubst = xyzstring.replace("\n","|")
	mol_script_string = r"data 'model X'|" + xyzsubst + r"|end 'model X';show data" 
	return mol_script_string 

def xyz_from_atoms(atomblock,comment=""):
	'''Generate a XYZ block from a block of text containing atoms and positions separated
	by newlines, as stored in the "geometry" property in RXReader'''
	Nat = len(atomblock.split("\n"))
	block = "%d\n%s\n" % (Nat,comment)
	block += atomblock
	return block

def generate_applet():
	'''JSMol applet generation for a molecule with a given XYZ block, borrowed
	from example in original repo of jsmol_bokeh_extension'''

	# We need three main elements: a dictionary with basic information to run JSMol
	# a Bokeh ColumnDataSource object for reactivity (which we will alter later on), and the applet	
	# Here we also start with a default molecule
	#mol_script = generate_inline_mol(molstring)
	info_dict = dict(
		height="100%",
		width="100%",
		serverURL="https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php",
		use="HTML5",
		j2sPath="https://chemapps.stolaf.edu/jmol/jsmol/j2s",
		script="background white;" 
	)
	# Instantiate a new source
	script_source = bkm.ColumnDataSource()
	# Build applet
	applet = JSMol(
		width=600,
		height=600,
		script_source=script_source,
		info=info_dict,
	)
	return info_dict,applet,script_source


def bokeh_network_view(G,positions=None):
	'''
	Use the from_networkx() method in Bokeh to build a visualization of the reaction network as
	processed by RXReader, and add interactivity (hovering & picking)
	'''

	# Start by instantiating the Bokeh figure and the graph that we will append there
	bfig = bokeh.plotting.figure(title="Graph",width=800,height=600,tools="pan,wheel_zoom,box_zoom,reset")
	# Check if positions were passed, and create inside instead
	if (not positions):
		print("Generating network layout")
		positions = nx.spring_layout(G)

	Gbok = bokeh.plotting.from_networkx(G,layout_function=positions)
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

def full_view_layout(bokeh_figure,bokeh_graph):
	'''From the hover-able graph visualization, add the JSMol widget and the interaction buttons.
	Nested functions here allow interactivity'''

	def change_mol():
		nodesource = bokeh_graph.node_renderer.data_source
		edgesource = bokeh_graph.edge_renderer.data_source
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

	def click_callback(attr,old,new):
		# Just to allow to change the molecule by direct clicking
		change_mol()


	# Instantiate the applet
	info,app,script_source = generate_applet()

	# Add the second tap
	tap_change = bkm.TapTool(renderers=[bokeh_graph.node_renderer,bokeh_graph.edge_renderer],gesture="doubletap")

	bokeh_figure.add_tools(tap_change)
	# Set up the callback for nodes and edges (and then the change_mol will take care of which shall be displayed)
	for source in [bokeh_graph.node_renderer,bokeh_graph.edge_renderer]:
		source.data_source.selected.on_change('indices',click_callback) 

	layout = bokeh.layouts.gridplot([bokeh_figure,app],ncols=2)

	return layout