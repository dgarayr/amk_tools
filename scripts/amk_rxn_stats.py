#!/usr/bin/env python3
'''Read a reaction network and compute its statistics through arx.stat_generator'''
import RXReader as arx
import argparse
import sys
from pathlib import Path
argparser = argparse.ArgumentParser()
argparser.add_argument("finaldir",help="Directory with AutoMeKin FINAL calculations",type=str)
argparser.add_argument("--rxnfile",'-r',help="Name of the RXNet file to use. Default is CG. Options: RXNet, RXNet.cg, RXNet.rel",
					   type=str,default="RXNet.cg")
argparser.add_argument("--skip_barrierless",'-sb',help="Exclude barrierless routes from RXNet.barrless",
					   action='store_true')
argparser.add_argument("--Ngraphs",'-Ng',help="Number of random graphs to consider for statistics",
					   type=int,default=1000)

try:
	args = argparser.parse_args()
except:
	print("finaldir must be passed")
	argparser.print_help()
	sys.exit()

# Network reading
fol = args.finaldir
net = args.rxnfile
barr_path = Path(fol + "/RXNet.barrless")
data = arx.RX_parser(fol,net,check_connection=True)
if (not args.skip_barrierless and barr_path.is_file()):
	data_barrless = arx.RX_parser(fol,rxnfile="RXNet.barrless")
	joined_data = [data[ii]+data_barrless[ii] for ii in range(len(data))]
	data = joined_data
Gx = arx.RX_builder(fol,data)

arx.stat_generator(Gx,Ng=args.Ngraphs,gen_file=True)
