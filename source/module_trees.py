import os
# import re
import sys
import datetime
import subprocess as sp

# from time import sleep
# from hashlib import md5
from string import strip, lower
from socket import getfqdn
# from getpass import getuser
# from operator import itemgetter
from module_utils import lookForDirectory, lookForFile, splitSequence, \
  format_time

''' Module which implements the functionality for reconstructing phylogenetic
    trees. It contains wrappers to three different programs: PhyML, RAxML &
    FastTree
'''

## Exit code meanings
##  80: Not enough sequences
exit_codes = {
  "phyml":          96,    ##  96: Problems associated to PhyML execution
  "codonphyml":     97,    ##  97: Problems associated to CodonPhyML execution
  "raxml":          98,    ##  98: Problems associated to RAxML execution
  "fasttree":       99,    ##  99: Problems associated to FastTree execution

  "generic":        95,    ##  95: Program not supported
}

def phylogenetic_trees(parameters):
  ''' Phylogenetic trees are reconstructed according to the input parameters.
      Once the different files have been generated, the function moves those
      files into a pre-established filename schema
  '''

  ## Get output folder/generic filename
  oFile = os.path.join(parameters["out_directory"], parameters["prefix"])

  current_directory = os.getcwd()
  ## Change current directory to the output folder. Any temporary file will be
  ## generated therefore in this folder
  os.chdir(parameters["out_directory"])

  ## Set output filename and log file
  if parameters["replace"] and parameters["step"] == 0:
    logFile = open(oFile + ".log", "w")
  else:
    logFile = open(oFile + ".log", "a+")

  start = datetime.datetime.now()
  date = start.strftime("%H:%M:%S %m/%d/%y")
  print >> logFile, ("###\n###\tSTEP\tPhylogenetic Tree Reconstruction\tSTART\t"
    + "%s\n###") % (date)
  logFile.flush()

  ## Get which program will be used to reconstruct phylogenetic trees. Check
  ## such program is listed among the available binaries
  if not "tree" in parameters:
    sys.exit("ERROR: Check your configuration file. There is no definition for "
      + "the Phylogenetic TREE reconstruction step")

  prog = parameters["tree"][0]
  if not prog in parameters:
    sys.exit(("ERROR: Selected program '%s' is not available accordding to the "
      "the configuration file") % (prog))

  ## Get binary as well as any default parameters for the selected program
  binary = parameters[prog]
  key = ("%s_params") % (prog)
  progr_params = parameters[key] if key in parameters else ""

  if not "evol_models" in parameters:
    sys.exit("ERROR: Check your configuration file. There is no definition for "
      + "the <evol_models> parameter")

  ## If the evolutionary model list is not appropiately formated, do it
  if isinstance(parameters["evol_models"], basestring):
    parameters["evol_models"] = map(strip, parameters["evol_models"].split())

  ## Check if <numb_models parameters is defined and how many models are
  ## requested to be evaluated
  if not "numb_models" in parameters or parameters["numb_models"].lower() \
    == "all":
    parameters["numb_models"] = len(parameters["evol_models"])
  parameters["numb_models"] = int(parameters["numb_models"])

  if not parameters["numb_models"] in range(1,len(parameters["evol_models"])+1):
    sys.exit(("ERROR: Check how many evolutionary models has been asked to re"
      + "construct '%d'") % (parameters["numb_models"]))

  ## Check which approaches should be used for the phylogenetic reconstruction
  ## and whether there are specific program's parameters for them
  if not "tree_approach" in parameters:
    parameters["tree_approach"] = ["ml"]

  ## Remove potential duplicates and lowercase all approaches for the tree
  ## reconstruction
  parameters["tree_approach"] = set(map(lower, parameters["tree_approach"]))

  ## We will first loot for Neighbour Joining tree reconstruction, then for
  ## Maximum likelihood and then for any other approach defined in the config
  ## file
  tree_approaches = []
  if "nj" in parameters["tree_approach"]:
    tree_approaches.append("nj")
  if "ml" in parameters["tree_approach"]:
    tree_approaches.append("ml")
  others = parameters["tree_approach"] - set(["nj", "ml"])
  if others != set():
    tree_approaches += sorted(others)

  selected_models = parameters["evol_models"]
  ## Reconstruct trees for each approach considering evolutionary models order
  ## according their likelihood values
  for approach in tree_approaches:

    ## Save results - we will use such data for selecting the best -if required-
    ## models fitting to the input data
    results = {}

    ## Format the choosen program's parameters according to the default ones and
    ## the specific ones for the current approach
    exec_params = ("%s ") % (progr_params)
    exec_params += parameters[approach] if approach in parameters else ""

    for model in selected_models:
      out_file = ("%s.%s.tree.%s.%s.nw") % (oFile, prog, approach, model)
      stats_file = ("%s.%s.tree.%s.%s.st") % (oFile, prog, approach, model)

      if prog in ["phyml", "codonphyml"]:
        exec_params = ("%s -m %s") % (exec_params, model)

      elif prog in ["fasttree"]:
        ## On FastTree is selected by default JTT model for AAs - so we don't
        ## set-up that model
        if model.lower() != "jtt":
          exec_params = ("%s -%s") % (exec_params, model)

      elif prog in ["raxml"]:
        ## -n run name - it couldnot contain "/"
        ## -p randon number
        ## -s input alignment name
        exec_params = ("%s%s") % (exec_params, model)



      if perform_tree(prog, binary, exec_params, parameters["in_file"],
        out_file, stats_file, logFile, parameters["replace"]):
          parameters["replace"] = True



  #~ results = []
  #~ ## Explore the different models and get the likelihood of each of them
  #~ for model in evolutionary_models:
#~
    #~ ## Call to the appropiate wrapper depending on the selected
    #~ if program == "phyml":
      #~ cmd = ("%s -i %s %s -m %s") % (parameters["phyml"], parameters["inFile"],\
        #~ parameters[approach], model)
#~
      #~ lk = wrapperPhyML(cmd, parameters["inFile"], statsFile, outFile,
        #~ parameters["log"], parameters["replace"])
#~
      #~ results.append((lk, model))
#~
  #~ results = [(pair[1], pair[0]) for pair in sorted(results, reverse = True)]
#~
  #~ outFile = ("%s.%s.tree.rank.%s") % (common, program, approach)
  #~ if parameters["replace"] or not utils.lookForFile(outFile):
    #~ ranking = "\n".join(["\t".join(map(str, pair)) for pair in results])
    #~ print >> open(outFile, "w"), ranking
#~
  #~ return [pair[0] for pair in results]

def perform_tree(label, binary, parameters, in_file, out_file, stats_file, \
  logFile, replace):

  '''
  Function to format the command-line of different phylogenetic tree reconstruc-
  tion programs and execute such command lines.
  '''

  ## Check whether the output file already exists. If it is not set to replace
  ## it, just return to the calling function
  if lookForFile(out_file) and not replace:
    return False

  if label in ["phyml", "codonphyml"]:
    cmd = ("%s -i %s %s") % (binary, in_file, parameters)

  elif label in ["fasttree"]:
    cmd = ("%s %s -log %s -out %s %s") % (binary, parameters, out_file, \
      stats_file, in_file)

  else:
    sys.exit(exit_codes["generic"])

  ## Record the time and precise command-line
  name = getfqdn()
  start = datetime.datetime.now()
  date = start.strftime("%H:%M:%S %m/%d/%y")

  print >> logFile, ("###\n###\t%s - Phylogenetic Trees\t") % (label.upper()),
  print >> logFile, ("%s\n###\t[%s]\tCommand-line\t%s\n###") % (date, name, cmd)
  logFile.flush()

  try:
    proc = sp.Popen(cmd, shell = True, stderr = logFile, stdout = logFile)
  except OSError, e:
    print >> sys.stderr, "ERROR: Execution failed: " + str(e)
    sys.exit(exit_codes[label])

  ## Press a key until any situation to continue the execution
  proc.communicate("Y\n")
  if proc.wait() != 0:
    print >> sys.stderr, ("ERROR: Execution failed: %s") % (label.upper())
    sys.exit(exit_codes[label])

  final = datetime.datetime.now()
  total = format_time((final - start).seconds if start else 0)
  print >> logFile, ("###\tTime\t%s\n###") % (total)
  logFile.flush()

  ## Process program's output and rename output files according to our own
  ## scheme
  if label in ["phyml", "codonphyml"]:
    try:
      sp.call(("mv %s_%s_tree.txt %s") % (in_file, label, out_file), shell = \
        True)
      sp.call(("mv %s_%s_stats.txt %s") % (in_file, label, stats_file), shell =
        True)
    except OSError:
      print >> sys.stderr, ("ERROR: Impossible to rename %s output files") \
        % (label.upper())
      sys.exit(exit_codes[label])

  elif label in ["raxml"]:
    try:
      sp.call(("mv RAxML_bestTree%s %s") % (suffix, outFile), shell = True)
      sp.call(("mv RAxML_info%s %s") % (suffix, statsFile), shell = True)
    except OSError:
      print >> sys.stderr, ("ERROR: Impossible to rename RAxML output files")
      sys.exit(exit_codes[label])

    oFile = open(statsFile, "a+")
    for oth_file in utils.listDirectory(outdirec, suffix):
      fileName = os.path.split(oth_file)[1]
      hz_line = "#" * (len(fileName) + 4)
      print >> oFile, ("%s\n%s\n%s") % (hz_line, fileName, hz_line)
      print >> oFile, ("%s") % ("".join(open(oth_file, "rU").readlines()))
      sp.call(("rm -f %s") % (oth_file), shell = True)
    oFile.close()

  return True

def get_likelihood(label):


  ## PHYML
  ## Independtly whether the files are exist or not, get the likelihood associate

  logLK = None
  ## Get likelihood value for the execution
  try:
    pipe = sp.Popen(("grep 'Log-likelihood:' %s") % (statsFile), shell = True, \
      stdout = sp.PIPE)
    logLK = float(map(strip, pipe.stdout.readlines()[0].split())[2])
  except OSError, e:
    sys.exit("ERROR: Impossible to get LK values for PhyML run")

  ## RAXML
  ## Parse output to get the likelihood of the best tree
  for line in open(statsFile, "rU"):
    if line.lower().startswith("final") and line.lower().find("score") != -1:
      logLK = float(map(strip, line.split())[-1])

  ## FastTree
    ## Get likelihood values for the current run
  for line in open(statsFile, "rU"):
    f = map(strip, line.split("\t"))
    try:
      value = float(f[2])
    except:
      continue
    logLK = value if not logLK or value < logLK else logLK

  ## Return the likelihood value for the current tree
  return logLK
