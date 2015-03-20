#! /usr/bin/python -u
from string import *
import sys
from commands import *
sys.path.insert(0, "/users/rg/mmariotti/libraries/")
sys.path.append('/users/rg/mmariotti/scripts')
from MMlib import *

help_msg=""" Script that wraps the scripts of the Phylomizer pipeline by Salvador Capella-Gutierrez. See https://github.com/Gabaldonlab/phylomizer
Usage:   phylomizer_wrapper.py  [specify routines and input]  -o output_folder  [other options]

Routines:  each require the options listed just below, unless coming from a preceding step. 
Note that activating a routine does not automatically activate the next ones.
____________________________________________________________________________________________
-H    build Homology group for a single sequence in input, using blast against a database.
   -q  +      protein query single-fasta file for homology search
   -db +      protein database for homology search. NOTE: it must include the query sequence
____________________________________________________________________________________________
-A    build Alignments
   -f  +      unaligned multi-fasta file containing an homology group
____________________________________________________________________________________________
-P    run Phyml including model optimization
   -a  +      aligned multi-fasta file containing an homology group
____________________________________________________________________________________________

### other options:
-o           +      output folder. By default it is derived from input name
-r                  force overwrite of existing files (reccomended for any second run)
-sp                 filter the input sequences to keep only the selenoprofiles titles (makes sense only if -A or -P)
-ALL                run all routines; like -H -A -P. Same requirements as for -H
-bin         +      folder where the phylomizer python scripts are searched
-config      +      config file for phylomizer
-python      +      path to python executable used to run the phylomizer scripts
-print_opt          print currently active options
-h OR --help        print this help and exit"""

command_line_synonyms={'c':'config'}

def_opt= {'temp':'/users/rg/mmariotti/temp', 
'H':0, 'A':0, 'P':0,
'config':'/users/tg/scapella/000.007.NewPipeline/config.pipeline.dc', 
'bin':'/users/tg/scapella/000.007.NewPipeline/',
'python':'python',
'q':0, 'db':0, 'f':0, 'a':0,
'ALL':0,
'o':0, 'r':0,
'v':0, 'sp':0
}

def is_selenoprofiles_title(title):  return ' chromosome:' in title

gawk_parser_line= """gawk 'BEGIN {i=1} NR>1 && /./ {if (! titles_over){titles[i] = ">" $1 "\\n" $2; i+=1}else{titles[i] = titles[i] $1; i+=1}}NR>1 && !/./{titles_over=1; i=1}END{ for (x =1; x<=length(titles); x++){print titles[x]}}' """ 
#########################################################
###### start main program function

def main(args={}):
#########################################################
############ loading options
  global opt
  if not args: opt=command_line(def_opt, help_msg, 'io', synonyms=command_line_synonyms )
  else:  opt=args
  set_MMlib_var('opt', opt)
  #global temp_folder; temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)
  #global split_folder;    split_folder=Folder(opt['temp']);               test_writeable_folder(split_folder); set_MMlib_var('split_folder', split_folder) 
  #checking input

  bin_folder=    Folder(opt['bin'])
  python_executable=    opt['python']
  config_file=          opt['config']
  output_folder= Folder(opt['o'])

  ### check input options; H, A, P and dependencies
  if opt['ALL']: opt['H'], opt['A'], opt['P']=1, 1, 1
  q_control, db_control, f_control, a_control=0, 0, 0, 0
  query_file, db_file, unaligned_file, aligned_file = None, None, None, None

  write('*** Active routines: '+ 'Homology'*int(bool(opt['H']))  +'  '+'Alignment'*int(bool(opt['A']))+'  '+'Phylogeny'*int(bool(opt['P'])) +' ****', 1)
  if opt['H']:    
    q_control, db_control, f_control, a_control=1, 1, 0, 0
    executable_name='11.IndividualStep.HomologySearch.py'
    main_input_name=opt['q']
    if opt['A']:       
      executable_name='21.CombinedSteps.HomologySearch_Alignments.py'      
      if opt['P']:    executable_name='10.pipeline.py'
    elif opt['P']:    raise Exception, "ERROR you cannot specify option -H (homology search) and -P (phylogeny) without running -A as well (alignment)."
  elif opt['A']:
    q_control, db_control, f_control, a_control=0, 0, 1, 0
    executable_name='12.IndividualStep.Alignments.py'
    if opt['P']:      executable_name='22.CombinedSteps.Alignments_PhylogeneticTrees.py'
  elif opt['P']:      
    q_control, db_control, f_control, a_control=0, 0, 0, 1
    executable_name='13.IndividualStep.PhylogeneticTrees.py'
  else:               raise Exception, "ERROR no routines requested! You must specify at least one option among: -H -A -P   (see program help with -h)"
    
  ## determining automatic output folder name from input
  if not opt['o']: 
    if q_control:     output_base= opt['q']
    elif f_control:   output_base= opt['f']
    elif a_control:   output_base= opt['a']
    output_folder=Folder(  join(output_base.split('.')[:-1], '.')+'.phylomizer')  
  else: output_folder=Folder(opt['o'])
  ## prefix name (fixed)
  prefix_output='sequences'
  options=''
  if opt['r']: options+=' -r '

  # sp option, for inputing results by selenoprofiles
  if opt['sp']:
    if opt['H']: raise Exception, "ERROR option -sp not compatible with -H"
    elif opt['A']: 
      filtered_ali_file=output_folder+prefix_output+'.selenoprofiles.fasta'
      input_file=opt['f']
      opt['f']=filtered_ali_file
    elif opt['P']:
      filtered_ali_file=output_folder+prefix_output+'.ali.selenoprofiles.fasta'
      input_file=opt['a']
      opt['a']=filtered_ali_file      
    filtered_ali_file_h=open(filtered_ali_file, 'w')
    sp_titles=0; non_sp_titles=0
    for title, seq in parse_fasta(input_file):
      if is_selenoprofiles_title(title):        
        sp_titles+=1
        print >> filtered_ali_file_h, ">" + title +'\n'+seq
      else: 
        non_sp_titles+=1
    filtered_ali_file_h.close()
    write('Option -sp:   accepted '+str(sp_titles)+' selenoprofiles titles from input ('+str(non_sp_titles)+' were discarded)', 1)

  id2title={}
  changed_titles=False
  ## controls
  if q_control and (not opt['q'] or not is_file(opt['q'])) :   raise Exception, "ERROR with these routines, it is required a single protein sequence fasta file with option -q"
  elif not q_control and opt['q']:                             printerr("WARNING option -q was specified but it is not utilized with these routines!", 1)
  elif q_control:                                              options+=' -i '+opt['q'] 
  if db_control and (not opt['db'] or not is_file(opt['db'])): raise Exception, "ERROR with these routines, it is required a protein database file with option -db"
  elif not db_control and opt['db']:                           printerr("WARNING option -db was specified but it is not utilized with these routines!", 1)    
  elif db_control:                                             options+=' -d '+opt['db'] 
  if f_control and (not opt['f'] or not is_file(opt['f'])):    raise Exception, "ERROR with these routines, it is required an unaligned fasta file with homologous sequences with option -f"
  elif not f_control and opt['f']:                             printerr("WARNING option -f was specified but it is not utilized with these routines!", 1)    
  elif f_control:                                              
    seq_collection=alignment(opt['f'])  #actually not aligned, but won't crash here
    if  any(   [ len(t.split()[0])>=40   for t in seq_collection.titles()]  ):  #at least one title too long for phyml/t_coffee; we should correct them now!
      changed_titles=True
      id2title_file= opt['f']+'.id2title.tab'
      new_fasta_file= opt['f']+'.w_ids'; new_fasta_file_h=open(new_fasta_file, 'w')
      printerr('WARNING one or more of the titles in the sequence file (-f) are too long! Producing version with short ids: '+new_fasta_file, 1)
      for t in seq_collection.titles():
        index_to_make_uniq=0 #used first time when it is 1
        new_id=t.split()[0][:40]
        while new_id in id2title:
          index_to_make_uniq+=1
          new_id = t.split()[0][  : (40-1-len(str(index_to_make_uniq)) )   ]+'_'+str(index_to_make_uniq)    ## final: longid_1  or longid_2   or longi_12  (but long 40 chars)
        id2title[new_id]=t
        print >> new_fasta_file_h, ">"+new_id+'\n'+fasta(  seq_collection.seq_of(t) )
      new_fasta_file_h.close()
      opt['f']=new_fasta_file
    options+=' -i '+opt['f'] 

  if a_control and (not opt['a'] or not is_file(opt['a'])):    raise Exception, "ERROR with these routines, it is required an aligned fasta file with homologous sequences with option -a"
  elif not a_control and opt['a']:                             printerr("WARNING option -a was specified but it is not utilized with these routines!", 1)    
  elif a_control:                                                  
    phylip_format_aligned_file=output_folder+prefix_output+'.phylip'
    a=alignment(  opt['a']  )   
    if any(   [ len(t.split()[0])>=40   for t in a.titles()]  ): #at least one title too long for phyml; we should correct them now!
      changed_titles=True
      id2title_file= phylip_format_aligned_file+'.id2title.tab'
      printerr('WARNING one or more of the titles in the aligned sequences file (-a) are too long! Producing version with short ids', 1)
      for t in a.titles():
        index_to_make_uniq=0 #used first time when it is 1
        new_id=t[:40]
        while new_id in id2title:
          index_to_make_uniq+=1
          new_id = t[  : (40-1-len(str(index_to_make_uniq)) )   ]+'_'+str(index_to_make_uniq)    ## final: longid_1  or longid_2   or longi_12  (but long 40 chars)
        id2title[new_id]=t
      for new_id in id2title:
        a.change_title(id2title[new_id], new_id)
    if not a.check_length(): raise Exception, "ERROR file provided with -a is not aligned fasta!"
    write_to_file( a.phylip_format(),  phylip_format_aligned_file )    
    options+=' -i '+ phylip_format_aligned_file

  if changed_titles: 
    write('Writing correspondances of titles with shortened ids --> '+id2title_file, 1)
    id2title_file_h=open(id2title_file, 'w')
    for new_id in id2title:
      print >> id2title_file_h, new_id+'\t'+id2title[new_id]
    id2title_file_h.close()

  ##### RUNNING 
  write('Output folder: '+output_folder, 1)  
  cmnd= python_executable+' '+ bin_folder + executable_name + ' -c '+config_file+' -o '+output_folder+' -p '+prefix_output  + options ### r replace
  write( "Running ... ", 1); write(cmnd, 1 )
  bbash(cmnd)
  
  if opt['H']: 
    #bbash(gawk_parser_line+' '+output_folder+prefix_output+'.seqs  > '+output_folder+prefix_output+'.seqs.fasta' ) 
    write('--> Homology group fasta sequence file: '+output_folder+prefix_output+'.seqs', 1)
  if opt['A']:
    final_ali_file= output_folder+prefix_output+'.alg.clean.fasta'
    bbash(gawk_parser_line+' '+output_folder+prefix_output+'.alg.clean  > '+final_ali_file)
    if changed_titles: 
      write('Restoring titles in output alignment ... ', 1)
      a=alignment(final_ali_file)
      for new_id in list(a.titles()): #making sure it's copied, thus fixed
        title=id2title[new_id]
        a.change_title(new_id, title)
      a.display(final_ali_file)
    write('--> Alignment file (converted to fasta): '+final_ali_file, 1)
  if opt['P']:
    possible_ml_tree_files= bbash('ls '+output_folder+prefix_output+'.tree.phyml.ml.*.nw', 1).split('\n')
    if len(possible_ml_tree_files)>1:
      rank_ml_file=output_folder+prefix_output+'.tree.phyml.rank.ml'
      rank_ml_file_h=open(rank_ml_file); first_line=rank_ml_file_h.readline(); rank_ml_file_h.close()
      best_ml_tree_file= output_folder+prefix_output+'.tree.phyml.ml.'+  first_line.split()[0]  +'.nw'      
    else:
      best_ml_tree_file= possible_ml_tree_files[0]

    if changed_titles: 
      from ete2 import Tree
      tree=Tree(best_ml_tree_file)
      for new_id in id2title: 
        node=tree&new_id 
        node.name=    id2title[new_id].split()[0]
      write('Restoring titles in output tree ... ', 1)
      write_to_file(  tree.write(format=0), best_ml_tree_file+'.restored_titles' )  
      best_ml_tree_file= best_ml_tree_file+'.restored_titles'

    bbash('cd '+output_folder+'; ln -fs '+base_filename(best_ml_tree_file)+' '+prefix_output+'.best_tree')  #linking best tree
    write('--> Best tree topology file (linked): '+output_folder+prefix_output+'.best_tree', 1)

 ### converter for ali files: gawk 'BEGIN{index=1}NR>1 && /./{if ( ! titles_over ){titles[index] = ">" titles[index] "\n" $2; index+=1 } else{titles[index] = titles[index] $1; index+=1} }  NR>1 && !/./ {titles_over=1; index=1} END { for (x =1; x<=length(titles); x++){ print titles[x]} }  ' sequences.alg.clean    WORK IN PROGRESS

  

  ###############



#######################################################################################################################################

def close_program():
  if 'temp_folder' in globals() and is_directory(temp_folder):
    bbash('rm -r '+temp_folder)
  try:
    if get_MMlib_var('printed_rchar'): 
      printerr('\r'+printed_rchar*' ' ) #flushing service msg space       
  except:
    pass

  if 'log_file' in globals(): log_file.close()


if __name__ == "__main__":
  try:
    main()
    close_program()  
  except Exception:
    close_program()
    raise 
