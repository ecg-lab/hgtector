#!/usr/bin/env python
"""Python implementation of HGTector method

Required python packages:
Biopython
Pandas
RPY2

By Tim Straub
Zhaxybayeva Lab
Dartmouth College
(C) 2014
timothy.j.straub.gr@dartmouth.edu

Original method developed by:
Qiyun Zhu, Michael Kosoy, and Katharina Dittmar
HGTector: an automated method facilitating genome-wide
discovery of putative horizontal gene transfers
BMC Genomics 2014, 15: 717
doi: 10.1186/1471-2164-15-717
"""

import sys, imp
try:
    imp.find_module("Bio")
    imp.find_module("rpy2")
    imp.find_module("pandas")
    imp.find_module("numpy")
    # imp.find_module("oursql")
except ImportError as e:
    print e
    print "Make sure you have the following modules installed:"
    print "rpy2, pandas, numpy, Biopython"
    sys.exit(1)
except Exception as e:
    print e
    sys.exit(1)

import os, glob, pandas, multiprocessing, time, numpy

import rpy2.robjects as ro
import pandas.rpy.common as com

# import oursql

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
from subprocess import Popen, list2cmdline, PIPE
from optparse import OptionParser

# taxonomic dump direcotry
TAXDUMP_DIRECTORY = "/home/tstraub/data/taxdmp_2014-11/"

# use this blast database
BLAST_DB = "/home/tstraub/data/db/nr_2014-11/nr"

# use this file to define multispecies RefSeq entries
MULTISPECIES_FILE = "/home/tstraub/data/release68.MultispeciesAutonomousProtein2taxname"

# ignore blast results with these in the name
IGNORE_TAXA = ["unknown",
               "uncultured",
               "unidentified",
               "unclassified",
               "environmental",
               "plasmid",
               "vector",
               "synthetic",
               "phage",
               "artificial",
              ]

# default config values
DEFAULT_CONFIG = {"title": "",
                  "interactive": 1,
                  "taxdump": TAXDUMP_DIRECTORY,
                  "nThreads": 1,
                  "dbBLAST": BLAST_DB,
                  "nHits": 500,
                  "evalue": 1e-5,
                  "ignoreTaxa": IGNORE_TAXA,
                  "selfGroup": None,
                  "closeGroup": None,
                  "modKCO": 1,
                  "selfLow": 0,
                  "modKCOSelf": None,
                  "modKCOClose": None,
                  "modKCODistal": None,
                  "topPercent": None,
                 }

def call_and_print(cmd):
    """Call a command with Popen and print its output
    """
    try:
        proc = Popen(cmd, stdout=PIPE, stderr=PIPE)
        (output, error) = proc.communicate()
        if output:
            return "%s\n%s" % (list2cmdline(cmd), output)
        elif error:
            return error
        else:
            return list2cmdline(cmd)
    except KeyboardInterrupt:
        if proc: proc.terminate()
    except Exception as e:
        print e
        if proc: proc.terminate()

def single_thread_cmds(cmds):
    """Call each command one by one
    """
    for cmd in cmds:
        try:
            print call_and_print(cmd)

        except KeyboardInterrupt:
            break
        except Exception as e:
            print e
            break

def __multi_call(args):
    """Wrapper for call_and_print with a queue
    """
    cmd, q = args

    output = call_and_print(cmd)
    q.put(output)
    return output

def thread_cmds(cmds, threads=1, quiet=False):
    """Multi-threaded call commands
    """
    if not (cmds and len(cmds)):
        print "Must specify commands"
        return None

    if not threads:
        print "Invalid thread count, set to 1"
        threads = 1

    elif type(threads) != int:
        try:
            threads = int(threads)
        except:
            print "Invalid thread count, set to 1"
            threads = 1

    elif threads < 1:
        print "Will use all available CPUs"
        threads = None

    if threads == 1:
        #print "Single thread specified. Will not multi-thread."
        single_thread_cmds(cmds)
    else:
        p = multiprocessing.Pool(threads)
        m = multiprocessing.Manager()
        q = m.Queue()
        args = [(cmd, q) for cmd in cmds]
        egm2 = p.map_async(__multi_call, args)
        i = 0
        results = None
        while True:
            try:
                if egm2.ready():
                    egm2.get()
                    p.close()
                    break
                else:
                    if i % 60 == 0 and not quiet:
                        print "There are %i remaining BLASTs to run" % (len(cmds) - q.qsize())
                    i += 1
                    time.sleep(1)
            except KeyboardInterrupt:
                p.terminate()
                break
            except Exception as e:
                print e
                p.terminate()
                break

def basename(name):
    """Get the basename of a file without its extension
    """

    if "." in name:
        return ".".join(os.path.basename(name).split(".")[:-1])
    else:
        return os.path.basename(name)

def get_gene_name(gene):
    """Get GI and RefSeq name for a gene text entry
    """

    if not gene:
        print "No gene"
        return ('', '')

    if '|' not in gene:
        if len(gene) > 10:
            temp = gene.split()
            return (temp[0], "")
        else:
            return (gene, "")

    temp = gene.split('|')
    if len(temp) >= 4:
        gi = temp[1]
        if '.' in temp[3]:
            ref = temp[3].split('.')[0]
        else:
            ref = temp[3]
    else:
        gi = temp[1]
        ref = ""

    return (gi, ref)

def read_config_file(file):
    """Read in the configuration file for HGTector
    """
    if not file:
        print "No config file"
        return

    config = DEFAULT_CONFIG

    with open(file, "rb") as f:
        for line in f:
            if line[0] == "#":
                continue
            if "=" not in line:
                continue
            temp = line.strip().split("=")
            if "," in temp[1]:
                temp[1] = temp[1].split(",")
            config[temp[0]] = temp[1]

    return config

def split_fasta(file, write_sequences=True):
    """Split a fasta file into one file per sequence
    """
    if not file:
        print "No file"
        return

    if write_sequences:
        if not os.path.isdir("blast"):
            os.mkdir("blast")

        workdir = os.path.join("blast", basename(file))

        if not os.path.isdir(workdir):
            os.mkdir(workdir)

    with open(file, "rU") as f:

        info = {}
        self_ids = {}
        seqs = SeqIO.parse(f, "fasta")
        for i, seq in enumerate(seqs):
            if seq.id:
                gene, ref_id = get_gene_name(seq.id)
                self_ids[gene] = ref_id
            elif seq.name:
                gene, ref_id = get_gene_name(seq.name)
                self_ids[gene] = ref_id
            else:
                gene = str(i)
            product = ' '.join(seq.description.split()[1:]).split('[')[0].strip()
            length = len(seq)
            info[gene] = (product, length)
            if write_sequences:
                name = "%s.fa" % gene
                with open(os.path.join(workdir, name), "wb") as w:
                    SeqIO.write(seq, w, "fasta")

        print "Found %i genes in %s" % (len(info), file)

        return info, self_ids

def get_blast_cmds(directory, blast_db=BLAST_DB, n_hits=500, evalue=1e-5, nucl=False):
    """Return blastp commands for a given directory of fasta files
    """
    if not directory:
        print "No directory"
        return

    if type(n_hits) != int:
        try:
            n_hits = int(n_hits)
        except:
            print "Using default nHits of 500"
            n_hits = 500

    if not (type(evalue) == float or type(evalue) == int):
        try:
            evalue = float(evalue)
        except:
            print "Using default evalue 1e-5"
            evalue = 1e-5

    commands = []
    in_progress = False
    for file in glob.glob(os.path.join(directory, "*.fa")):
        if not os.path.isfile(os.path.join(directory, "%s.bla" % basename(file))):
            # get command line
            if nucl:
                cline = blastn(query=file, db=blast_db, num_alignments=n_hits, evalue=evalue, dust='yes', out=os.path.join(directory, "%s.bla" % basename(file)), outfmt="6")
            else:
                cline = blastp(query=file, db=blast_db, num_alignments=n_hits, evalue=evalue, seg='yes', out=os.path.join(directory, "%s.bla" % basename(file)), outfmt="6")
            # convert to list for Popen
            commands.append(str(cline).split())
        else:
            # already a result file, so the BLAST was in progress
            # do not redo this BLAST command
            in_progress = True

    if in_progress:
        if not len(commands):
            print "BLAST results were complete"
            return []
        else:
            print "Found some BLAST results. Will not redo them."

    return commands

def load_taxids(taxdump=TAXDUMP_DIRECTORY, nucl=False):
    """Load information about GIs and TaxIDs to memory
    """
    if not taxdump:
        taxdump=TAXDUMP_DIRECTORY
    
    if nucl:
        file = os.path.join(taxdump, "gi_taxid_nucl.dmp")
    else:
        file = os.path.join(taxdump, "gi_taxid_prot.dmp")

    taxids = {}
    with open(file, "rb") as f:
        for line in f:
            (gi, taxid) = line.split()
            taxids[gi] = taxid
    
    return taxids

def load_nodes_file(taxdump=TAXDUMP_DIRECTORY):
    """Load parent/children information about TaxIDs to memory
    """
    if not taxdump:
        taxdump=TAXDUMP_DIRECTORY

    nodes = os.path.join(taxdump, "nodes.dmp")

    if not os.path.isfile(nodes):
        print "No nodes file"
        return

    with open(nodes, "rb") as f:
        parent_nodes = {}
        children_nodes = {}
        for line in f:
            temp = [x.strip() for x in line.split("|")]
            if len(temp) < 2:
                continue
            taxid = temp[0]
            parent = temp[1]

            if taxid == parent:
                continue

            if parent not in parent_nodes:
                parent_nodes[parent] = []

            parent_nodes[parent].append(taxid)

            if taxid in children_nodes:
                print "Warning, duplicate entry: %s" % taxid
            else:
                children_nodes[taxid] = parent

        return (parent_nodes, children_nodes)

def load_names_db(taxdump=TAXDUMP_DIRECTORY):
    """Load names.dmp into memory
    """
    if not taxdump:
        taxdump = TAXDUMP_DIRECTORY

    names_dmp = os.path.join(taxdump, "names.dmp")

    if not os.path.isfile(names_dmp):
        print "Cannot find names.dmp file"
        return

    names_db = {}
    with open(names_dmp, "rb") as f:
        for line in f:
            temp = [x.strip() for x in line.split("|")]
            if len(temp) < 2:
                continue

            taxid = temp[0]
            name = temp[1]

            if len(temp) >= 4:
                tax_class = temp[3]
            else:
                tax_class = None

            if taxid not in names_db:
                names_db[taxid] = []

            names_db[taxid].append((name, tax_class))

    return names_db

def read_multispecies(file=MULTISPECIES_FILE):
    """Read the multispecies RefSeq file
    """
    with open(file, "rb") as f:
        multi = {}
        for line in f:
            line = line.strip()
            (wp, gi, taxid, name) = line.split("\t")
            if wp[-2:] == ".1":
                wp = wp[:-2]
            if wp not in multi:
                multi[wp] = [(gi, taxid, name)]
            else:
                multi[wp].append((gi, taxid, name))

        return multi

def lookup_taxids(gis=[], taxids_db={}):
    """Get TaxIDs for a list of GIs
    """
    if not gis:
        print "No GIs"
        return

    if not taxids_db:
        print "No taxids loaded"
        return

    taxids = []
    for gi in gis:
        if gi in taxids_db:
            taxids.append(taxids_db[gi])
        else:
            taxids.append(None)

    return taxids

def find_self_taxid(self_gis={}, taxids_db={}, children_nodes={}, multispecies={}):
    """Find the lowest taxonomic unit that is shared across all self genes
    """

    if not self_gis:
        print "No self GIs"
        return

    if not children_nodes:
        print "No children nodes database"
        return

    gis = []
    ref_ids = []

    for gi in self_gis.keys():
        ref_id = self_gis[gi]
        if multispecies and ref_id in multispecies:
            ref_ids.append(ref_id)
        else:
            gis.append(gi)

    self_taxids = set(lookup_taxids(gis, taxids_db))
    
    if multispecies and ref_ids:
        for ref_id in ref_ids:
            ref_taxids = set()
            for (gi, taxid, name) in multispecies[ref_id]:
                ref_taxids.add(taxid)

            if not ref_taxids & self_taxids:
                # no common TaxIDs for this multi-species entry
                print "Found multispecies entry that does not contain other self TaxIDs"

    if len(self_taxids) == 1:
        # if only one self TaxID
        return self_taxids.pop()

    common = set()
    for taxid in self_taxids:
        temp = taxid
        parents = []
        while temp in children_nodes:
            parents.append(temp)
            if temp in children_nodes:
                temp = children_nodes[temp]
        if len(common) == 0:
            # initial
            common = set(parents)
        else:
            common = common.intersection(parents)

    if len(common) > 1:
        # this is so bad
        # returns the first TaxID in the last parents list
        # that is also in the common list
        # this will return the lowest taxon
        for taxid in parents:
            if taxid in common:
                return taxid
    else:
        # return the common TaxID
        return common.pop()

def determine_tax_groups(self_gis={}, taxids_db={}, children_nodes={}, multispecies={}):
    """Determine the taxonomic groupings for self and close
    """
    if not self_gis:
        print "No self GIs"
        return

    if not taxids_db:
        print "No taxids database"
        return

    if not children_nodes:
        print "No nodes database"
        return

    self = find_self_taxid(self_gis, taxids_db, children_nodes, multispecies)

    if not self:
        print "Could not determine self group"
        return (None, None)

    if self in children_nodes:
        close = children_nodes[self]
    else:
        print "Self group (%s) not in database" % self
        close = None


    print "Found self group TaxID: %s" % self
    print "Found close group TaxID: %s" % close

    return ([self], [close])

def _ignore_taxa_list_filter(taxid=None, ignore_list=[], names_db={}):
    """Searches for ignore list in TaxID name
    """

    if not taxid:
        return True
    if not ignore_list:
        return
    if not names_db:
        return

    if taxid in names_db:
        for word in ignore_list:
            if word.lower() in " ".join([x[0] for x in names_db[taxid]]).lower():
                return True
    
    if taxid in ignore_list:
        return True

def _filter_coverage(filename="", blast={}, info={}, min_coverage=70):
    """Filter a BLAST result file for sequence coverage
    
    This function does not remove subject entries from a BLAST result
    if you feed in a dictionary of blast results. It copies it yay!
    """
    
    if not filename:
        print "No filename"
        return {}
    
    if not info:
        print "No gene info"
        return {}
    
    # copy yay
    results = blast.copy()
    
    coverage = {}
    
    if results:
        for subject in sorted(results, key=lambda x: results[x]['lines'][0][0]):
            try:
                base_name = basename(filename)
                if base_name in info:
                    seq_length = info[base_name][1]
                
                if seq_length:
                    if len(results[subject]['lines']) > 1:
                        temp = [0 for x in xrange(seq_length)]
                        for i, start in enumerate(results[subject]['qstart']):
                            end = results[subject]['qend'][i]
                            for j in range(start-1, end):
                                if j >= 0 and j < len(temp):
                                    temp[j] = 1

                        # this is needed for overlapping BLAST hits in a query
                        fixed_aln_length = float(sum(temp))
                    else:
                        fixed_aln_length = float(results[subject]['qend'][0] - results[subject]['qstart'][0] + 1)
                    
                    p_coverage = fixed_aln_length / seq_length * 100.
                    coverage[subject] = p_coverage
                    if p_coverage < min_coverage:
                        del results[subject]
                
            except Exception as e:
                print e
    
    try:
        all_lines = []
        if results:
            for subject in results:
                all_lines.extend(results[subject]['lines'])
            
            all_lines = sorted(all_lines, key=lambda x: x[0])
        with open("%s.filtered" % filename, "wb") as w_fil:
            for i, line in all_lines:
                w_fil.write(line)
                
    except Exception as e:
        print e
    
    return coverage

def _parse_blast_file(file, pident_threshold=None, evalue_cutoff=None):
    """Parse a blast file into a dictionary
    """
    if not file:
        return ({}, 0)
    
    blast = {}
    removed = 0
    with open(file, 'rb') as f:
        for i, line in enumerate(f):
            try:
                temp = line.strip().split("\t")
                if len(temp) != 12:
                    removed += 1
                    continue

                pident = float(temp[2])
                if pident_threshold and pident < pident_threshold:
                    removed += 1
                    continue

                evalue = float(temp[-2])
                if evalue_cutoff and evalue > evalue_cutoff:
                    removed += 1
                    continue

                subject = temp[1]
                bitscore = float(temp[-1])
                aln_length = float(temp[3])
                qstart = int(temp[6])
                qend = int(temp[7])

                if subject not in blast:
                    blast[subject] = {}
                    blast[subject]['bitscore'] = []
                    blast[subject]['aln_length'] = []
                    blast[subject]['evalue'] = []
                    blast[subject]['pident'] = []
                    blast[subject]['qstart'] = []
                    blast[subject]['qend'] = []
                    blast[subject]['lines'] = []

                blast[subject]['bitscore'].append(bitscore)
                blast[subject]['aln_length'].append(aln_length)
                blast[subject]['evalue'].append(evalue)
                blast[subject]['pident'].append(pident)
                blast[subject]['qstart'].append(qstart)
                blast[subject]['qend'].append(qend)
                blast[subject]['lines'].append((i, line))
            except Exception as e:
                print e
    
    return (blast, removed)

def parse_blast_result(file, config={}, taxids_db={}, names_db={}, multispecies={}, info={}, nucl=False, n_hits=None, top_percent=None):
    """Parse a BLAST result file to include TaxIDs, also expanding if there are

    Fixed for having multiple alignments in BLAST result
    """
    if not file:
        print "No result file"
        return
    if not taxids_db:
        print "No TaxIDs database"
        return
    if not names_db:
        print "No names database"
        return
    # if not multispecies:
        # print "No multispecies database"
        # return
    if not info:
        print "No gene information"
        return

    if config and "percIdent" in config:
        pident_threshold = float(config["percIdent"])
    else:
        pident_threshold = None

    if config and "evalue" in config:
        evalue_cutoff = float(config["evalue"])
    else:
        evalue_cutoff = None

    print "Parsing BLAST results in %s" % os.path.basename(file)
    (blast, removed) = _parse_blast_file(file, pident_threshold, evalue_cutoff)
    
    if nucl:
        min_coverage = 80
    else:
        min_coverage = 70
    coverage = _filter_coverage(file, blast, info, min_coverage)
    
    results = {}
    with open("%s.out" % file, "wb") as w:
        w.write("sseqid\ttaxid\tpident\tcoverage\tevalue\tscore\n")
        if blast:
            
            # sort by highest to lowest score = bitscore / aln_length
            if top_percent:
                filter_score = None
            
            for i, subject in enumerate(sorted(blast.keys(), key=lambda x: sum(blast[x]['bitscore']) / sum(blast[x]['aln_length']), reverse=True)):
                try:
                    if n_hits and i >= n_hits:
                        break
                    
                    total_bitscore = sum(blast[subject]['bitscore'])
                    total_aln_length = (sum(blast[subject]['aln_length']))
                    score = total_bitscore / total_aln_length / 2
                    evalue = max(blast[subject]['evalue'])
                    pident = 0
                    for i, p in enumerate(blast[subject]['pident']):
                        pident += p * blast[subject]['aln_length'][i]
                    pident /= total_aln_length
                    
                    if top_percent:
                        if i == 0:
                            filter_score = score * (1 - (top_percent / 100.))
                        else:
                            if filter_score and score < filter_score:
                                break
                    
                    p_coverage = coverage.get(subject, 0)
                    if p_coverage < min_coverage:
                        removed += 1
                        # percent coverage < 70%
                        # print "Removed %s due to < 70 percent coverage" % subject
                        # del blast[subject]
                        continue
                    
                    (sub_gi, sub_ref) = get_gene_name(subject)
                    
                    if not nucl and multispecies and sub_ref in multispecies:
                        # remove old lines so you don't keep the multispecies ones
                        # lines = blast[subject]['lines'][:]
                        # blast[subject]['lines'] = []
                        for j, (gi, taxid, name) in enumerate(multispecies[sub_ref]):
                            if config and "ignoreTaxa" in config and _ignore_taxa_list_filter(taxid, config["ignoreTaxa"], names_db):
                                # removed += 1
                                continue

                            w.write("%s\t%s\t%.2f\t%.2f\t%g\t%f\n" % (gi, taxid, pident, p_coverage, evalue, score))
                            
                            if gi in results:
                                #print "Duplicate GI found %s" % gi
                                gi += '%i.%i' % (i, j)
                            
                            results[gi] = (taxid, score)
                    
                    else:
                        if sub_gi in taxids_db:
                            taxid = taxids_db[sub_gi]
                            if config and "ignoreTaxa" in config and _ignore_taxa_list_filter(taxid, config["ignoreTaxa"], names_db):
                                removed += 1
                                # del blast[subject]
                                continue
                            w.write("%s\t%s\t%.2f\t%.2f\t%g\t%f\n" % (sub_gi, taxid, pident, p_coverage, evalue, score))
                            if sub_gi in results:
                                sub_gi += ".%i" % i
                            results[sub_gi] = (taxid, score)

                except Exception as e:
                    print e
    
    if blast:
        print "Filtered out %i BLAST hits from %s" % (removed, os.path.basename(file))
    else:   
        print "No BLAST results for %s" % file
        
    return results

def load_results_file(file, config={}, lines=None):
    """Load previously parsed results from a text file
    """
    if not file:
        print "No result file"
        return
    
    if config and "percIdent" in config:
        pident_threshold = float(config["percIdent"])
    else:
        pident_threshold = None

    if config and "evalue" in config:
        evalue_cutoff = float(config["evalue"])
    else:
        evalue_cutoff = None
    
    results = {}
    with open(file, 'rb') as f:
        f.readline() # don't need the header, just skip it
        for i, line in enumerate(f):
            if lines and i >= lines:
                    break
            temp = line.strip().split('\t')
            pident = float(temp[2])
            p_coverage = float(temp[3])
            evalue = float(temp[4])
            if pident_threshold and pident < pident_threshold:
                # pident too low
                continue
            if evalue_cutoff and evalue > evalue_cutoff:
                # evalue too high
                continue

            # key is GI, value is tuple of (taxid, score)
            if temp[0] in results:
                temp[0] += ".%i" % i
            results[temp[0]] = (temp[1], float(temp[-1]))

    return results

def expand_taxids(taxids=[], parent_nodes={}):
    """Expand TaxIDs to a list of it and its children
    """
    expanded = []
    for taxid in taxids:
        if taxid in parent_nodes:
            expanded.append(taxid)
            expanded.extend(expand_taxids(parent_nodes[taxid], parent_nodes))
        else:
            expanded.append(taxid)

    return expanded

def get_weights(blast_results={}, self_ids={}, config={}, taxids_db={}, children_nodes={}, parent_nodes={}, multispecies={}):
    """Get the self, close, and distal weights of genes
    """
    if not blast_results:
        print "No BLAST results"
        return
    if not config:
        config = DEFAULT_CONFIG
    if not taxids_db:
        print "No TaxIDs database"
        return
    if not children_nodes:
        print "No children nodes database"
        return
    if not parent_nodes:
        print "No parent nodes database"
        return
    # if not multispecies:
        # print "No multispecies database"
        # return

    if config['selfGroup']:
        self = config['selfGroup']
        if type(self) != list:
            self = [self]
        self = set(self)
        if config['closeGroup']:
            close = (config['closeGroup'])
            if type(close) != list:
                close = [close]
            close = set(close)
        else:
            print "Please specify close group if specifying self group!"
            return

    else:
        (self, close) = determine_tax_groups(self_ids, taxids_db, children_nodes, multispecies)
        if not (self and close):
            print "Could not determine self & close groups"
            return
        self = set(self)
        close = set(close)

    self = set(expand_taxids(self, parent_nodes))
    close = set(expand_taxids(close, parent_nodes))
    
    weights = {}
    for gene in blast_results:
        self_weight = 0.
        close_weight = 0.
        distal_weight = 0.
        for gi in blast_results[gene]:
            (taxid, score) = blast_results[gene][gi]
            if self != close:
                if taxid in self:
                    self_weight += score
                elif taxid in close:
                    close_weight += score
                else:
                    distal_weight += score
            else:
                # ignore self
                if taxid in close:
                    close_weight += score
                else:
                    distal_weight += score

        weights[gene] = (self_weight, close_weight, distal_weight)

    return weights, self, close

def _get_cutoff(values, relaxed=False):
    """Use R density and pastecs to determine cutoff values
    """
    if not values:
        print "No values"
        return

    try:
        rvals = com.convert_to_r_dataframe(pandas.DataFrame(values))
        ro.globalenv['values'] = rvals
        ro.r('smooth <- density(t(values))')
        ro.r('library(pastecs)')
        ro.r('tps <- turnpoints(smooth$y)')
        ro.r('peaks <- smooth$x[tps$peaks]')
        ro.r('pits <- smooth$x[tps$pits]')

        if relaxed:
            ro.r('cutoff <- pits[1]')
        else:
            ro.r('cutoff <- median(c(peaks[1], pits[1]))')

        cutoff = com.load_data('cutoff')[0]
        if cutoff != None and cutoff >= 0:
            return cutoff

    except Exception as e:
        print e

    print "Could not use R to determine cutoff."
    cutoff = numpy.percentile(values, 25)
    print "Using 25th percentile instead: %f" % cutoff
    return cutoff

def get_cutoffs(weights={}, self_relaxed=False, close_relaxed=False, distal_relaxed=False):
    """Return cutoff weights for each group
    """

    self = []
    close = []
    distal = []

    for (x, y, z) in weights.values():
        self.append(x)
        close.append(y)
        distal.append(z)

    self_cutoff = _get_cutoff(self, self_relaxed)
    close_cutoff = _get_cutoff(close, close_relaxed)
    distal_cutoff = _get_cutoff(distal, distal_relaxed)

    return (self_cutoff, close_cutoff, distal_cutoff)

def assign_hgt(weights={}, self=None, close=None, distal=None, useSelf=False):
    """Assign whether genes are HGT derived or not based on cutoffs of weights
    """

    if not weights:
        print "No weights"
        return

    if useSelf and self == None:
        print "No self cutoff, but useSelf set to True"
        return
    if close == None:
        print "No close cutoff"
        return
    if distal == None:
        print "No distal cutoff"
        return

    hgt_results = {}

    # if useSelf, be below self cutoff
    # be below close cutoff
    # be above distal cutoff
    # then you are transferred
    for gene in weights:
        hgt_results[gene] = ((not useSelf or weights[gene][0] < self) and weights[gene][1] < close and weights[gene][2] > distal)

    return hgt_results

def _get_sci_name(taxid=None, names_db={}):
    """Return scientific name of a TaxID
    """

    if not taxid:
        return ''

    if not names_db:
        return ''

    if taxid in names_db:
        for (name, tax_class) in names_db[taxid]:
            if tax_class and tax_class == 'scientific name':
                return name

        return names_db[taxid][0][0]


    return ''

def _get_best_match(blast_result={}, self=[]):
    """Returns the best BLAST match and its TaxID
    """
    if not blast_result:
        print "No blast result"
        return

    non_self = {}
    for gi in blast_result:
        taxid = blast_result[gi][0]
        if taxid not in self:
            non_self[gi] = blast_result[gi]

    if not non_self:
        #print "No non-self hits"
        return '', ''

    best = max(non_self.keys(), key=lambda x: non_self[x][1])
    return (best, non_self[best][0])

def _write_blast_matches(file, blast_result={}, self=[], names_db={}):
    """Write non-self blast matches in order of highest score
    """
    if not file:
        print "No file to write"
        return

    if not blast_result:
        print "No blast result"
        return

    if not names_db:
        print "No names database"
        return

    non_self = {}
    for gi in blast_result:
        taxid = blast_result[gi][0]
        if taxid not in self:
            non_self[gi] = blast_result[gi]

    if not non_self:
        print "No non-self hits"
        return

    sorted_seqs = sorted(non_self.keys(), key=lambda x: non_self[x][1], reverse=True)

    with open(file, 'wb') as w:
        w.write('Sequence\tTaxID\tTaxa\tScore\n')
        for seq in sorted_seqs:
            (taxid, score) = non_self[seq]
            w.write('%s\t%s\t%s\t%f\n' % (seq, taxid, _get_sci_name(taxid, names_db), score))


def write_results(file, info={}, blast_results={}, weights={}, hgt_results={}, names_db={}, self_taxids=[], write_matches=False, close_taxids=[]):
    """Write results from pipeline to a file
    """

    if not file:
        print "No file"
        return

    if not weights:
        print "No weights to write"
        return

    if not blast_results:
        print "No blast results"
        return

    if not hgt_results:
        print "No HGT results"
        return

    if not names_db:
        print "No names database"
        return

    genes = sorted(set(weights.keys()) & set(blast_results.keys()) & set(hgt_results.keys()))

    if write_matches:
        match_dir = os.path.join('result', 'matches')
        if not os.path.isdir(match_dir):
            os.mkdir(match_dir)

    with open(file, 'wb') as w:
        w.write("Query\tLength\tProduct\tHits\tSelf\tClose\tDistal\tHGT\tRecent\tMatchTaxID\tMatchTaxName\tMatchSeqID\n")
        for query in genes:
            try:
                (product, length) = info.get(query, ('', 0))
                length = "%i" % length
                hits = "%i" % len(blast_results.get(query, 0))
                (self, close, distal) = ["%.3f" % x for x in weights.get(query, (0,0,0))]
                recent = '0'
                if hgt_results.get(query):
                    hgt = '1'
                    if write_matches:
                        match_file = os.path.join(match_dir, "%s.txt" % query)
                        _write_blast_matches(match_file, blast_results[query], self_taxids, names_db)
                else:
                    hgt = '0'
                if blast_results.get(query):
                    (match_seq, match_tax) = _get_best_match(blast_results[query], self_taxids)
                    match_tax_name = _get_sci_name(match_tax, names_db)
                    # fix if there were duplicate GIs
                    if '.' in match_seq:
                        match_seq = match_seq.split('.')[0]
                    if hgt == '1' and match_tax and match_tax not in close_taxids:
                        recent = '1'                        
                    
                else:
                    match_seq = ""
                    match_tax = ""
                    match_tax_name = ""
                w.write("\t".join([query, length, product, hits, self, close, distal, hgt, recent, match_tax, match_tax_name, match_seq])+"\n")
            except Exception as e:
                print e
                break

def _load_dbs(taxdump=TAXDUMP_DIRECTORY, multispecies_file=MULTISPECIES_FILE, nucl=False):

    if not taxdump:
        taxdump = TAXDUMP_DIRECTORY

    if not multispecies_file:
        multispecies_file = MULTISPECIES_FILE

    print "Loading TaxIDs...",
    taxids_db = load_taxids(taxdump, nucl)
    print "Done."

    print "Loading children and parent nodes...",
    (parent_nodes, children_nodes) = load_nodes_file(taxdump)
    print "Done."

    print "Loading names database...",
    names_db = load_names_db(taxdump)
    print "Done."
    
    if nucl:
        print "Not loading multispecies file, since nucleotide sequences."
        multispecies = {}
    else:
        print "Loading multispecies database...",
        multispecies = read_multispecies(multispecies_file)
        print "Done."

    return taxids_db, parent_nodes, children_nodes, names_db, multispecies

def pipeline(directory=None, override_config={}, multispecies_file=MULTISPECIES_FILE, write_matches=False, nucl=False):
    """Main pipeline for running HGTector method
    """

    try:
        sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

        if not directory:
            directory = os.path.curdir
        else:
            os.chdir(directory)

        if not os.path.isdir("input"):
            print "No input folder"
            return

        if os.path.isfile("config.txt"):
            config = read_config_file("config.txt")
        else:
            print "No config.txt file. Using default parameters."
            config = DEFAULT_CONFIG

        # if script was ran with some override parameters
        if override_config:
            config.update(override_config)

        if int(config['modKCO']):
            self_relaxed = False
            close_relaxed = False
            distal_relaxed = False
        else:
            self_relaxed = True
            close_relaxed = True
            distal_relaxed = True

        # new configuration variables

        # modKCOSelf allows user to specify self cutoff
        if config['modKCOSelf'] != None:
            if int(config['modKCOSelf']):
                self_relaxed = False
            else:
                self_relaxed = True

        # modKCOClose allows user to specify close cutoff
        if config['modKCOClose'] != None:
            if int(config['modKCOClose']):
                close_relaxed = False
            else:
                close_relaxed = True

        # modKCODistal allows user to specify distal cutoff
        if config['modKCODistal'] != None:
            if int(config['modKCODistal']):
                distal_relaxed = False
            else:
                distal_relaxed = True
        
        if config['topPercent']:
            try:
                top_percent = float(config['topPercent'])
            except:
                top_percent = None
        else:
            top_percent = None
        
        
        if int(config['selfLow']):
            useSelf = True
        else:
            useSelf = False

        if int(config['interactive']):
            interactive = True
        else:
            interactive = False
                
        # split fasta files into sequences
        input_files = sorted(glob.glob(os.path.join("input", "*")))
        info = {}
        self_ids = {}

        if len(glob.glob(os.path.join('blast', '*', '*.fa'))):
            if interactive:
                print "It appears that there are already fasta sequence files"
                answer = raw_input("Do you want to copy sequences anyway? ")
                if 'y' in answer:
                    print "Will copy all sequences"
                    write_sequences = True
                elif 'n' in answer:
                    print "Will not copy sequences"
                    write_sequences = False
            else:
                print "Found fasta files already in blast folder. Will not copy over sequences."
                write_sequences = False
        else:
            write_sequences = True

        for file in input_files:
            info[basename(file)], self_ids[basename(file)] = split_fasta(file, write_sequences)

        folders = sorted([dir for dir in glob.glob(os.path.join("blast", "*")) if os.path.isdir(dir)])

        # run blast on all sequences
        blastcmds = []

        blast_db = config["dbBLAST"]
        n_hits = int(config["nHits"])
        evalue = float(config["evalue"])
        for folder in folders:
            blastcmds.extend(get_blast_cmds(folder, blast_db=blast_db, n_hits=n_hits, evalue=evalue, nucl=nucl))

        if len(blastcmds):
            print "Running BLAST on %i sequences..." % len(blastcmds)
            thread_cmds(blastcmds, int(config["nThreads"]))
        elif glob.glob(os.path.join("blast", "*", "*.bla")):
            if interactive:
                print "There are previously created blast results."
                answer = raw_input("Do you wish to continue? ")
                if 'y' in answer:
                    print "Continuing..."
                elif 'n' in answer:
                    print "Okay, quitting"
                    return
            print "Using previously created blast results"
        else:
            print "Nothing to BLAST. No previous BLAST results."
            return
        
        print "Loading databases..."
        
        taxids_db, parent_nodes, children_nodes, names_db, multispecies = _load_dbs(config["taxdump"], multispecies_file, nucl)
        
        blast_results = {}

        # parse BLAST results
        for folder in folders:
            blast_results[folder] = {}
            
            if len(glob.glob(os.path.join(folder, '*.bla.out'))) == len(glob.glob(os.path.join(folder, "*.bla"))):
                print "Found previous output files"
                if interactive:
                    answer = raw_input("Do you wish to use these files? ")
                    if 'y' in answer:
                        print "Okay. Using these files."
                        analyze_blast = False
                    else:
                        print "Okay. Will redo analysis."
                        analyze_blast = True
                else:
                    print "Will use previously analyzed output files."
                    analyze_blast = False
            elif glob.glob('*.bla.out'):
                print "Found some previous output files, but incomplete."
                print "Will redo analysis"
                analyze_blast = True
            else:
                analyze_blast = True
            
            if analyze_blast:
                print "Parsing BLAST results for %s..." % folder
                for file in sorted(glob.glob(os.path.join(folder, "*.bla")), key=lambda x: os.path.getsize(x) if os.path.isfile(x) else -1):
                    if os.path.isfile(file):
                        blast_results[folder][basename(file)] = parse_blast_result(file, config, taxids_db, names_db, multispecies, info[basename(folder)], nucl, n_hits=n_hits, top_percent=top_percent)
                    else:
                        print "Warning: %s no longer there!" % file
                
            else:
                print "Loading previously analyzed results in %s..." % folder
                for file in sorted(glob.glob(os.path.join(folder, '*.bla.out')), key=lambda x: os.path.getsize(x) if os.path.isfile(x) else -1):
                    if os.path.isfile(file):
                        blast_results[folder][os.path.basename(file).replace('.bla.out', '')] = load_results_file(file, config, lines=n_hits)
                    else:
                        print "Warning: %s no longer there!" % file
            
            print "Done. Found %i BLAST results" % len(blast_results[folder])


        if not os.path.isdir('result'):
            os.mkdir('result')

        for folder in blast_results:
            if not blast_results[folder]:
                print "No BLAST results for %s" % folder
                continue

            input_file = basename(folder)

            # calculate weights
            print "Calculating weights for %s..." % folder,
            weights, self_taxids, close_taxids = get_weights(blast_results[folder], self_ids[input_file], config, taxids_db, children_nodes, parent_nodes, multispecies)
            print "Done."

            if not weights:
                print "No weights calculated for %s" % folder
                continue

            # get the cutoff values for each weight
            print "Obtaining cutoffs for each taxonomic group...",

            (self, close, distal) = get_cutoffs(weights, self_relaxed, close_relaxed, distal_relaxed)
            print "Done."
            print "Self cutoff: %g" % self
            if interactive:
                answer = raw_input('Is this value okay to use? ')
                if 'y' in answer:
                    pass
                else:
                    try:
                        self = float(raw_input('Enter a value for self: '))
                    except:
                        pass
            print "Close cutoff: %g" % close
            if interactive:
                answer = raw_input('Is this value okay to use? ')
                if 'y' in answer:
                    pass
                else:
                    try:
                        close = float(raw_input('Enter a value for close: '))
                    except:
                        pass
            print "Distal cutoff: %g" % distal
            if interactive:
                answer = raw_input('Is this value okay to use? ')
                if 'y' in answer:
                    pass
                else:
                    try:
                        close = float(raw_input('Enter a value for distal: '))
                    except:
                        pass

            # assign HGT or not to each gene based on cutoff of weights
            print "Assigning HGT to genes...",
            hgt_results = assign_hgt(weights, self, close, distal, useSelf)
            print "Done."

            # write results to a text file
            fileout = os.path.join('result', '%s.txt' % input_file)
            print "Writing results to %s..." % fileout

            write_results(fileout, info[input_file], blast_results[folder], weights, hgt_results, names_db, self_taxids, write_matches, close_taxids)
            print "Done. %s is finished. Find results in the result folder." % input_file

        print "HGTector method complete."

    except Exception as e:
        print e
        return

def main():
    parser = OptionParser(usage='Usage: %prog directory [...]')
    parser.add_option('--blastn', dest='blastn', action='store_true', default=False,
                      help='Perform blastn on DNA sequences instead of blastp')
    parser.add_option('-n', '--nhits', dest='nhits', type='int',
                      help='Override nHits in config.txt')
    parser.add_option('-e', '--evalue', dest='evalue', type='float',
                      help='Override evalue in config.txt')
    parser.add_option('-r', '--relaxed', dest='relaxed', action='store_true', default=False,
                      help='Use relaxed cutoff')
    parser.add_option('-c', '--conservative', dest='conservative', action='store_true', default=False,
                      help='Use conservative cutoff')
    parser.add_option('-t', '--threads', dest='threads', type='int',
                      help='Override nThreads in config.txt')
    parser.add_option('--multispecies', dest='multispecies', default=MULTISPECIES_FILE,
                      help='Specify RefSeq multispecies definition file')
    parser.add_option('--db', dest='db', default=BLAST_DB,
                      help='Specify different nr database to use')
    parser.add_option('--taxdump', dest='taxdump', default=TAXDUMP_DIRECTORY,
                      help='Specify different taxdmp directory to use')
    parser.add_option('-w', '--writematches', dest='writematches', action='store_true', default=False,
                      help='Write BLAST matches for detected HGT genes to a text file')
    (options, args) = parser.parse_args()

    override_config = {}

    if options.nhits:
        print "Will use %i hits instead of value in config.txt" % options.nhits
        override_config['nHits'] = options.nhits

    if options.evalue:
        print "Will use %f as evalue cutoff instead of config.txt" % options.evalue
        override_config['evalue'] = options.evalue

    if options.relaxed and options.conservative:
        print "Specified both cutoffs. Please specify only one"
        return

    if options.relaxed:
        print "Will use relaxed cutoff"
        override_config['modKCO'] = 0

    if options.conservative:
        print "Will use conservative cutoff"
        override_config['modKCO'] = 1

    if options.threads:
        print "Will use %i threads instead of value from config.txt" % options.threads
        override_config['nThreads'] = options.threads

    if options.multispecies and os.path.isfile(options.multispecies):
        print "Will use %s as the definitions for multispecies RefSeq entries" % options.multispecies
    else:
        print "Will use default multispecies file"
        options.multispecies = MULTISPECIES_FILE

    if options.db and options.db != BLAST_DB:
        print "Will use %s as BLAST database" % options.db
        override_config['dbBLAST'] = options.db

    if options.taxdump and options.taxdump != TAXDUMP_DIRECTORY:
        print "Will use directory %s for taxonomic data"% options.taxdump
        override_config['taxdump'] = options.taxdump

    if args:
        # do for each specified directory
        for dir in args:
            pipeline(dir, override_config=override_config, multispecies_file=options.multispecies, write_matches=options.writematches, nucl=options.blastn)

    else:
        # just do for current directory
        pipeline(override_config=override_config, multispecies_file=options.multispecies, write_matches=options.writematches, nucl=options.blastn)

if __name__ == "__main__":
    main()
