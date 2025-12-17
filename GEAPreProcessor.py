import copy
import os
import shutil
import time
import xml.dom.minidom
from os.path import exists
import pandas as pd
from franklin_model.experiment_model.python_experiment_model import Context, User, Institution, Comparison

from experiment.RosalindExperiment import known_orgs, ensembl_genomes, read_exp_meta, RosalindExperiment
# from GEODatasetsImportController import clean_title
# from LocalConfig import gea_dir
# from ArrayDataPreprocessor import process_array_dataframe
from experiment.geo.ArrayDataPreprocessor import process_array_dataframe
from experiment.geo.GEODatasetsImportController import clean_title
from LocalConfig import ROOT_DIR
from experiment.util.converter_utils import get_immediate_subdirectories, get_files
from strings2things.TextMatcher import TextMatcher
gea_dir = os.path.join(ROOT_DIR, 'data', 'experiments', 'raw', 'gea')
"""
This script contains methods for converting experiments in the Gene Expression Atlas (GEA)
into experiments ready for upload to the Rosalind platform.  GEA experiment data may be found in:
http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/   or in 
ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/atlas-latest-data.tar.gz 
The latter is a single giant archive that causes most methods to fail when unzipping it. 
(Used https://www.urbanophile.com/arenn/coding/gzrt/ to get through opening it.) 
"""


def make_attributes_file_from_condensed_sdrf(sdrf_file, outfile, include_uris=False):
    """
    The EBI uses the SDRF file format to describe experimental metadata for ArrayExpress, Gene Expression Atlas etc.
    https://www.ebi.ac.uk/arrayexpress/help/creating_a_sdrf.html
    This method extracts a simple tab-delimited 'attributes file' suitable for loading into Rosalind.
    """
    orgs = set()
    # sdrf_list, header, header_dict = read_sdrf_file(sdrf_file)
    sdrf_headers = ['Accession', 'Array', 'Sample', 'Annot_type', 'Annot', 'Annot_value', 'Annot_ont_URI']
    df = pd.read_csv(sdrf_file, sep="\t", names=sdrf_headers)
    # make header row - including all factors we want
    character_headers = sorted(df.loc[df['Annot_type'] == 'characteristic'].Annot.unique())
    factor_headers = sorted(df.loc[df['Annot_type'] == 'factor'].Annot.unique())
    sample_ids = sorted(df.Sample.unique())

    hrow = 'Sample'
    for factor in factor_headers:
        hrow = hrow + "\tfactor:" + factor
        if include_uris:
            hrow = hrow + "\t" + factor + "_factor_URI"
    for character in character_headers:
        hrow = hrow + "\tcharacteristic:" + character
        if include_uris:
            hrow = hrow + "\t" + character + "_character_URI"
    with open(outfile, "w") as out:
        out.write(hrow + "\n")
    for sample in sample_ids:
        # clean sample ids of certain characters
        # entries may only contain letters numbers spaces and dashes
        row = clean_title(sample) + "\t"
        for factor in factor_headers:
            factor_values = df.loc[(df['Sample'] == sample) & (df['Annot'] == factor)].Annot_value.unique()
            if len(factor_values) == 0:
                factor_text_value = 'na'
            else:
                factor_text_value = str(factor_values[0])
                factor_text_value = clean_title(factor_text_value)
            if factor == 'organism':
                orgs.update(factor_values)
            if include_uris:
                factor_uri_values = df.loc[(df['Sample'] == sample) & (df['Annot'] == factor)].Annot_ont_URI.unique()
                factor_uri_value = 'na'
                if len(factor_uri_values) == 1:
                    if str(factor_uri_values[0]) != 'nan':
                        factor_uri_value = str(factor_uri_values[0])
                row = row + factor_text_value + "\t" + factor_uri_value + "\t"
            else:
                row = row + factor_text_value + "\t"
        for character in character_headers:
            character_values = df.loc[(df['Sample'] == sample) & (df['Annot'] == character)].Annot_value.unique()
            if len(character_values) == 0:
                character_text_value = 'na'
            else:
                character_text_value = str(character_values[0])
                character_text_value = clean_title(character_text_value)
            if character == 'organism':
                orgs.update(character_values)
            if include_uris:
                character_uri_values = df.loc[
                    (df['Sample'] == sample) & (df['Annot'] == character)].Annot_ont_URI.unique()
                character_uri_value = 'na'
                if len(character_uri_values) == 1:
                    if str(character_uri_values[0]) != 'nan':
                        character_uri_value = str(character_uri_values[0])
                row = row + character_text_value + "\t" + character_uri_value + "\t"
            else:
                row = row + character_text_value + "\t"
        row = row.rstrip()
        row = row.rstrip("\t")
        with open(outfile, "a") as out:
            out.write(row + "\n")
    return True, orgs


def prepare_data(is_raw=False, inputfile='', outputfile='', exp_id='', is_ensembl=False):
    '''
    This method reads a tab-separated microarray file and gets it ready for upload to Rosalind.
    Most of the work is done by ArrayDataPreprocessor.process_array_dataframe
    '''
    data = pd.read_csv(inputfile, sep="\t")
    data, worked = process_array_dataframe(data, genekey='Gene ID', is_raw_count=is_raw, exp_id=exp_id)
    if not worked:
        return False, is_ensembl
    test_id1 = data.iloc[0, 0]
    test_id2 = data.iloc[1, 0]
    if not is_ensembl:
        if (test_id1.startswith('ENS') or test_id1.startswith('AT') or test_id1.startswith('FB') or test_id1.startswith('WB')) \
            and (test_id2.startswith('ENS') or test_id2.startswith('AT') or test_id2.startswith('FB') or test_id2.startswith('WB')):
            is_ensembl = True
        else:
            print("really not ensembl")
    colnames = list(data.columns.values)
    # todo make sure ensembl is being detected properly!!!!
    if is_ensembl and 'Gene ID' in colnames and 'Gene Name' in colnames:
        data.drop(['Gene Name', ], axis=1, inplace=True)
    #TODO heavily dependent on their 'gene name' matching internal 'symbol'.  Would be nice to check.
    elif not is_ensembl and 'Gene ID' in colnames and 'Gene Name' in colnames:
        data.drop(['Gene ID', ], axis=1, inplace=True)
        data.rename(columns={'Gene Name': 'symbol'}, inplace=True)
    if 'DesignElementAccession' in colnames:
        data.drop(['DesignElementAccession', ], axis=1, inplace=True)
    #clean sample ids of certain characters
    #entries may only contain letters numbers spaces and dashes
    data.columns = [clean_title(s) for s in data.columns]
    data.to_csv(outputfile, sep='\t', index=False)
    return True, is_ensembl


def get_comparisons(config_xml):
    '''
    The EBI uses a standardized xml file to describe sample groupings for experiments.
    This method converts that xml file into the equivalent json file in the Rosalind platform.
    '''
    doc = xml.dom.minidom.parse(config_xml)
    # get assay groups
    assay_groups = doc.getElementsByTagName("assay_group")
    groups = {}
    comparisons = []
    for assay_group in assay_groups:
        if assay_group.hasAttribute("id"):
            gid = assay_group.getAttribute("id")
            if assay_group.hasAttribute("label"):
                name = assay_group.getAttribute("label")
            else:
                name = gid
            sample_names = set()
            assays = assay_group.getElementsByTagName("assay")
            for assay in assays:
                sample_names.add(clean_title(assay.childNodes[0].data))
            g = Context(name, is_control=False, sample_names=list(sample_names), local_id=gid, category="Context")
            groups[gid] = g
    for contrast in doc.getElementsByTagName("contrast"):
        if contrast.hasAttribute("id"):
            comparison_id = contrast.getAttribute("id")
            comparison_name = contrast.getElementsByTagName("name")[0].childNodes[0].data
            control_group_id = contrast.getElementsByTagName("reference_assay_group")[0].childNodes[0].data
            test_group_id = contrast.getElementsByTagName("test_assay_group")[0].childNodes[0].data
            # copy sometimes one comparison uses a group as a control while another uses it as a the test
            control_group = copy.deepcopy(groups[control_group_id])
            control_group.isControl = True
            test_group = copy.deepcopy(groups[test_group_id])
            test_group.isControl = False
            c = Comparison(name=comparison_name, id=comparison_id, category="Comparison")
            #groups = [control_group, test_group],
            c.has_control_context = control_group
            c.has_test_context = test_group
            comparisons.append(c)
    return comparisons


def ensembl_to_ncbi(ensembl_gene, gene_matcher):
    matches = gene_matcher.match_text("ENSEMBL:"+str(ensembl_gene))
    if len(matches) > 0:
        return matches[0]
    else:
        return None


def prepare_experiment_for_load(gea_acc, gea_experiments_dir, rosalind_experiments_dir, include_att_uris=True,
                                errorfile=None, gene_matcher=None):
    """
    Prepare all of the data files needed to load the experiment into Rosalind.
    For GEA, this is mostly finding the right files from their export, parsing and reformatting, and then
    doing some mild remapping and bioinformatics checks on the array data files.
    """
    exp = RosalindExperiment(experiment_id="exp:"+str(gea_acc), public_id=gea_acc)
    exp.user = User(id="user:gea", name="GEA", category="User")
    # start with metadata
    # get the experiment name and description
    gea_exp_folder = os.path.join(gea_experiments_dir, gea_acc)
    rosalind_exp_folder = os.path.join(rosalind_experiments_dir, gea_acc)
    if not os.path.exists(rosalind_exp_folder):
        os.makedirs(rosalind_exp_folder)
    idf_file = os.path.join(gea_exp_folder, gea_acc + '.idf.txt')
    if not exists(idf_file):
        exp.log_message("gea metadata parse error. couldn't find idf file", errorfile)
        return None
    with open(idf_file, "r") as idf:
        for line in idf:
            row = line.split("\t")
            if row[0] == 'Investigation Title':
                exp.experiment.name = row[1]
                exp.experiment.name = clean_title(exp.experiment.name)
            elif row[0] == 'Experiment Description':
                exp.experiment.description = row[1]
            elif row[0] == 'Person Email':
                if row[1] != 'geo@ncbi.nlm.nih.gov' and row[1] != '':
                    exp.user = User(id="gea_user:"+str(row[1]), email=row[1], category="User")
            elif row[0] == 'Person Last Name':
                exp.user.last_name = row[1]
            elif row[0] == 'Person First Name':
                exp.user.first_name = row[1]
            elif row[0] == 'Person Affiliation':
                exp.institution = Institution(id="institution:"+str(row[1].replace(" ", "_")), name=row[1], category="Institution")
            elif row[0] == 'Publication Title':
                exp.publication_title = row[1]
            elif row[0] == 'PubMed ID':
                exp.pubmed_id = row[1]
        if exp.institution and exp.user:
            exp.user.member_of = [exp.institution]
        exp.experiment.created_by = exp.user
        if not exp.experiment.name or not exp.experiment.description:
            exp.log_message("gea metadata parse error. couldn't get name from idf file", errorfile)
            return None

    # make the attributes file
    sdrf = os.path.join(gea_experiments_dir, gea_acc, gea_acc + '.condensed-sdrf.tsv')
    exp.experiment.attribute_file_path = os.path.join(rosalind_exp_folder, gea_acc + "-attributes.tsv")
            # folder = gea_data_for_import + gea_acc + '/'
    if not exists(sdrf):
        exp.log_message("SDRF missing: "+sdrf, errorfile)
        return None
    try:
        atts_made, orgs = make_attributes_file_from_condensed_sdrf(sdrf, exp.experiment.attribute_file_path, include_uris=include_att_uris)
    except:
        exp.log_message("gea metadata parse error. couldn't make attributes from sdrf file", errorfile)
        return None
    if not atts_made:
        exp.log_message("gea metadata parse error. couldn't make attributes from sdrf file", errorfile)
        return None
    if not orgs or len(orgs) != 1:
        exp.log_message("gea metadata parse error. couldn't couldn't parse organism from sdrf file", errorfile)
        return None
    organism = orgs.pop()
    if organism not in known_orgs:
        exp.log_message("gea metadata parse error. organism unknown: "+organism, errorfile)
        return None
    exp.experiment.genome = known_orgs[organism]
    exp.experiment.is_ensembl = False
    if exp.experiment.genome in ensembl_genomes:
        exp.experiment.is_ensembl = True

    # make the data file
    exp.experiment.data_file_path = None
    is_raw = False
    # find the source data..
    source_datafile = None
    for f in os.listdir(gea_exp_folder):
        f = str(f)
        if f.endswith('normalized-expressions.tsv') or f.endswith('normalized-expression.tsv'):
            source_datafile = os.path.join(gea_exp_folder, f)
            is_raw = False
            exp.experiment.exp_type = "RNA-counts-normalized"
            break
                #to use undecorated, need to actually work with ensembl id
                # elif f.endswith('normalized-expressions.tsv.undecorated') or f.endswith('normalized-expression.tsv.undecorated'):
                #     datafile = folder + f
                #     is_raw = False
                #     exp_type = "RNA-counts-normalized"
                #     break
        elif f.endswith('raw-counts.tsv'):
            source_datafile = os.path.join(gea_exp_folder, f)
            is_raw = True
            exp.experiment.exp_type = 'RNA-counts-raw'
            break
                # elif f.endswith('raw-counts.tsv.undecorated'):
                #     datafile = folder + f
                #     is_raw = True
                #     exp_type = 'RNA-counts-raw'
                #     break
    if not source_datafile:
        exp.log_message("gea parse error. couldn't find source data file ", errorfile)
        return None

    exp.experiment.data_file_path = os.path.join(rosalind_exp_folder, exp.experiment.public_id + "-data.tsv")
    print("preparing data for ", exp.experiment.public_id)
    data_ready, is_really_ensembl = prepare_data(is_raw=is_raw, inputfile=source_datafile, outputfile=exp.experiment.data_file_path,
                              exp_id=exp.experiment.public_id, is_ensembl=exp.experiment.is_ensembl)
    if not data_ready:
        exp.log_message("couldn't produce data file", errorfile)
        return None
    exp.experiment.is_ensembl = is_really_ensembl
    if exp.experiment.is_ensembl:
        if organism == 'Homo sapiens':
            exp.experiment.genome = 'GRCh38'
        if organism == 'Mus musculus':
            exp.experiment.genome = 'GRCm38'
        if organism == 'Rattus norvegicus':
            exp.experiment.genome = 'Rnor_6.0'
    else:
        print("not ensembl")
    # get the comparisons
    compare_xml = os.path.join(gea_exp_folder, gea_acc + '-configuration.xml')
    if exists(compare_xml):
        exp.experiment.has_comparison = get_comparisons(compare_xml)
    else:
        exp.log_message("couldn't find comparisons configuration.xml file", errorfile)
        return None
    # check them and remove and with problems (like sample names that aren't in the data)
    exp.experiment.has_comparison = exp.clean_out_invalid_comparisons()
    # write them out to a file
    exp.experiment.comparisons_file_path = os.path.join(rosalind_exp_folder, exp.experiment.public_id + "-comparisons.jlist")
    exp.write_comparisons()
    # write out all experiment metadata to file
    exp.experiment.metadata_file_path = os.path.join(rosalind_exp_folder, exp.experiment.public_id + "-meta.json")
    exp.status = "ready"
    exp.save()
    # if present, capture the GEA analysis results
    results_folder = os.path.join(rosalind_experiments_dir, gea_acc + "/comp_data/")
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    for f in get_files(os.path.join(gea_experiments_dir, gea_acc)):
        f = str(f)
        if f.endswith("-analytics.tsv") and f.startswith(gea_acc):
            shutil.copyfile(os.path.join(gea_experiments_dir, gea_acc, f), results_folder + "all-gene-results.tsv")
        if f.endswith(".gsea.tsv") and f.startswith(gea_acc):
            shutil.copy(os.path.join(gea_experiments_dir, gea_acc, f), results_folder)
    if os.path.exists(os.path.join(results_folder, "all-gene-results.tsv")):
        # split the all-gene-results.tsv into separate files for each comparison
        for comparison in exp.experiment.has_comparison:
            comparison_id = str(comparison.id)
            comparison_file = os.path.join(results_folder, "comparison_"+str(comparison_id) + ".txt")
            # open the all gene results file as a data frame
            # filter the data frame to only include rows with the comparison id
            # write the filtered data frame to the comparison file
            columns_to_keep = ['Gene ID', 'Gene Name', comparison_id+".p-value", comparison_id+".log2foldchange"]
            df = pd.read_csv(os.path.join(results_folder, "all-gene-results.tsv"), sep='\t')
            skip = False
            for c in columns_to_keep:
                if c not in df.columns:
                    exp.log_message("couldn't find column " + c + " in " + results_folder + "all-gene-results.tsv", errorfile)
                    skip = True
            if skip:
                continue
            df = df[columns_to_keep]
            df = df.rename(columns={'Gene ID': 'GeneID', 'Gene Name': 'Symbol', comparison_id+".p-value": 'pvalue',
                                    comparison_id+".log2foldchange": 'log2FoldChange'})
            if exp.experiment.is_ensembl:
                df['EnsemblID'] = df['GeneID']
                # convert the ensembl ids to ncbi ids and put the ncbi ids in the GeneID column
                df['GeneID'] = df['GeneID'].apply(lambda x: ensembl_to_ncbi(x, gene_matcher))
                # remove rows where GeneId is None
                df = df[df['GeneID'].notnull()]
                # remove unuseful data
                df.drop_duplicates(subset=['GeneID'], keep='first', inplace=True)
                df.dropna(how='any', axis=0, inplace=True)
                # filter out genes that are not differentially expressed
                df = df.loc[(df['pvalue'] <= .05) & (df['log2FoldChange'].abs() >= 1)]

            df.to_csv(comparison_file, sep='\t', index=False)
    else:
        exp.log_message("couldn't find all-gene-results.tsv file", errorfile)
    return True


def prepare_all_gea_experiments_for_upload(gea_experiments_dir=None, rosalind_experiments=None, gene_matcher=None):
    """
    Run through all experiments downloaded from the Gene Expression Atlas and prepare files for upload to Rosalind.
    prepare_experiment_for_load should produce 4 files in a new directory under the passed rosalind_experiments dir
    """
    experiments = get_immediate_subdirectories(gea_experiments_dir)
    timestamp = time.strftime('%Y-%m-%d %a %H:%M:%S')
    exps_ready = set()
    exps_problems = set()
    errors = os.path.join(rosalind_experiments, "preprocessing_errors_gea.txt")
    for exp_file_name in experiments:
        print("preparing ", exp_file_name)
        #already processed ?
        full_path_to_exp = os.path.join(rosalind_experiments, exp_file_name, exp_file_name + "-meta.json")
        exp_file_exists = exists(full_path_to_exp)
        if(exp_file_exists):
            exp = read_exp_meta(full_path_to_exp)
            #processing okay?  (if not try again)
            if exp.status == 'ready' or exp.status == 'pre':
                exps_ready.add(exp_file_name)
                #already good to go, skip
                continue
            else:
                exps_problems.add(exp_file_name)
                # take this out to try again.
                continue
        try:
           ready = prepare_experiment_for_load(exp_file_name, gea_experiments_dir, rosalind_experiments,
                                            errorfile=errors, gene_matcher=gene_matcher)
        except Exception as e:
            with open(errors, "a") as error_log:
                error_log.write(exp_file_name + "\t" + str(e) + "\n")
            continue
        with open(os.path.join(rosalind_experiments, "load_status_gea.txt"), "a") as status_log:
            if ready:
                exps_ready.add(exp_file_name)
                status_log.write(exp_file_name+"\tready\t"+timestamp+"\n")
            else:
                exps_problems.add(exp_file_name)
                status_log.write(exp_file_name + "\tfailed preprocessing\t" + timestamp + "\n")
    print("experiment directories processed:", str(len(experiments)))
    print("ready: ", str(len(exps_ready)))
    print("problems: ", str(len(exps_problems)))


if __name__ == '__main__':
    print("hello GEA")
    gea_data_for_import = gea_dir
    gea_acc = "E-GEOD-15389" #"'E-GEOD-48459' #'E-GEOD-21293' #
    rosalind_experiments = os.path.join(ROOT_DIR, 'data', 'experiments', 'transformed', 'gea')
 #   prepare_experiment_for_load(gea_acc, gea_data_for_import, rosalind_experiments)
    gene_matcher = TextMatcher()
    keep_categories = ['Gene']
    node_file = '/Users/bgood/Documents/GitHub/franklin/etl/data/merged/merged-kg_nodes.tsv'
    gene_matcher.make_string_to_thing_maps_from_franklin_nodes(node_file, keep_categories=keep_categories)
    prepare_all_gea_experiments_for_upload(gea_experiments_dir=gea_data_for_import,
                                           rosalind_experiments=rosalind_experiments,
                                           gene_matcher=gene_matcher)
    # experiment directories processed: 3010
    # ready:  2635
    # problems:  375
    print("goodbye GEA")

#     experiment directories processed: 4465
# ready:  3603
# problems:  861
# goodbye GEA
