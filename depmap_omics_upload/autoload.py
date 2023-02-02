import json
import os
import numpy as np
import pandas as pd
import dalmatian as dm
import datetime
import logging
import sys

import tracker as track
from genepy import sequencing as seq
from genepy.utils import helper as h
from genepy.google import gcp

#####################
# Loading Functions
#####################

def loadFromMultipleWorkspaces(
    gumbo_env,
    wsnames,
    wsidcol,
    gumboidcol,
    stype,
    extract,
    maxage,
    minsizes_bam,
    minsizes_cram,
    only_load_mapped=False,
    bamcol="cram_or_bam_path",
    load_undefined=False,
    accept_unknowntypes=True,
    addonly=[],
):
    """
    Load new RNAseq samples from multiple terra workspace and attempt to map them to existing profiles in gumbo.

    Args:
    -----
        wsnames (list): list of tuples mapping source to workspace name ([(source, ws name)])
        wsidcol (str): name of column in wsname's sample table that contains IDs to map to ProfileIDs by
        gumboidcol (str): name of column in gumbo's profile table that contains IDs to map to ProfileIDs by
        stype (str): type of the data (wgs, rna, etc.)
        extract: if you want to specify what values should refer to which column names

    Returns:
    --------
        samples: pd dataframe the filtered sample list
        unmapped_new_lines: list of lines that don't have PR ids assigned in gumbo
    """
    samples = []
    for s, wsname, ftype in wsnames:
        logging.info("loading " + stype + " samples from terra workspace: " + wsname)
        bamcol = "cram_or_bam_path" if ftype == "bam" else "cram_path"
        samples_per_ws = loadFromTerraWorkspace(
            gumbo_env,
            wsname,
            wsidcol,
            gumboidcol,
            s,
            stype,
            ftype,
            bamcol,
            only_load_mapped,
            load_undefined,
            extract,
            accept_unknowntypes=accept_unknowntypes,
            addonly=addonly,
            maxage=maxage,
            minsizes_bam=minsizes_bam,
            minsizes_cram=minsizes_cram
        )
        samples.append(samples_per_ws)
    return pd.concat(samples)


def loadFromTerraWorkspace(
    gumbo_env,
    wsname,
    wsidcol,
    gumboidcol,
    source,
    stype,
    ftype,
    bamcol,
    only_load_mapped=False,
    load_undefined=False,
    extract={},
    accept_unknowntypes=True,
    addonly=[],
    maxage="",
    minsizes_bam={},
    minsizes_cram={},
):
    """
    Load new samples from a terra workspace and attempt to map them to existing profiles in gumbo.

    Args:
    -----
        wsname (str): name of terra workspace to load data from
        wsidcol (str): name of column in wsname's sample table that contains IDs to map to ProfileIDs by
        gumboidcol (str): name of column in gumbo's profile table that contains IDs to map to ProfileIDs by
        source (str): source of the data delivered (DEPMAP, IBM, etc.)
        stype (str): type of the data (wgs, rna, etc.)
        ftype (str): type of the sequenced file (bam/cram)
        extract: if you want to specify what values should refer to which column names
        dict{
        'name':

    Returns:
    --------
        samples: pd dataframe the filtered sample list
        unmapped_new_lines: list of lines that don't have PR ids assigned in gumbo
    """
    mytracker = track.SampleTracker(gumbo_env=gumbo_env)
    seq_table = mytracker.read_seq_table()
    pr_table = mytracker.read_pr_table()
    wm = dm.WorkspaceManager(wsname).disable_hound()
    samples_in_ws = wm.get_samples().replace(np.nan, "", regex=True).reset_index()

    # check broken bam file paths
    logging.info("checking if there are any broken file paths")
    foundfiles = gcp.lsFiles(samples_in_ws[bamcol])
    broken_bams = set(samples_in_ws[bamcol]) - set(foundfiles)
    logging.warning(
        "These "
        + str(len(broken_bams))
        + " bam file path do not exist: "
        + str(broken_bams)
    )
    samples = samples_in_ws[(~samples_in_ws[bamcol].isin(broken_bams))]

    logging.info(
        "extracting information from workspace including hash, size, update time, etc."
    )
    samples = extractFromWorkspace(samples, stype, ftype, bamcol, extract=extract, minsizes_bam=minsizes_bam, minsizes_cram=minsizes_cram)
    logging.info("generating CDS-ids, annotating source, and renaming columns")
    samples = mapSamples(samples, source, ftype, extract=extract)
    logging.info("checking for duplicates in the workspace by comparing file sizes")
    samples = resolveFromWorkspace(
        samples,
        seq_table[seq_table[extract["expected_type"]] == stype],
        wsidcol,
        ftype,
        accept_unknowntypes,
        addonly,
        extract,
    )
    samples = samples[samples[extract["update_time"]] > maxage]
    if len(samples) == 0:
        logging.warning("no new data available")
        return None

    # init profile id column
    samples[extract["profile_id"]] = ""
    samples[extract["version"]] = 1
    samples[extract["expected_type"]] = stype
    # FOR NOW, assume the id col to map by is in profile table (might change??)
    for k, v in samples.iterrows():
        if not pd.isnull(v[wsidcol]):
            # different datatypes from the same line might share the same SM-ID,
            # so mapping should condition on datatype as well;
            # sometimes there will be multiple SM-ids associated with one bam on Terra
            # in the format of "SM-xxxxx,SM-yyyyy", so split them
            pr_id = pr_table[
                (~pr_table[gumboidcol].isnull()) & (pr_table[gumboidcol].isin(set(v[wsidcol].split(",")))) & (pr_table.Datatype == stype) & (pr_table.BlacklistOmics != True)
            ].index.tolist()
            if len(pr_id) > 1:
                raise ValueError(
                    "multiple profile ids mapped to the same validation id. check with ops!"
                )
            elif len(pr_id) == 1:
                samples.loc[k, extract["profile_id"]] = pr_id[0]

    return samples


def extractFromWorkspace(
    samples,
    stype,
    ftype,
    bamcol,
    recomputeTime=True,
    recomputesize=True,
    recomputedate=True,
    recompute_hash=True,
    extract={},
    minsizes_bam={},
    minsizes_cram={}
):
    """
    Extract more information from a list of samples found on GP workspaces

    Args:
    -----
        samples: pd dataframes of samples with at least arxspan ids and sizes
        stype: str sequencing type
        recomputeTime: bool whether to recompute the date of upload of the bam file
        recomputesize: bool whether to recompute the of the bam file
        recomputehash: bool whether to recompute the of the bam file
        extract: if you want to specify what values should refer to which column names
        dict{
        'name':
        'bai':
        'bam':
        'source':
        'from_arxspan_id':
        ...} (see extract_defaults)

    Returns:
    --------
        samples: pd dataframe the filtered sample list
    """
    if extract["legacy_hash"] not in samples.columns or recompute_hash:
        samples[extract["hash"]] = [
            gcp.extractHash(val)
            for val in gcp.lsFiles(samples[bamcol].tolist(), "-L", 200)
        ]
    lis = gcp.lsFiles(samples[bamcol].tolist(), "-al", 200)
    if extract["legacy_size"] not in samples.columns or recomputesize:
        samples[extract["legacy_size"]] = [gcp.extractSize(i)[1] for i in lis]
    if extract["update_time"] not in samples.columns or recomputeTime:
        samples[extract["update_time"]] = [gcp.extractTime(i) for i in lis]
    todrop = []
    for k, val in samples.iterrows():
        if ftype == "bam":
            if val[extract["legacy_size"]] < minsizes_bam[stype]:
                todrop.append(k)
                logging.warning(
                    "too small size, removing sample: "
                    + str(val[extract["root_sample_id"]])
                )
        elif ftype == "cram":
            if val[extract["legacy_size"]] < minsizes_cram[stype]:
                todrop.append(k)
                logging.warning(
                    "too small size, removing sample: "
                    + str(val[extract["root_sample_id"]])
                )
    samples = samples.drop(index=todrop)
    # getting the date released
    if len(samples) == 0:
        return None
    if extract["release_date"] not in samples.columns or recomputedate:
        samples[extract["release_date"]] = seq.getBamDate(samples[bamcol])
    return samples


def mapSamples(samples, source, ftype, extract={}):
    """
    Convert samples from a list of GP workspaces to something being able to be merged with the sample tracker

    Args:
    -----
        samples: pd dataframes of samples with at least arxspan ids and sizes
        extract: if you want to specify what values should refer to which column names
        dict{
        'name':
        'bai':
        'bam':
        'source':
        'from_arxspan_id':
        ...} (see extract_defaults)
        source:

    Returns:
    --------
        samples: pd dataframe the filtered sample list
    """
    # creating unique ids
    samples[extract["ref_id"]] = [
        "CDS-" + h.randomString(stringLength=6, stype="all", withdigits=True)
        for _ in range(len(samples))
    ]
    samples.reset_index(drop=True, inplace=True)
    samples[extract["source"]] = source

    # renamings
    samples = samples.rename(
        columns={
            extract["bam"]: extract["ref_bam"],
            extract["bai"]: extract["ref_bai"],
            extract["root_sample_id"]: extract["sm_id"],
            # extract["PDO_id_terra"]: extract["PDO_id_gumbo"],
            extract["cram"]: extract["ref_cram"],
            extract["crai"]: extract["ref_crai"],
        }
    ).set_index(extract["ref_id"], drop=True)
    # subsetting
    if ftype == "bam":
        samples = samples[
            list(
                set(
                    [
                        extract["ref_bam"],
                        extract["ref_bai"],
                        extract["release_date"],
                        extract["legacy_size"],
                        # extract["PDO_id_gumbo"],
                        extract["sm_id"],
                        extract["update_time"],
                        extract["source"],
                    ]
                )
            )
        ]
    else:
        samples = samples[
            list(
                set(
                    [
                        extract["ref_cram"],
                        extract["ref_crai"],
                        extract["release_date"],
                        extract["legacy_size"],
                        # extract["PDO_id_gumbo"],
                        extract["sm_id"],
                        extract["update_time"],
                        extract["source"],
                    ]
                )
            )
        ]
    return samples


def resolveFromWorkspace(
    samples, refsamples, wsidcol, ftype, accept_unknowntypes=True, addonly=[], extract={},
):
    """
    Filters our list by trying to find duplicate in our dataset and remove any sample that isn't tumor

    Args:
    -----
        match: list[str]|str the possible values that a sample id need to contain to be considered valid
        participantslicepos: int the length of the sample id string
        accept_unknowntypes: bool whether or not the sample type column for that sample can be different from "Tumor"
        refsamples: pd dataframe representing a sample tracker
        samples: pd dataframes of samples with at least arxspan ids and sizes
        extract: if you want to specify what values should refer to which column names
        dict{
        'name':
        'bai':
        'bam':
        'source':
        'from_arxspan_id':
        ...} (see extract_defaults)

    Returns:
    --------
        samples: pd dataframe the filtered sample list
    """
    if ftype == "bam":
        sample_size = {
            gcp.extractSize(val)[1]: gcp.extractSize(val)[0]
            for val in gcp.lsFiles(samples[extract["ref_bam"]].tolist(), "-la")
        }
    else:
        sample_size = {
            gcp.extractSize(val)[1]: gcp.extractSize(val)[0]
            for val in gcp.lsFiles(samples[extract["ref_cram"]].tolist(), "-la")
        }
    dups_to_remove = [
        sample_size[a]
        for a in set(sample_size.keys()) & set(refsamples[extract["legacy_size"]])
    ]
    # remove the duplicates from consideration
    logging.info("Len of samples before removal: " + str(len(samples)))
    logging.info(
        "Dups from this workspace has len "
        + str(len(dups_to_remove))
        + ":\n "
        + str(dups_to_remove)
    )
    # remove the samples with broken bam filepaths from consideration
    if ftype == "bam":
        samples = samples[~samples[extract["ref_bam"]].isin(dups_to_remove)]
    else:
        samples = samples[~samples[extract["ref_cram"]].isin(dups_to_remove)]

    logging.info("Len of samples after removal: " + str(len(samples)))
    if len(samples) == 0:
        return None

    # if only add some samples
    if len(addonly) > 0:
        samples = samples[samples[wsidcol].isin(addonly)]

    # unknown types
    if "sample_type" in samples.columns:
        if not accept_unknowntypes:
            samples = samples[samples["sample_type"].isin(["Tumor"])]
    return samples

if __name__ == '__main__':
    config_filename = sys.argv[1]
    with open(config_filename, "r") as f:
        config = json.load(f)
    
    today = str(datetime.datetime.now().replace(microsecond=0)).replace(" ", "_").replace(":", "-")

    if not os.path.exists(config["savefolder"]):
        os.makedirs(config["savefolder"])

    logging.basicConfig(filename=config["savefolder"] + today + "_" + "-".join(config["datatypes"]) + '.log', level=logging.INFO)

    if "rna" in config["datatypes"]:
        rnasamples = loadFromMultipleWorkspaces(config['gumbo_env'], config['rnaworkspaces'], config["extract_defaults"]["sm_id"], "SMIDOrdered", "rna", config["extract_defaults"], config['maxage'], config['minsizes_bam'], config['minsizes_cram'], bamcol="cram_or_bam_path")
        rnasamples[rnasamples[config["extract_defaults"]["profile_id"]] != ""].to_csv(config["savefolder"] + today + "_" + "mappedRNAsamples.csv")
        rnasamples[rnasamples[config["extract_defaults"]["profile_id"]] == ""].to_csv(config["savefolder"] + today + "_" + "unmappedRNAsamples.csv")

    if "wgs" in config["datatypes"]:
        wgssamples = loadFromMultipleWorkspaces(config['gumbo_env'], config['wgsworkspaces'], config["extract_defaults"]["sm_id"], "SMIDOrdered", "wgs", config["extract_defaults"], config['maxage'], config['minsizes_bam'], config['minsizes_cram'], bamcol="cram_path")
        wgssamples[wgssamples[config["extract_defaults"]["profile_id"]] != ""].to_csv(config["savefolder"] + today + "_" + "mappedWGSsamples.csv")
        wgssamples[wgssamples[config["extract_defaults"]["profile_id"]] == ""].to_csv(config["savefolder"] + today + "_" + "unmappedWGSsamples.csv")