from __future__ import print_function
import pandas as pd
from datetime import date

from genepy.utils import helper as h
from depmap_omics_upload import tracker as track
from depmapomics.config import *
from taigapy import TaigaClient


def getPRToRelease(today=None, portals=DATASETS, date_col_dict=DATE_COL_DICT):
    """generate lists of profiles to release based on date for all portals

    Returns:
        prs (dict{(portal: list of PRs)}): for each portal, list of profile IDs
    """
    if today == None:
        today = date.today()
    mytracker = track.SampleTracker()
    pr_table = mytracker.read_pr_table()
    mytracker.close_gumbo_client()
    pr_table = pr_table[pr_table.BlacklistOmics != True]
    prs = dict()
    for k, v in date_col_dict.items():
        prs_with_date = pr_table[~(pr_table[v].isnull())]
        if k in portals:
            prs[k] = prs_with_date[
                (prs_with_date[v] <= today)
                & (prs_with_date.ProfileSource != "taiga")
                & (~prs_with_date.MainSequencingID.isnull())
            ].index.tolist()
    assert (len(set(prs['dmc']) - set(prs['internal'])) == 0), "Lines with DMC release dates missing internal release dates: " + str(set(prs['dmc']) - set(prs['internal']))
    return prs


def makeAchillesChoiceTable(
    prs, source_priority=SOURCE_PRIORITY, colnames=ACH_CHOICE_TABLE_COLS,
):
    """generate a table for each portal that indicates which profiles are released corresponding to which MC

    Args:
        prs (list): list of profile IDs to be released
        source_priority (list, optional): ordered list of how different data sources should be prioritized

    Returns:
        ach_table (pd.DataFrame): a df containing MC-PR mapping
    """
    mytracker = track.SampleTracker()
    pr_table = mytracker.read_pr_table()
    seq_table = mytracker.read_seq_table()
    source_priority = {source_priority[i]: i for i in range(len(source_priority))}
    rows = []
    subset_pr_table = pr_table.loc[prs]
    subset_pr_table = subset_pr_table[subset_pr_table.BlacklistOmics != 1]
    mcs = set(subset_pr_table["ModelCondition"])
    mytracker.close_gumbo_client()
    # we're only picking one PR per datatype (rna/dna) for each MC
    for mc in mcs:
        prs_in_mc = subset_pr_table[(subset_pr_table.ModelCondition == mc)]
        # rna
        # if there is only one rna PR associated with this MC, pick that PR
        if len(prs_in_mc[prs_in_mc.Datatype == "rna"]) == 1:
            pr = prs_in_mc[prs_in_mc.Datatype == "rna"].index[0]
            rows.append((mc, pr, "rna"))
        # if there are more than one PRs associated with this MC,
        # pick the PR from the most prioritized source according to the
        # ranking in source_priority
        elif len(prs_in_mc[prs_in_mc.Datatype == "rna"]) > 1:
            cds_ids = prs_in_mc[prs_in_mc.Datatype == "rna"].MainSequencingID.tolist()
            # at this point it is guaranteed that all cds_ids have different sources
            subset_seq_table = seq_table[seq_table.index.isin(cds_ids)]
            subset_seq_table.source = subset_seq_table.source.replace(source_priority)
            latest_cds_id = subset_seq_table.loc[cds_ids, "source"].idxmin()
            pr = subset_pr_table[subset_pr_table.MainSequencingID == latest_cds_id].index[0]
            rows.append((mc, pr, "rna"))
        # dna
        # if there is only one dna (wes + wgs) PR associated with this MC, pick that PR
        if (
            len(
                prs_in_mc[(prs_in_mc.Datatype == "wgs") | (prs_in_mc.Datatype == "wes")]
            )
            == 1
        ):
            pr = prs_in_mc[
                (prs_in_mc.Datatype == "wgs") | (prs_in_mc.Datatype == "wes")
            ].index[0]
            rows.append((mc, pr, "dna"))
        # if there are more than one PRs associated with this MC
        elif (
            len(
                prs_in_mc[(prs_in_mc.Datatype == "wgs") | (prs_in_mc.Datatype == "wes")]
            )
            > 1
        ):
            cds_ids_wgs = prs_in_mc[prs_in_mc.Datatype == "wgs"].MainSequencingID.tolist()
            cds_ids_wes = prs_in_mc[prs_in_mc.Datatype == "wes"].MainSequencingID.tolist()
            pr = ""
            # if this MC doesn't have any wgs PRs
            if len(cds_ids_wgs) == 0:
                # out of the wes PRs, select the PR from the most prioritized source
                # according to the ranking in source_priority
                subset_seq_table = seq_table[seq_table.index.isin(cds_ids_wes)]
                subset_seq_table.source = subset_seq_table.source.replace(
                    source_priority
                )
                latest_cds_id_wes = subset_seq_table.loc[cds_ids_wes, "source"].idxmin()
                pr = subset_pr_table[subset_pr_table.MainSequencingID == latest_cds_id_wes].index[
                    0
                ]
            # if this MC has wgs PR(s)
            else:
                # ignore the wes PRs, select the PR from the most prioritized source
                # according to the ranking in source_priority from the wgs PR(s)
                subset_seq_table = seq_table[seq_table.index.isin(cds_ids_wgs)]
                subset_seq_table.source = subset_seq_table.source.replace(
                    source_priority
                )
                latest_cds_id_wgs = subset_seq_table.loc[cds_ids_wgs, "source"].idxmin()
                pr = subset_pr_table[subset_pr_table.MainSequencingID == latest_cds_id_wgs].index[
                    0
                ]
            rows.append((mc, pr, "dna"))
    ach_table = pd.DataFrame(rows, columns=colnames)

    return ach_table


def makeDefaultModelTable(
    prs, source_priority=SOURCE_PRIORITY, colnames=DEFAULT_TABLE_COLS,
):
    """generate a table that indicates which profiles are released corresponding to which modelID

    Args:
        prs (list): list of profile IDs to be released
        source_priority (list, optional): ordered list of how different data sources should be prioritized

    Returns:
        default_table (pd.DataFrame): a df containing Model-PR mapping
    """
    mytracker = track.SampleTracker()
    pr_table = mytracker.read_pr_table()
    pr_table = pr_table[pr_table.BlacklistOmics != True]
    mc_table = mytracker.read_mc_table()
    seq_table = mytracker.read_seq_table()
    source_priority = {source_priority[i]: i for i in range(len(source_priority))}
    rows = []
    subset_pr_table = pr_table.loc[prs]
    subset_pr_table = subset_pr_table[subset_pr_table.BlacklistOmics != 1]
    mcs = set(subset_pr_table["ModelCondition"])
    models = set(mc_table.loc[mcs].ModelID)
    mytracker.close_gumbo_client()
    # we're only picking one PR per datatype (rna/dna) for each Model
    for m in models:
        subset_mc_table = mc_table[mc_table.ModelID == m]
        mcs_in_model = subset_mc_table.index.tolist()
        prs_in_model = subset_pr_table[
            (subset_pr_table.ModelCondition.isin(mcs_in_model))
        ]
        # rna
        # if there is only one rna PR associated with this Model, pick that PR
        if len(prs_in_model[prs_in_model.Datatype == "rna"]) == 1:
            pr = prs_in_model[prs_in_model.Datatype == "rna"].index[0]
            rows.append((m, pr, "rna"))
        # if there are more than one PRs associated with this Model,
        # pick the PR from the most prioritized source according to the
        # ranking in source_priority
        elif len(prs_in_model[prs_in_model.Datatype == "rna"]) > 1:
            cds_ids = prs_in_model[prs_in_model.Datatype == "rna"].MainSequencingID.tolist()
            subset_seq_table = seq_table[seq_table.index.isin(cds_ids)]
            subset_seq_table.source = subset_seq_table.source.replace(source_priority)
            latest_cds_id = subset_seq_table.loc[cds_ids, "source"].idxmin()
            pr = subset_pr_table[subset_pr_table.MainSequencingID == latest_cds_id].index[0]
            rows.append((m, pr, "rna"))
        # dna
        # if there is only one dna (wes + wgs) PR associated with this Model, pick that PR
        if (
            len(
                prs_in_model[
                    (prs_in_model.Datatype == "wgs") | (prs_in_model.Datatype == "wes")
                ]
            )
            == 1
        ):
            pr = prs_in_model[
                (prs_in_model.Datatype == "wgs") | (prs_in_model.Datatype == "wes")
            ].index[0]
            rows.append((m, pr, "dna"))
        # if there are more than one PRs associated with this Model
        elif (
            len(
                prs_in_model[
                    (prs_in_model.Datatype == "wgs") | (prs_in_model.Datatype == "wes")
                ]
            )
            > 1
        ):
            cds_ids_wgs = prs_in_model[prs_in_model.Datatype == "wgs"].MainSequencingID.tolist()
            cds_ids_wes = prs_in_model[
                (prs_in_model.Datatype == "wes") & (prs_in_model.MainSequencingID != "")
            ].MainSequencingID.tolist()  # MainSequencingID is '' when the profile is in legacy
            pr = ""
            # if no wgs, look at MC table and select the most prioritized source
            # according to the ranking in source_priority
            if len(cds_ids_wgs) == 0:
                subset_seq_table = seq_table[seq_table.index.isin(cds_ids_wes)]
                subset_seq_table.source = subset_seq_table.source.replace(
                    source_priority
                )
                latest_cds_id_wes = subset_seq_table.loc[cds_ids_wes, "source"].idxmin()
                pr = subset_pr_table[subset_pr_table.MainSequencingID == latest_cds_id_wes].index[
                    0
                ]
            # if there is wgs, always select wgs
            else:
                subset_seq_table = seq_table[seq_table.index.isin(cds_ids_wgs)]
                subset_seq_table.source = subset_seq_table.source.replace(
                    source_priority
                )
                latest_cds_id_wgs = subset_seq_table.loc[cds_ids_wgs, "source"].idxmin()
                pr = subset_pr_table[subset_pr_table.MainSequencingID == latest_cds_id_wgs].index[
                    0
                ]
            rows.append((m, pr, "dna"))
    default_table = pd.DataFrame(rows, columns=colnames)
    return default_table

def makeProfileTable(prs, columns=PROFILE_TABLE_COLS):
    """subset gumbo profile table both column- and row-wise for the release

    Args:
        prs (list): list of profile IDs to be released
        columns (list, optional): list of columns to be included in the table

    Returns:
        pr_table (pd.DataFrame): a df containing omics profile information
    """
    mytracker = track.SampleTracker()
    pr_table = mytracker.add_model_cols_to_prtable(["ModelID"])
    pr_table = pr_table.loc[prs, columns]
    pr_table = pr_table[pr_table["Datatype"].isin(["rna", "wgs", "wes", "SNParray"])]
    mytracker.close_gumbo_client()
    return pr_table

def initVirtualDatasets(
    samplesetname, taiga_folder_id=VIRTUAL_FOLDER, portals=DATASETS
):
    """initialize both PR- and Model-level taiga virtual datasets for all 4 portals by uploading an empty dummy file
    """
    virutal = dict()
    tc = TaigaClient()
    for p in portals:
        virutal[p] = tc.create_dataset(
            p + "_" + samplesetname,
            dataset_description=samplesetname
            + " release of the DepMap dataset for the DepMap Portal. Please look at the README file for additional information about this dataset. ",
            upload_files=[
                {
                    "path": "/dev/null",
                    "name": "init",
                    "format": "Raw",
                    "encoding": "utf-8",
                }
            ],
            folder_id=taiga_folder_id,
        )
    return virutal


def uploadPRMatrix(
    prs,
    taiga_latest,
    taiga_virtual,
    latest_fn,
    virtual_fn,
    matrix_format,
    pr_col="index",
    folder=WORKING_DIR,
    change_desc="",
):
    """subset, save and upload to taiga PR-level matrix

    Args:
        prs (list): list of PR-ids to release
        taiga_latest (str): which dataset the matrices to be subsetted are being read from
        taiga_virtual (str): which dataset the matrices are being uploaded to
        latest_fn (str): file name on taiga latest
        virtual_fn (str): file name on taiga virtual
        matrix_format (str): which format this matrix should be uploaded in (NumericMatrixCSV, TableCSV, etc)
        pr_col (str): which column in the matrix contains pr ids
        folder (str): where the file should be stores before uploading to virtual
        change_desc (str): change description on taiga virtual
    """
    print("loading ", latest_fn, " from latest")
    tc = TaigaClient()
    to_subset = tc.get(name=taiga_latest, file=latest_fn)

    print("subsetting ", latest_fn)
    if pr_col == "index":
        subset_mat = to_subset[to_subset.index.isin(prs)]
        subset_mat.to_csv(folder + virtual_fn + ".csv")
    else:
        subset_mat = to_subset[to_subset[pr_col].isin(prs)]
        subset_mat = subset_mat.rename(columns={pr_col: "ProfileID"})
        subset_mat.to_csv(folder + virtual_fn + ".csv", index=False)

    print("uploading ", virtual_fn, " to virtual")
    tc.update_dataset(
        dataset_id=taiga_virtual,
        changes_description=change_desc,
        upload_files=[
            {
                "path": folder + virtual_fn + ".csv",
                "name": virtual_fn,
                "format": matrix_format,
                "encoding": "utf-8",
            },
        ],
        add_all_existing_files=True,
    )


def uploadModelMatrix(
    pr2model_dict,
    taiga_latest,
    taiga_virtual,
    latest_fn,
    virtual_fn,
    matrix_format,
    pr_col="index",
    folder=WORKING_DIR,
    change_desc="",
):
    """subset, rename, save and upload to taiga model-level matrix

    Args:
        pr2model_dict (dict): dictionary mapping profile ids to model ids
        taiga_latest (str): which dataset the matrices to be subsetted are being read from
        taiga_virtual (str): which dataset the matrices are being uploaded to
        latest_fn (str): file name on taiga latest
        virtual_fn (str): file name on taiga virtual
        matrix_format (str): which format this matrix should be uploaded in (NumericMatrixCSV, TableCSV, etc)
        pr_col (str): which column in the matrix contains pr ids
        folder (str): where the file should be stores before uploading to virtual
        change_desc (str): change description on taiga virtual
    """
    print("loading ", latest_fn, " from latest")
    tc = TaigaClient()
    to_subset = tc.get(name=taiga_latest, file=latest_fn)

    print("subsetting ", latest_fn)
    if pr_col == "index":
        subset_mat = to_subset[to_subset.index.isin(set(pr2model_dict.keys()))].rename(
            index=pr2model_dict
        )
        subset_mat.to_csv(folder + virtual_fn + ".csv")
    else:
        subset_mat = to_subset[
            to_subset[pr_col].isin(set(pr2model_dict.keys()))
        ].replace({SAMPLEID: pr2model_dict})
        subset_mat.to_csv(folder + virtual_fn + ".csv", index=False)

    print("uploading ", virtual_fn, " to virtual")
    tc.update_dataset(
        dataset_id=taiga_virtual,
        changes_description=change_desc,
        upload_files=[
            {
                "path": folder + virtual_fn + ".csv",
                "name": virtual_fn,
                "format": matrix_format,
                "encoding": "utf-8",
            },
        ],
        add_all_existing_files=True,
    )


def uploadBinaryGuideMutationMatrixModel(
    pr2model_dict, portal, taiga_latest=TAIGA_MUTATION, fn_mapping=VIRTUAL_FILENAMES_GUIDEMUT, taiga_virtual="",
):
    """subset, rename, save and upload to taiga germline binary matrix
    
    Args:
        pr2model_dict (dict): renaming scheme mapping from PR-id to model id
        portal (str): which portal the data is being uploaded to
        taiga_latest (str): which dataset the matrices to be subsetted are being read from
        taiga_virtual (str): which dataset the matrices are being uploaded to
    """
    folder = WORKING_DIR + portal + "/model/"
    h.createFoldersFor(folder)
    for latest_fn, virtual_fn in fn_mapping.items():
        # load pr-id indexed matrices for the current quarter
        print("Guide mutation matrix: loading from taiga latest")
        tc = TaigaClient()
        germline = tc.get(name=taiga_latest, file=latest_fn)

        # subset and rename
        print("Guide mutation matrix: subsetting and renaming")
        whitelist = [x for x in germline.columns if x in pr2model_dict]
        whitelist_germline = germline[whitelist]
        whitelist_germline = whitelist_germline.rename(columns=pr2model_dict)
        whitelist_germline = whitelist_germline.astype(bool).astype(int)
        sorted_mat = germline.iloc[:, :4].join(whitelist_germline)
        sorted_mat["end"] = sorted_mat["end"].astype(int)
        # ONE-OFF:
        sorted_mat = sorted_mat.rename(columns={"foldchange": "sgRNA"})
        sorted_mat.to_csv(folder + virtual_fn + ".csv", index=False)

        # upload to taiga
        print("Guide mutation: uploading to taiga")
        tc.update_dataset(
            dataset_id=taiga_virtual,
            changes_description="adding model-level germline matrix",
            upload_files=[
                {
                    "path": folder + virtual_fn + ".csv",
                    "name": virtual_fn,
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
            ],
            add_all_existing_files=True,
        )


def uploadAuxTables(
    taiga_ids=VIRTUAL,
    ach_table_name=ACH_CHOICE_TABLE_NAME,
    default_table_name=DEFAULT_TABLE_NAME,
    release_pr_table_name=RELEASE_PR_TABLE_NAME,
    folder=WORKING_DIR + SAMPLESETNAME + "/",
    today=None,
):
    """upload achilles choice and default model table to all portals
    Args:
    
        taiga_ids (dict, optional): dict mapping portal name to taiga virtual dataset id
        folder (str, optional): where the tables are saved
    """
    prs_allportals = getPRToRelease(portals=taiga_ids.keys(), today=today)
    for portal, prs in prs_allportals.items():
        achilles_table = makeAchillesChoiceTable(prs)
        default_table = makeDefaultModelTable(prs)
        profile_table = makeProfileTable(prs)
        achilles_table.to_csv(
            folder + portal + "_" + ach_table_name + ".csv", index=False
        )
        default_table.to_csv(
            folder + portal + "_" + default_table_name + ".csv", index=False
        )
        profile_table.to_csv(
            folder + portal + "_" + release_pr_table_name + ".csv"
        )
        tc = TaigaClient()
        tc.update_dataset(
            dataset_id=taiga_ids[portal],
            changes_description="adding mapping tables",
            upload_files=[
                {
                    "path": folder + "/" + portal + "_" + ach_table_name + ".csv",
                    "name": "OmicsDefaultModelConditionProfiles",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "/" + portal + "_" + default_table_name + ".csv",
                    "name": "OmicsDefaultModelProfiles",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "/" + portal + "_" + release_pr_table_name + ".csv",
                    "name": "OmicsProfiles",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
            ],
            add_all_existing_files=True,
        )


def makePRLvMatrices(virtual_ids=VIRTUAL, files_nummat=LATEST2FN_NUMMAT_PR, folder=WORKING_DIR + SAMPLESETNAME, files_table=LATEST2FN_TABLE_PR, today=None):
    """for each portal, save and upload profile-indexed data matrices
    
    Args:
        taiga_ids (dict): dictionary that maps portal name to virtual taiga dataset id

    Returns:
        prs (dict{(portal: list of PRs)}): for each portal, list of profile IDs
    """
    prs_allportals = getPRToRelease(portals=virtual_ids.keys(), today=today)
    for portal, prs_to_release in prs_allportals.items():
        print("uploading profile-level matrices to ", portal)
        for latest_id, fn_dict in files_nummat.items():
            for latest, virtual in fn_dict.items():
                uploadPRMatrix(
                    prs_to_release,
                    latest_id,
                    virtual_ids[portal],
                    latest,
                    virtual,
                    "NumericMatrixCSV",
                    pr_col="index",
                    folder=folder+ "/",
                    change_desc="adding " + virtual,
                )
        for latest_id, fn_dict in files_table.items():
            for latest, virtual in fn_dict.items():
                uploadPRMatrix(
                    prs_to_release,
                    latest_id,
                    virtual_ids[portal],
                    latest,
                    virtual,
                    "TableCSV",
                    pr_col=SAMPLEID,
                    folder=folder+ "/",
                    change_desc="adding " + virtual,
                )


def makeModelLvMatrices(virtual_ids=VIRTUAL, folder=WORKING_DIR + SAMPLESETNAME, files_nummat=LATEST2FN_NUMMAT_MODEL, files_table=LATEST2FN_TABLE_MODEL, upload_guide_matrices=True, today=None):
    """for each portal, save and upload profile-indexed data matrices
    
    Args:
        taiga_ids (dict): dictionary that maps portal name to virtual taiga dataset id

    Returns:
        prs (dict{(portal: list of PRs)}): for each portal, list of profile IDs
    """
    prs_allportals = getPRToRelease(portals=virtual_ids.keys(), today=today)
    for portal, prs_to_release in prs_allportals.items():
        default_table = makeDefaultModelTable(prs_to_release)
        pr2model_dict = dict(list(zip(default_table.ProfileID, default_table.ModelID)))
        h.dictToFile(pr2model_dict, folder + "/" + portal + "_pr2model_renaming.json")
        print("uploading model-level matrices to", portal)
        for latest_id, fn_dict in files_nummat.items():
            for latest, virtual in fn_dict.items():
                uploadModelMatrix(
                    pr2model_dict,
                    latest_id,
                    virtual_ids[portal],
                    latest,
                    virtual,
                    "NumericMatrixCSV",
                    pr_col="index",
                    folder=folder+ "/",
                    change_desc="adding " + virtual,
                )
        for latest_id, fn_dict in files_table.items():
            for latest, virtual in fn_dict.items():
                uploadModelMatrix(
                    pr2model_dict,
                    latest_id,
                    virtual_ids[portal],
                    latest,
                    virtual,
                    "TableCSV",
                    pr_col=SAMPLEID,
                    folder=folder+ "/",
                    change_desc="adding " + virtual,
                )
        if upload_guide_matrices:
            uploadBinaryGuideMutationMatrixModel(
                pr2model_dict, portal, taiga_virtual=virtual_ids[portal]
            )


def findLatestVersion(dataset, approved_only=True):
    highest = 0
    latest_version = 0
    tc = TaigaClient()
    data = tc.get_dataset_metadata(dataset)
    for val in data["versions"]:
        if val["state"] == "approved" or not approved_only:
            if int(val["name"]) > highest:
                highest = int(val["name"])
                latest_version = highest
    if latest_version == 0:
        raise ValueError("could not find a version")
    return data["permanames"][0] + "." + str(latest_version)


def updateEternal(
    eternal_id=TAIGA_ETERNAL_UPLOAD, virtual=VIRTUAL, samplesetname=SAMPLESETNAME
):
    """update taiga eternal dataset by linking to latest virtual internal dataset"""
    latest_version = findLatestVersion(virtual["internal"])

    files = [
        VIRTUAL_FILENAMES_NUMMAT_EXP.values()
        + VIRTUAL_FILENAMES_NUMMAT_CN.values()
        + VIRTUAL_FILENAMES_NUMMAT_MUT.values()
        + VIRTUAL_FILENAMES_GERMLINE.values()
        + VIRTUAL_FILENAMES_TABLE_FUSION.values()
        + VIRTUAL_FILENAMES_TABLE_CN.values()
        + VIRTUAL_FILENAMES_TABLE_MUT.values()
    ]

    tc = TaigaClient()
    tc.update_dataset(
        eternal_id,
        changes_description="new " + samplesetname + " omics dataset.",
        add_taiga_ids=[
            {"taiga_id": latest_version + "/" + file, "name": file} for file in files
        ],
        add_all_existing_files=True,
    )


def CCLEupload(taiga_ids=""):
    if taiga_ids == "":
        taiga_ids = initVirtualDatasets(samplesetname=SAMPLESETNAME)

    makePRLvMatrices(virtual_ids=taiga_ids)
    makeModelLvMatrices(virtual_ids=taiga_ids)
    uploadAuxTables(taiga_ids=taiga_ids)
    updateEternal(virtual=taiga_ids)
