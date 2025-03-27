from __future__ import print_function
import pandas as pd
from datetime import date

from depmap_omics_upload.mgenepy.utils import helper as h
from depmap_omics_upload import tracker as track
from taigapy import TaigaClient
import json
import pkgutil

# loading config

configdata = pkgutil.get_data(__name__, "config.json")
config = json.loads(configdata)  # type: ignore

config["latest2fn_nummat_model"] = {
    config["taiga_cn"]: config["virtual_filenames_nummat_cn_model"],
    config["taiga_expression"]: config["virtual_filenames_nummat_exp_model"],
    config["taiga_mutation"]: config["virtual_filenames_nummat_mut_model"],
}

config["latest2fn_table_model"] = {
    config["taiga_cn"]: config["virtual_filenames_table_cn_model"],
    config["taiga_fusion"]: config["virtual_filenames_table_fusion_model"],
    config["taiga_mutation"]: config["virtual_filenames_table_mut_model"],
}

config["latest2fn_nummat_pr"] = {
    config["taiga_cn"]: config["virtual_filenames_nummat_cn_pr"],
    config["taiga_expression"]: config["virtual_filenames_nummat_exp_pr"],
    config["taiga_mutation"]: config["virtual_filenames_nummat_mut_pr"],
}

config["latest2fn_table_pr"] = {
    config["taiga_cn"]: config["virtual_filenames_table_cn_pr"],
    config["taiga_fusion"]: config["virtual_filenames_table_fusion_pr"],
    config["taiga_mutation"]: config["virtual_filenames_table_mut_pr"],
}
config["latest2fn_raw_pr"] = {
    config["taiga_mutation"]: config["virtual_filenames_raw_mut_pr"]
}


def getPRToRelease(
    today=None, portals=config["datasets"], date_col_dict=config["date_col_dict"]
):
    """generate lists of profiles to release based on date for all portals

    Returns:
        prs (dict{(portal: list of PRs)}): for each portal, list of profile IDs
    """
    if today is None:
        today = date.today()
    mytracker = track.SampleTracker()
    pr_table = mytracker.read_pr_table()
    pr_table = pr_table[
        (pr_table.BlacklistOmics != True)
        & (pr_table.Datatype.isin(["rna", "wgs", "wes"]))
    ]
    prs = dict()
    for k, v in date_col_dict.items():
        prs_with_date = pr_table[~(pr_table[v].isnull())]
        if k in portals:
            prs_to_release = prs_with_date[
                (prs_with_date[v] <= today) & (prs_with_date.ProfileSource != "taiga")
            ]
            if prs_to_release.MainSequencingID.isnull().values.any():
                raise Exception(
                    "No main SequencingID associated with the following profiles:"
                    + str(
                        set(
                            prs_to_release[
                                prs_to_release.MainSequencingID.isnull()
                            ].index
                        )
                    )
                    + ". Contact Ops to update release dates"
                )
            prs[k] = prs_with_date[
                (prs_with_date[v] <= today)
                & (prs_with_date.ProfileSource != "taiga")
                & (~prs_with_date.MainSequencingID.isnull())
            ].index.tolist()
    if "dmc" in portals and "internal" in portals:
        assert (
            len(set(prs["dmc"]) - set(prs["internal"])) == 0
        ), "Lines with DMC release dates missing internal release dates: " + str(
            set(prs["dmc"]) - set(prs["internal"])
        )
    return prs


def checkDataPermission(today=None, datecol=config["date_col_dict"]["internal"]):
    """confirm that all profiles slated as part of the release all have
    permission to release on the model level

    should be part of gumbo's sanity check but we should run it to be extra safe.

    doesn't rely on availability of data so can be run at any point."""
    mytracker = track.SampleTracker()
    pr_table = mytracker.add_model_cols_to_prtable(
        cols=["ModelID", "PermissionToRelease"]
    )
    problematic_prs = pr_table[
        ~(pr_table[datecol].isnull()) & (~pr_table["PermissionToRelease"])
    ]
    if len(problematic_prs) > 0:
        raise Exception(
            "The following profiles do not have permission to release: "
            + str(set(problematic_prs.index))
            + "The corresponding ModelIDs are: "
            + str(set(problematic_prs.ModelID))
        )


def checkDataAvailability(
    today=None,
    exptaigaid=config["taiga_expression"],
    exptaigafn="proteinCoding_genes_tpm_logp1_profile",
    cntaigaid=config["taiga_cn"],
    cntaigafn="merged_gene_cn_profile",
):
    """confirm that all profiles that are part of the release actually have
    valid data. If not, notify ops so they can update release dates accordingly.

    we only need to check internal PRs because it should be a superset of all portals.

    should be run after data generation for the corresponding release.

    """
    tc = TaigaClient()
    prs = getPRToRelease(today=today)
    mytracker = track.SampleTracker()
    pr_table = mytracker.read_pr_table()
    all_rna_prs = set(pr_table[pr_table.Datatype == "rna"].index)
    all_dna_prs = set(pr_table[pr_table.Datatype.isin(["wgs", "wes"])].index)

    exp = tc.get(name=exptaigaid, file=exptaigafn)
    exp_prs_avail = set(exp.index)
    cn = tc.get(name=cntaigaid, file=cntaigafn)
    cn_prs_avail = set(cn.index)
    unavail_rna = set(prs["internal"]).intersection(all_rna_prs) - exp_prs_avail
    unavail_dna = set(prs["internal"]).intersection(all_dna_prs) - cn_prs_avail

    print("No data available for the following RNAseq profiles: ", unavail_rna)
    print("No data available for the following WES/WGS profiles: ", unavail_dna)

    assert (
        len(unavail_rna) == 0
    ), "missing data in omics latest datasets for profiles above"
    assert (
        len(unavail_dna) == 0
    ), "missing data in omics latest datasets for profiles above"


def makeAchillesChoiceTable(
    prs,
    source_priority=config["source_priority"],
    colnames=config["ach_choice_table_cols"],
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
            pr = subset_pr_table[
                subset_pr_table.MainSequencingID == latest_cds_id
            ].index[0]
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
            cds_ids_wgs = prs_in_mc[
                prs_in_mc.Datatype == "wgs"
            ].MainSequencingID.tolist()
            cds_ids_wes = prs_in_mc[
                prs_in_mc.Datatype == "wes"
            ].MainSequencingID.tolist()
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
                pr = subset_pr_table[
                    subset_pr_table.MainSequencingID == latest_cds_id_wes
                ].index[0]
            # if this MC has wgs PR(s)
            else:
                # ignore the wes PRs, select the PR from the most prioritized source
                # according to the ranking in source_priority from the wgs PR(s)
                subset_seq_table = seq_table[seq_table.index.isin(cds_ids_wgs)]
                subset_seq_table.source = subset_seq_table.source.replace(
                    source_priority
                )
                latest_cds_id_wgs = subset_seq_table.loc[cds_ids_wgs, "source"].idxmin()
                pr = subset_pr_table[
                    subset_pr_table.MainSequencingID == latest_cds_id_wgs
                ].index[0]
            rows.append((mc, pr, "dna"))
    ach_table = pd.DataFrame(rows, columns=colnames)

    return ach_table


def makeDefaultModelTable(
    prs,
    source_priority=config["source_priority"],
    colnames=config["default_table_cols"],
):
    """generate a table that indicates which profiles are selected to represent each Model

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
    # add drug column from MC table to seq table
    seq_table["ModelCondition"] = seq_table["ProfileID"].map(
        dict(zip(pr_table.index, pr_table.ModelCondition))
    )
    seq_table["Drug"] = seq_table["ModelCondition"].map(
        dict(zip(mc_table.index, mc_table.Drug))
    )
    source_priority = {source_priority[i]: i for i in range(len(source_priority))}
    rows = []
    subset_pr_table = pr_table.loc[prs]
    subset_pr_table = subset_pr_table[subset_pr_table.BlacklistOmics != 1]
    mcs = set(subset_pr_table["ModelCondition"])
    # MCs treated with a Drug that's not DMSO can never be considered default basel omics
    # If there's no non-drug treated MC for a model, the model should not have any default profiles
    nondrugged_mc_table = mc_table[
        (mc_table.index.isin(mcs))
        & ((mc_table.Drug.isna()) | (mc_table.Drug == "DMSO"))
    ]
    models = set(mc_table.loc[list(mcs), "ModelID"])
    # we're only picking one PR per datatype (rna/dna) for each Model
    for m in models:
        subset_mc_table = nondrugged_mc_table[nondrugged_mc_table.ModelID == m]
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
            cds_ids = prs_in_model[
                prs_in_model.Datatype == "rna"
            ].MainSequencingID.tolist()
            subset_seq_table = seq_table[seq_table.index.isin(cds_ids)].copy()
            subset_seq_table["source"] = subset_seq_table["source"].map(source_priority)
            # None > DMSO > drug. If Drug == DMSO, it should be ranked lower than Drug == NaN
            # for example, a Sanger no-drug should have higher priority than Broad DMSO
            # add 100 to make sure it overrides source priority
            subset_seq_table.loc[subset_seq_table.Drug == "DMSO", "source"] += 100
            latest_cds_id = subset_seq_table.loc[cds_ids, "source"].idxmin()
            pr = subset_pr_table[
                subset_pr_table.MainSequencingID == latest_cds_id
            ].index[0]
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
            cds_ids_dna = prs_in_model[
                (prs_in_model.Datatype == "wgs")
                | (
                    (prs_in_model.Datatype == "wes")
                    & (prs_in_model.MainSequencingID != "")
                )
            ].MainSequencingID.tolist()  # MainSequencingID is '' when the profile is in legacy
            pr = ""
            # 3 things to look at to assign priority, ranked from most important to least important
            # DMSO vs. Drug == NaN, wes vs. wgs, and source
            # for example, a sanger drug==NaN WES gets picked over broad DMSO WGS
            subset_seq_table = seq_table[seq_table.index.isin(cds_ids_dna)].copy()
            subset_seq_table["source"] = subset_seq_table["source"].map(source_priority)
            subset_seq_table.loc[subset_seq_table.Drug == "DMSO", "source"] += 100
            subset_seq_table.loc[
                subset_seq_table.expected_type == "wes", "source"
            ] += 50
            latest_cds_id_dna = subset_seq_table.loc[cds_ids_dna, "source"].idxmin()
            pr = subset_pr_table[
                subset_pr_table.MainSequencingID == latest_cds_id_dna
            ].index[0]
            rows.append((m, pr, "dna"))
    default_table = pd.DataFrame(rows, columns=colnames)
    return default_table


def makeProfileTable(prs, columns=config["profile_table_cols"]):
    """subset gumbo profile table both column- and row-wise for the release

    Args:
        prs (list): list of profile IDs to be released
        columns (list, optional): list of columns to be included in the table

    Returns:
        pr_table (pd.DataFrame): a df containing omics profile information
    """
    mytracker = track.SampleTracker()
    pr_table = mytracker.add_model_cols_to_prtable(["ModelID"])
    seq_table = mytracker.read_seq_table()
    mc_table = mytracker.read_mc_table()
    pr_table["Stranded"] = pr_table["MainSequencingID"].map(
        dict(zip(seq_table.index, seq_table.stranded))
    )
    pr_table["Source"] = pr_table["ModelCondition"].map(
        dict(zip(mc_table.index, mc_table.Source))
    )
    pr_table = pr_table.loc[prs, columns].rename(
        columns={
            "Baits": "WESKit",
            "actual_seq_technology": "Product",
            "shared_to_dbgap": "SharedToDbgap",
        }
    )
    pr_table = pr_table[pr_table["Datatype"].isin(["rna", "wgs", "wes", "SNParray"])]
    pr_table.loc[
        pr_table[pr_table["Datatype"].isin(["rna", "wgs", "SNParray"])].index.tolist(),
        "WESKit",
    ] = ""

    return pr_table


def initVirtualDatasets(
    samplesetname, taiga_folder_id=config["virtual_folder"], portals=config["datasets"]
):
    """initialize both PR- and Model-level taiga virtual datasets for all 4 portals by uploading an empty dummy file"""
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
    folder=config["working_dir"],
    change_desc="",
    save_format=".csv",
    save_sep=",",
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

    if "EntrezGeneID" in set(to_subset.columns):
        print("making sure Entrez column is Int64")
        to_subset["EntrezGeneID"] = to_subset["EntrezGeneID"].fillna(0)
        to_subset["EntrezGeneID"] = (
            to_subset["EntrezGeneID"].astype("Int64").astype(str)
        )
        to_subset["EntrezGeneID"] = to_subset["EntrezGeneID"].replace({"0": ""})

    print("subsetting ", latest_fn)
    if pr_col == "index":
        subset_mat = to_subset[to_subset.index.isin(prs)]
        subset_mat.to_csv(folder + virtual_fn + save_format)
    elif pr_col == "Tumor_Sample_Barcode":
        subset_mat = to_subset[to_subset[pr_col].isin(prs)]
        subset_mat.to_csv(folder + virtual_fn + save_format, sep=save_sep, index=False)
    else:
        subset_mat = to_subset[to_subset[pr_col].isin(prs)]
        subset_mat = subset_mat.rename(columns={pr_col: "ProfileID"})
        subset_mat.to_csv(folder + virtual_fn + save_format, sep=save_sep, index=False)

    print("uploading ", virtual_fn, " to virtual")
    tc.update_dataset(
        dataset_id=taiga_virtual,
        changes_description=change_desc,
        upload_files=[
            {
                "path": folder + virtual_fn + save_format,
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
    folder=config["working_dir"],
    change_desc="",
    sampleid=config["sample_id"],
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

    if "EntrezGeneID" in set(to_subset.columns):
        print("making sure Entrez column is Int64")
        to_subset["EntrezGeneID"] = to_subset["EntrezGeneID"].fillna(0)
        to_subset["EntrezGeneID"] = (
            to_subset["EntrezGeneID"].astype("Int64").astype(str)
        )
        to_subset["EntrezGeneID"] = to_subset["EntrezGeneID"].replace({"0": ""})

    print("subsetting ", latest_fn)
    if pr_col == "index":
        subset_mat = to_subset[to_subset.index.isin(set(pr2model_dict.keys()))].rename(
            index=pr2model_dict
        )
        subset_mat.to_csv(folder + virtual_fn + ".csv")
    else:
        subset_mat = to_subset[
            to_subset[pr_col].isin(set(pr2model_dict.keys()))
        ].replace({sampleid: pr2model_dict})
        subset_mat = subset_mat.rename(columns={sampleid: "ModelID"})
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
    pr2model_dict,
    portal,
    taiga_latest=config["taiga_mutation"],
    fn_mapping=config["virtual_filenames_guidemut"],
    taiga_virtual="",
    working_dir=config["working_dir"],
):
    """subset, rename, save and upload to taiga germline binary matrix

    Args:
        pr2model_dict (dict): renaming scheme mapping from PR-id to model id
        portal (str): which portal the data is being uploaded to
        taiga_latest (str): which dataset the matrices to be subsetted are being read from
        taiga_virtual (str): which dataset the matrices are being uploaded to
    """
    folder = working_dir + portal + "/model/"
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


def uploadMSRepeatProfile(
    prs,
    taiga_virtual,
    taiga_latest=config["taiga_cn"],
    fn_mapping=config["virtual_filenames_ms_repeat_pr"],
    folder=config["working_dir"],
    save_format=".csv",
    save_sep=",",
    num_static_cols=5,
):
    """subset, save and upload to taiga PR-level microsatellite repeat matrix

    Args:
        prs (list): list of PR-ids to release
        taiga_latest (str): which dataset the matrices to be subsetted are being read from
        taiga_virtual (str): which dataset the matrices are being uploaded to
        fn_mapping (dict): dictionary where key = filename in latest and value = filename in virtual
        folder (str): where the file should be stores before uploading to virtual
        num_static_cols (int): number of columns in the df that are static/not profiles
    """
    for latest_fn, virtual_fn in fn_mapping.items():
        print("loading ", latest_fn, " from latest")
        tc = TaigaClient()
        to_subset = tc.get(name=taiga_latest, file=latest_fn)

        print("subsetting ", latest_fn)
        subset_mat = to_subset.iloc[:, :num_static_cols].join(
            to_subset.iloc[:, num_static_cols:][list(set(to_subset.columns) & set(prs))]
        )
        subset_mat.to_csv(folder + virtual_fn + save_format, sep=save_sep, index=False)

        print("uploading ", virtual_fn, " to virtual")
        tc.update_dataset(
            dataset_id=taiga_virtual,
            changes_description="adding " + virtual_fn,
            upload_files=[
                {
                    "path": folder + virtual_fn + save_format,
                    "name": virtual_fn,
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
            ],
            add_all_existing_files=True,
        )


def uploadAuxTables(
    taiga_ids,
    ach_table_name=config["ach_choice_table_name"],
    default_table_name=config["default_table_name"],
    release_pr_table_name=config["release_pr_table_name"],
    folder=config["working_dir"] + config["sampleset"] + "/",
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
        profile_table.to_csv(folder + portal + "_" + release_pr_table_name + ".csv")
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
                    "path": folder
                    + "/"
                    + portal
                    + "_"
                    + release_pr_table_name
                    + ".csv",
                    "name": "OmicsProfiles",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
            ],
            add_all_existing_files=True,
        )


def makePRLvMatrices(
    virtual_ids,
    files_nummat=config["latest2fn_nummat_pr"],
    folder=config["working_dir"] + config["sampleset"],
    files_table=config["latest2fn_table_pr"],
    files_raw=config["latest2fn_raw_pr"],
    today=None,
    sampleid=config["sample_id"],
    exclude=config["exclude"],
):
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
                if latest not in exclude[portal]:
                    uploadPRMatrix(
                        prs_to_release,
                        latest_id,
                        virtual_ids[portal],
                        latest,
                        virtual,
                        "NumericMatrixCSV",
                        pr_col="index",
                        folder=folder + "/",
                        change_desc="adding " + virtual,
                    )
        for latest_id, fn_dict in files_table.items():
            for latest, virtual in fn_dict.items():
                if latest not in exclude[portal]:
                    uploadPRMatrix(
                        prs_to_release,
                        latest_id,
                        virtual_ids[portal],
                        latest,
                        virtual,
                        "TableCSV",
                        pr_col=sampleid,
                        folder=folder + "/",
                        change_desc="adding " + virtual,
                    )
        for latest_id, fn_dict in files_raw.items():
            for latest, virtual in fn_dict.items():
                if latest not in exclude[portal]:
                    uploadPRMatrix(
                        prs_to_release,
                        latest_id,
                        virtual_ids[portal],
                        latest,
                        virtual,
                        "Raw",
                        pr_col="Tumor_Sample_Barcode",
                        folder=folder + "/",
                        change_desc="adding " + virtual,
                        save_format=".maf",
                        save_sep="\t",
                    )
        uploadMSRepeatProfile(prs_to_release, virtual_ids[portal], folder=folder + "/")


def makeModelLvMatrices(
    virtual_ids,
    folder=config["working_dir"] + config["sampleset"],
    files_nummat=config["latest2fn_nummat_model"],
    files_table=config["latest2fn_table_model"],
    upload_guide_matrices=True,
    today=None,
    sampleid=config["sample_id"],
    exclude=config["exclude"],
    omics_id_mapping_table_name="",
):
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
                if latest not in exclude[portal]:
                    uploadModelMatrix(
                        pr2model_dict,
                        latest_id,
                        virtual_ids[portal],
                        latest,
                        virtual,
                        "NumericMatrixCSV",
                        pr_col="index",
                        folder=folder + "/",
                        change_desc="adding " + virtual,
                    )
        for latest_id, fn_dict in files_table.items():
            for latest, virtual in fn_dict.items():
                if latest not in exclude[portal]:
                    uploadModelMatrix(
                        pr2model_dict,
                        latest_id,
                        virtual_ids[portal],
                        latest,
                        virtual,
                        "TableCSV",
                        pr_col=sampleid,
                        folder=folder + "/",
                        change_desc="adding " + virtual,
                    )
        if upload_guide_matrices:
            uploadBinaryGuideMutationMatrixModel(
                pr2model_dict, portal, taiga_virtual=virtual_ids[portal]
            )

def makeWESandWGSMatrices(virtual_ids,
    folder=config["working_dir"] + config["sampleset"],
    files_nummat_model={config["taiga_cn"]: config['virtual_filenames_nummat_cn_model_wgs_and_wes']},
    files_nummat_pr={config["taiga_cn"]: config['virtual_filenames_nummat_cn_pr_wgs_and_wes']},
    files_table_pr={config["taiga_cn"]: config['virtual_filenames_table_cn_pr_wgs_and_wes']},
    today=None,
    sampleid=config["sample_id"],
    exclude=config["exclude"],
    omics_id_mapping_table_name="",
    ):
    """for each portal, save and upload data matrices, each matrix has a WGS and a WES

    Args:
        taiga_ids (dict): dictionary that maps portal name to virtual taiga dataset id

    Returns:
        prs (dict{(portal: list of PRs)}): for each portal, list of profile IDs
    """
    tc = TaigaClient()
    for portal, taiga_id in virtual_ids.items():
        omics_id_mapping_table = tc.get(name=virtual_ids[portal], file=omics_id_mapping_table_name)
        prs_to_release_wes = omics_id_mapping_table[omics_id_mapping_table.datatype == 'wes'].profile_id.tolist()
        prs_to_release_wgs = omics_id_mapping_table[omics_id_mapping_table.datatype == 'wgs'].profile_id.tolist()
        
        default_table_wes = omics_id_mapping_table[(omics_id_mapping_table['is_default_entry'] == True) & (omics_id_mapping_table.datatype == 'wes')]
        default_table_wgs = omics_id_mapping_table[(omics_id_mapping_table['is_default_entry'] == True) & (omics_id_mapping_table.datatype == 'wgs')]
        
        pr2model_dict_wes = dict(list(zip(default_table_wes.profile_id, default_table_wes.model_id)))
        pr2model_dict_wgs = dict(list(zip(default_table_wgs.profile_id, default_table_wgs.model_id)))
        print("uploading respective WES/WGS matrices to", portal)
        for latest_id, fn_dict in files_nummat_model.items():
            for latest, virtual in fn_dict.items():
                if latest not in exclude[portal]:
                    uploadModelMatrix(
                        pr2model_dict_wes,
                        latest_id,
                        taiga_id,
                        latest,
                        virtual+"WES",
                        "NumericMatrixCSV",
                        pr_col="index",
                        folder=folder + "/",
                        change_desc="adding " + virtual+"WES",
                    )
                    uploadModelMatrix(
                        pr2model_dict_wgs,
                        latest_id,
                        taiga_id,
                        latest,
                        virtual+"WGS",
                        "NumericMatrixCSV",
                        pr_col="index",
                        folder=folder + "/",
                        change_desc="adding " + virtual+"WGS",
                    )
        for latest_id, fn_dict in files_nummat_pr.items():
            for latest, virtual in fn_dict.items():
                if latest not in exclude[portal]:
                    uploadPRMatrix(
                        prs_to_release_wes,
                        latest_id,
                        taiga_id,
                        latest,
                        virtual+"WES",
                        "NumericMatrixCSV",
                        pr_col="index",
                        folder=folder + "/",
                        change_desc="adding " + virtual+"WES",
                    )
                    uploadPRMatrix(
                        prs_to_release_wgs,
                        latest_id,
                        taiga_id,
                        latest,
                        virtual+"WGS",
                        "NumericMatrixCSV",
                        pr_col="index",
                        folder=folder + "/",
                        change_desc="adding " + virtual+"WGS",
                    )
        for latest_id, fn_dict in files_table_pr.items():
            for latest, virtual in fn_dict.items():
                if latest not in exclude[portal]:
                    uploadPRMatrix(
                        prs_to_release_wes,
                        latest_id,
                        taiga_id,
                        latest,
                        virtual+"WES",
                        "TableCSV",
                        pr_col=sampleid,
                        folder=folder + "/",
                        change_desc="adding " + virtual+"WES",
                    )
                    uploadPRMatrix(
                        prs_to_release_wgs,
                        latest_id,
                        taiga_id,
                        latest,
                        virtual+"WGS",
                        "TableCSV",
                        pr_col=sampleid,
                        folder=folder + "/",
                        change_desc="adding " + virtual+"WGS",
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
    eternal_id=config["taiga_eternal"],
    virtual=config["virtual"],
    samplesetname=config["sampleset"],
):
    """update taiga eternal dataset by linking to latest virtual internal dataset"""
    latest_version = findLatestVersion(virtual["internal"])

    files = (
        list(config["virtual_filenames_nummat_exp_model"].values())
        + list(config["virtual_filenames_nummat_exp_pr"].values())
        + list(config["virtual_filenames_nummat_cn_model"].values())
        + list(config["virtual_filenames_nummat_cn_pr"].values())
        + list(config["virtual_filenames_nummat_mut_model"].values())
        + list(config["virtual_filenames_guidemut"].values())
        + list(config["virtual_filenames_table_fusion_model"].values())
        + list(config["virtual_filenames_table_fusion_pr"].values())
        + list(config["virtual_filenames_table_cn_pr"].values())
        + list(config["virtual_filenames_table_mut_model"].values())
        + list(config["virtual_filenames_table_mut_pr"].values())
        + list(config["virtual_filenames_raw_mut_pr"].values())
        + list(config["virtual_filenames_ms_repeat_pr"].values())
    )

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
        taiga_ids = initVirtualDatasets(samplesetname=config["sampleset"])

    makePRLvMatrices(virtual_ids=taiga_ids)
    makeModelLvMatrices(virtual_ids=taiga_ids)
    uploadAuxTables(taiga_ids=taiga_ids)
    updateEternal(virtual=taiga_ids)
