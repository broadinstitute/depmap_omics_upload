# tracker.py
from depmap_omics_upload.mgenepy.utils import helper as h
import numpy as np
import os

from depmap_omics_upload.mgenepy import terra
from depmap_omics_upload.mgenepy.google import gcp
import dalmatian as dm
import signal
import gumbo_rest_client
import gumbo_rest_client.utils as gumbo_utils
from datetime import date
import json
import pkgutil

# loading config
configdata = pkgutil.get_data(__name__, "config.json")
config = json.loads(configdata)  # type: ignore


# condense all interactions with tracker (for emeril integration)
class SampleTracker:
    """
    interacts with (read + write) the sample tracker gsheet
    """

    # pylint: disable=too-many-instance-attributes

    def __init__(self, gumbo_env="production"):
        self.model_table_name = config["model_table_name"]
        self.mc_table_name = config["mc_table_name"]
        self.pr_table_name = config["pr_table_name"]
        self.seq_table_name = config["seq_table_name"]
        self.model_table_index = config["model_table_index"]
        self.mc_table_index = config["mc_table_index"]
        self.pr_table_index = config["pr_table_index"]
        self.seq_table_index = config["seq_table_index"]
        self.sample_table_name = config["sample_table_name"]
        self.screen_table_name = config["screen_table_name"]
        self.screen_table_index = config["screen_table_index"]
        self.str_table_name = config["str_table_name"]
        self.str_table_index = config["str_table_index"]
        self.client = gumbo_rest_client.Client(
            authed_session=gumbo_rest_client.create_authorized_session(
                use_default_service_account=False
            ),
            username="depmap_omics_upload",
            base_url=gumbo_rest_client.const.prod_url
            if gumbo_env == "production"
            else gumbo_rest_client.const.staging_url,
        )
        self.mapping_utils = gumbo_utils.NameMappingUtils()

    def read_model_table(self):
        model_table = self.client.get(self.model_table_name)
        model_table_camel = self.mapping_utils.rename_columns(
            self.model_table_name, model_table
        )
        model_table_camel = model_table_camel.set_index(self.model_table_index)
        return model_table_camel

    def read_mc_table(self):
        mc_table = self.client.get(self.mc_table_name)
        mc_table_camel = self.mapping_utils.rename_columns(self.mc_table_name, mc_table)
        mc_table_camel = mc_table_camel.set_index(self.mc_table_index)
        return mc_table_camel

    def read_pr_table(self):
        pr_table = self.client.get(self.pr_table_name)
        pr_table_camel = self.mapping_utils.rename_columns(self.pr_table_name, pr_table)
        pr_table_camel = pr_table_camel.set_index(self.pr_table_index)
        return pr_table_camel

    def read_seq_table(self):
        seq_table = self.client.get(self.seq_table_name)
        seq_table_camel = self.mapping_utils.rename_columns(
            self.seq_table_name, seq_table
        )
        seq_table_camel = seq_table_camel.set_index(self.seq_table_index)
        return seq_table_camel

    def read_screen_table(self):
        screen_table = self.client.get(self.screen_table_name)
        screen_table_camel = self.mapping_utils.rename_columns(
            self.screen_table_name, screen_table
        )
        screen_table_camel = screen_table_camel.set_index(self.screen_table_index)
        return screen_table_camel

    def read_str_table(self):
        str_table = self.client.get(self.str_table_name)
        str_table = str_table.set_index(self.str_table_index)
        return str_table

    def write_mc_table(self, df):
        # assumes df's columns are camelCase, and converts it back to snake_case
        df = df.reset_index(level=0)
        self.mapping_utils.rename_columns(
            self.mc_table_name, df, convert_to_custom_names=False, inplace=True
        )
        self.client.update_only(self.mc_table_name, df)

    def write_pr_table(self, df):
        # assumes df's columns are camelCase, and converts it back to snake_case
        df = df.reset_index(level=0)
        self.mapping_utils.rename_columns(
            self.pr_table_name, df, convert_to_custom_names=False, inplace=True
        )
        self.client.update_only(self.pr_table_name, df)

    def write_seq_table(self, df):
        # assumes df's columns are camelCase, and converts it back to snake_case
        df.index.name = self.seq_table_index
        df = df.reset_index(level=0)
        self.mapping_utils.rename_columns(
            self.seq_table_name, df, convert_to_custom_names=False, inplace=True
        )
        self.client.update_only(self.seq_table_name, df)

    def write_str_table(self, df):
        df.index.name = self.str_table_index
        df = df.reset_index(level=0)
        self.client.update_only(self.str_table_name, df)

    def insert_to_seq_table(self, df):
        df.index.name = self.seq_table_index
        df = df.reset_index(level=0)
        self.mapping_utils.rename_columns(
            self.seq_table_name, df, convert_to_custom_names=False, inplace=True
        )
        self.client.insert_only(self.seq_table_name, df)

    def insert_to_str_table(self, df):
        df.index.name = self.str_table_index
        df = df.reset_index(level=0)
        self.client.insert_only(self.str_table_name, df)

    def get_participant_id(self, seqid, seq_table, pr_table, mc_table, model_table):
        # assumes all tables are camelCase
        pr = seq_table.loc[seqid, self.pr_table_index]
        if pr != "":
            mc = pr_table.loc[pr, self.mc_table_index]
            model = mc_table.loc[mc, self.model_table_index]
            participant = model_table.loc[model, "PatientID"]
            return participant
        else:
            return "NA"

    def add_model_cols_to_seqtable(self, cols=["ModelID"]):
        # add columns from model table to seq table
        seq_table = self.read_seq_table()
        pr_table = self.read_pr_table()
        mc_table = self.read_mc_table()
        model_table = self.read_model_table().reset_index(level=0)
        for c in cols:
            assert c in model_table.columns, c + " is not a column in model table"
        for seq_id in seq_table.index:
            pr = seq_table.loc[seq_id, self.pr_table_index]
            if pr is not None:
                mc = pr_table.loc[pr, "ModelCondition"]
                model = mc_table.loc[mc, self.model_table_index]
                for c in cols:
                    seq_table.loc[seq_id, c] = model_table[
                        model_table[self.model_table_index] == model
                    ][c].values[0]
            else:
                for c in cols:
                    seq_table.loc[seq_id, c] = "nan"
        return seq_table

    def add_model_cols_to_prtable(self, cols=["ModelID"]):
        # add columns from model table to seq table
        pr_table = self.read_pr_table()
        pr_table = pr_table[~pr_table.ModelCondition.isnull()]
        mc_table = self.read_mc_table()
        model_table = self.read_model_table().reset_index(level=0)
        for c in cols:
            assert c in model_table.columns, c + " is not a column in model table"
        for pr_id in pr_table.index:
            if pr_id is not None:
                mc = pr_table.loc[pr_id, "ModelCondition"]
                model = mc_table.loc[mc, self.model_table_index]
                for c in cols:
                    pr_table.loc[pr_id, c] = model_table[
                        model_table[self.model_table_index] == model
                    ][c].values[0]
            else:
                for c in cols:
                    pr_table.loc[pr_id, c] = "nan"
        return pr_table

    def lookup_model_from_pr(self, pr_id, model_col):
        # given a profile ID, look for a column in the model table
        pr_table = self.read_pr_table()
        mc_table = self.read_mc_table()
        model_table = self.read_model_table().reset_index(level=0)
        mc = pr_table.loc[pr_id, "ModelCondition"]
        model = mc_table.loc[mc, self.model_table_index]
        return model_table[model_table[self.model_table_index] == model][
            model_col
        ].values[0]

    def update_pr_from_seq(
        self,
        datatype,
        cols={
            "bam_public_sra_path": "BamPublicSRAPath",
            "blacklist": "BlacklistOmics",
            "issue": "Issue",
            "prioritized": "Prioritized",
        },
        priority=None,
        dryrun=False,
    ):
        seq_table = self.read_seq_table()
        seq_table = seq_table[
            (seq_table.blacklist != True)
            & (seq_table.expected_type.isin(datatype))
            & (seq_table.processed_sequence == True)
        ]
        pr_table = self.read_pr_table()
        prs_in_seq_table = seq_table.ProfileID.unique()

        cds2pr_dict = {}
        for pr in prs_in_seq_table:
            if len(seq_table[seq_table.ProfileID == pr]) == 1:
                pr_table.loc[pr, "MainSequencingID"] = seq_table[
                    seq_table.ProfileID == pr
                ].index
            else:
                allv = seq_table[seq_table["ProfileID"] == pr]
                for k, val in allv.iterrows():
                    if priority is None:
                        if val["version"] == max(allv.version.values):
                            cds2pr_dict[k] = pr
                            break
                    else:
                        if val["version"] == max(allv.version.values):
                            cds2pr_dict[k] = pr
                        if val["priority"] == 1:
                            cds2pr_dict[k] = pr
                            break
        for k, v in cds2pr_dict.items():
            pr_table.loc[v, "MainSequencingID"] = k
        if not dryrun:
            self.write_pr_table(pr_table)
        return pr_table

    def shareCCLEbams(
        self,
        samples,
        users=[],
        groups=[],
        raise_error=True,
        arg_max_length=100000,
        bamcols=["bam_filepath", "bai_filepath"],
        unshare=False,
        requesterpays_project="",
    ):
        """
        same as shareTerraBams but is completed to work with CCLE bams from the CCLE sample tracker

        Args:
        ----
            users: list[str] of users' google accounts
            groups: list[str] of groups' google accounts
            samples list[str] of samples cds_ids for which you want to share data
            bamcols: list[str] list of column names where bams/bais are
            raise_error: whether or not to raise an error if we find blacklisted lines
            refsheet_url: the google spreadsheet where the samples are stored
            privacy_sheeturl: the google spreadsheet where the samples are stored
            requesterpays_project: the google project where the requester pays bucket is located

        Returns:
        --------
            a list of the gs path we have been giving access to
        """
        refdata = self.read_seq_table()
        pr_table = self.read_pr_table()
        blacklist = [i for i in refdata["blacklist"].values.tolist() if i == 1]
        blacklisted = set(blacklist) & set(samples)
        print("we have " + str(len(blacklisted)) + " files blacklisted for low quality")
        if len(blacklisted):
            print("these lines are blacklisted " + str(blacklisted))
            if raise_error:
                raise ValueError("blacklistedlines")
        if type(users) is str:
            users = [users]

        embargoed = []
        for s in samples:
            pr_id = refdata.loc[s, "ProfileID"]
            release_date = pr_table.loc[pr_id, "InternalReleaseDate"]
            today = date.today()
            if release_date == "" or release_date > today:
                embargoed.append(s)

        if len(embargoed) > 0:
            print(
                "the following lines are currently under embargo, can't share the files yet!"
            )
            print(embargoed)
            samples = [s for s in samples if s not in embargoed]

        togiveaccess = np.ravel(refdata[bamcols].loc[samples].values)
        usrs = ""
        for group in groups:
            usrs += (
                (" -d " if unshare else " -g") + group + (":R" if not unshare else "")
            )
        for user in users:
            usrs += (
                (" -d " if unshare else " -u ") + user + (":R" if not unshare else "")
            )
        requesterpays_cmd = (
            "" if requesterpays_project == "" else "-u " + requesterpays_project
        )
        cmd_prefix = "gsutil {} -m acl ch ".format(requesterpays_cmd) + usrs
        cmd = cmd_prefix
        for n, filename in enumerate(togiveaccess):
            if type(filename) is str and filename:
                oldcmd = cmd
                cmd += " " + filename
                if (len(cmd) > arg_max_length) | (n == len(togiveaccess) - 1):
                    if n < len(togiveaccess) - 1:
                        cmd = oldcmd
                    if unshare:
                        print("preventing access to {:d} files".format(n + 1))
                    else:
                        print("granting access to {:d} files".format(n + 1))
                    with open("/tmp/grantaccess{:d}.sh".format(n), "w") as f:
                        f.write(cmd)
                    code = os.system(cmd)
                    cmd = cmd_prefix + " " + filename
                    if code == signal.SIGINT:
                        print("Awakened")
                        return None

        print("\n\njust install and use gsutil to copy them")
        print("https://cloud.google.com/storage/docs/gsutil_install")
        print("https://cloud.google.com/storage/docs/gsutil/commands/cp")
        return togiveaccess


def findIssue(
    tracker,
    dup=[
        "age",
        "sex",
        "arxspan_id",
        "cellosaurus_id",
        "collection_site",
        "primary_disease",
        "subtype",
        "lineage",
        "stripped_cell_line_name",
    ],
):
    """
    findIssue looks at a couple metrics:

    'things that are from the same patient but don\'t have the same value'
    'things that have duplicate versions'
    'things that don\'t have their legacy bam file'
    'things that don\'t have their bam file path'

    Args:
        tracker (pandas.DataFrame): the tracker
        dup (list, optional): the list of columns to check for duplicates.
        Defaults to ['age', 'sex', 'arxspan_id', 'cellosaurus_id', 'primary_site', 'primary_disease',
        'subtype', 'lineage', 'stripped_cell_line_name']
    """
    print("things that are from the same patient but don't have the same value")
    dup = tracker[dup].set_index("arxspan_id").drop_duplicates()
    print(dup.loc[h.dups(dup.index)])
    print("things that have duplicate versions")
    print(
        h.dups(
            tracker["arxspan_id"]
            + "_"
            + tracker["datatype"]
            + "_"
            + tracker["version"].astype(str)
        )
    )
    print("things that don't have their legacy bam file")
    print(
        tracker[
            tracker["datatype"].isin(["rna", "wes", "wgs"])
            & (
                tracker["legacy_bam_filepath"].isna()
                | tracker["legacy_bai_filepath"].isna()
            )
        ].index
    )
    print("things that don't have their bam file path")
    print(
        tracker[
            (
                tracker["internal_bam_filepath"].isna()
                | tracker["internal_bai_filepath"].isna()
            )
        ].index
    )


def updateFromTracker(
    samples,
    ccle_refsamples,
    arxspan_id="arxspan_id",
    participant_id="participant_id",
    toupdate={},
):
    """update a list of samples' missing information from what is known in the ccle sample tracker

    Args:
        samples (pandas.DataFrame): the samples to update
        ccle_refsamples (pandas.DataFrame): the ccle sample tracker
        arxspan_id (str, optional): the name of the arxspan id column. Defaults to 'arxspan_id'.
        participant_id (str, optional): the name of the participant id column. Defaults to 'participant_id'.
        toupdate (dict(str, []), optional): the columns to update. Defaults to {}.

    Returns:
        (pandas.DataFrame): the updated samples
        (list(str)): the list of samples that were not found in the ccle sample tracker
    """
    # If I have a previous samples I can update unknown data directly
    index = []
    notfound = []
    if len(toupdate) == 0:
        toupdate = {
            "sex": [],
            "primary_disease": [],
            "cellosaurus_id": [],
            "age": [],
            "collection_site": [],
            "subtype": [],
            "subsubtype": [],
            "lineage": [],
            "parent_cell_line": [],
            "matched_normal": [],
            "comments": [],
            # "mediatype": [],
            "condition": [],
            "stripped_cell_line_name": [],
            "participant_id": [],
        }
    # import pdb;pdb.set_trace()
    for k, val in samples.iterrows():
        dat = ccle_refsamples[ccle_refsamples[arxspan_id] == val[arxspan_id]]
        if len(dat) > 0:
            index.append(k)
            for k, v in toupdate.items():
                toupdate[k].append(dat[k].tolist()[0])
        else:
            notfound.append(k)
    # doing so..
    for k, v in toupdate.items():
        samples.loc[index, k] = v
    # Commented out the following because it has no effect
    # len(samples.loc[notfound][participant_id]), samples.loc[notfound][
    #    participant_id
    # ].tolist()
    return samples, notfound


def removeOlderVersions(
    names, refsamples, arxspan_id="arxspan_id", version="version", priority=None
):
    """
    will set it to your sample_ids using the latest version available for each sample

    Given a dataframe containing ids, versions, sample_ids and you dataset df indexed
    by the same ids, will set it to your sample_ids using the latest version
    available for each sample

    Args:
    -----
        refsamples (pd.df): the reference metadata. should contain [id, version, arxspan_id,...]
        names (list[str]): only do it on this set of samples.
        arxspan_id (str, optional): the name of the arxspan id column. Defaults to 'arxspan_id'.
        version (str, optional): the name of the version column. Defaults to 'version'.

    Returns:
    --------
        (dict): the subseted samples

    """
    # pandas throws an error if index is unavailable
    names = [x for x in names if x in refsamples.index.tolist()]

    lennames = len(names)
    res = {}

    refsamples = refsamples.loc[names].copy()
    if lennames > len(refsamples):
        print(set(names) - set(refsamples.index))
        import pdb

        pdb.set_trace()
        raise ValueError(
            "we had some ids in our dataset not registered in this refsample dataframe"
        )
    for arxspan in set(refsamples[arxspan_id]):
        allv = refsamples[refsamples[arxspan_id] == arxspan]
        for k, val in allv.iterrows():
            if priority is None:
                if val[version] == max(allv.version.values):
                    res[k] = arxspan
                    break
            else:
                if val[version] == max(allv.version.values):
                    res[k] = arxspan
                if val[priority] == 1:
                    res[k] = arxspan
                    break

    print("removed " + str(lennames - len(res)) + " duplicate samples")
    # remove all the reference metadata columns except the arxspan ID
    return res


def updateIsogenecity(
    di,
    tracker,
    unset=False,
    toupdate=["participant_id", "age", "sex", "matched_normal"],
):
    """

    Args:
        di ([type]): [description]
        tracker ([type]): [description]
        unset (bool, optional): [description]. Defaults to False.

    Raises:
        ValueError: [description]

    Returns:
        [type]: [description]
    """
    tracker = tracker.copy()
    for k, v in di.items():
        print("________________________________")
        a = tracker[tracker.arxspan_id == k]
        b = tracker[tracker.arxspan_id == v]
        if len(a) == 0:
            print(v, "does not exist")
        if len(b) == 0:
            print(k, "does not exist")
        if len(set(a.participant_id)) > 1 or len(set(b.participant_id)) > 1:
            raise ValueError("not same participant for same cell line")
        if a.participant_id[0] == b.participant_id[0] or (
            unset and a.participant_id[0] != b.participant_id[0]
        ):
            print("already set")
            continue
        print("merging:")
        print(k, v)
        if unset:
            print("changing participant_id of ", v)
            tracker.loc[b.index, "participant_id"] = "PT-" + h.randomString()
        else:
            print("doing:")
            print(a.loc[a.index[0], toupdate].values)
            print("into")
            print(
                tracker.loc[
                    tracker[tracker.participant_id == b.participant_id[0]].index,
                    toupdate,
                ].values
            )
            tracker.loc[
                tracker[tracker.participant_id == b.participant_id[0]].index, toupdate
            ] = a.loc[a.index[0], toupdate].tolist()
    return tracker


def changeCellLineNameInNew(
    ref,
    new,
    datatype,
    dupdict,
    toupdate=[
        "stripped_cell_line_name",
        "arxspan_id",
        "patient_id",
        "sex",
        "primary_disease",
        "cellosaurus_id",
        "age",
        "collection_site",
        "subtype",
        "subsubtype",
    ],
):
    """
    Rename a/some line in a DF and takes care of corresponding metadata and versions from the sample tracker

    !!!! DOES NOT YET WORK !!!! version compute is wrong
    Args:
    -----
        new: change the cell line name in this dataframe
        dupdict: dict(tochange,newname) (arxspan_id:arxspan_id)
        datatype: str for a ref with many datatype (to get the right version number)

    Returns:
    --------
        the updated dataframe
    """
    for k, v in dupdict.items():
        new.loc[new[new.arxspan_id == k].index, toupdate] = ref[ref.arxspan_id == v][
            toupdate
        ].values[0]
        new.loc[new[new.arxspan_id == v].index, "version"] = (
            len(ref[(ref.arxspan_id == v) & (ref.datatype == datatype)]) + 1
        )
    return new


def changeCellLineName(
    tracker,
    datatype,
    dupdict,
    toupdate=[
        "stripped_cell_line_name",
        "participant_id",
        "cellosaurus_id",
        "sex",
        "arxspan_id",
        "parent_cell_line",
        "matched_normal",
        "age",
        "collection_site",
        "primary_disease",
        "subtype",
        "subsubtype",
        "lineage",
    ],
):
    """
    Rename a/some line in our sample tracker and takes care of corresponding metadata and versions from the sample tracker

    Args:
    -----
        dupdict: dict(tochange,newname): the dict of the new name for the cell line: cds-id: arxspan_id
        datatype: str for a tracker with many datatype (to get the right version number)

    Returns:
    --------
        the updated dataframe
    """
    for k, v in dupdict.items():
        try:
            tracker.loc[k, toupdate] = tracker[tracker.arxspan_id == v][
                toupdate
            ].values[0]
            tracker.loc[k, "version"] = len(
                tracker[(tracker.arxspan_id == v) & (tracker.datatype == datatype)]
            )
        except IndexError:
            raise IndexError(str(v) + " not found in tracker")
    return tracker


def cleanVersions(
    tracker,
    samplecol="arxspan_id",
    dryrun=False,
    datatypecol="datatype",
    versioncol="version",
):
    """
    updates and sorts the versions of samples in the sample tracker:

    checks that we get 1,2,3 instead of 2,4,5 when samples are renamed or removed

    Args:
    -----
        tracker: dataframe of the sample tracker
        samplecol: str colname of samples in the trackerr
        dryrun: bool whether or not to apply it or just print if there is issues
        datatypecol: str colname of the datatype values
        versioncol: str colname of the version values

    Returns:
    -------
        an updated sample tracker
    """
    tracker = tracker.copy()
    tracker["samval"] = tracker[samplecol] + tracker[datatypecol]
    for v in set(tracker["samval"]):
        sams = tracker[tracker["samval"] == v]
        vs = sams[versioncol].tolist()
        if max(vs) == len(vs):
            continue
        print("found issue")
        if dryrun:
            continue
        vs.sort()
        rn = {}
        for i, v in enumerate(vs):
            rn.update({v: i + 1})
        for k, val in sams.iterrows():
            tracker.loc[k, versioncol] = rn[val[versioncol]]
    tracker = tracker.drop(columns="samval")
    return tracker


def setRightName(tracker, name="stripped_cell_line_name", signs=["-", "_", ".", " "]):
    """
    cell line name need to be all uppercase and not contain special signs will take carre of ones that are not

    BE CARREFUL, it does not solve issues of duplicate lines (diff sample id, same new name)
    see findLikelyDup for that

    Args:
    -----
        tracker: dataframe of the sample tracker

    Returns:
    -----
        an updated sample tracker
    """
    new = []
    for val in tracker[name]:
        for s in signs:
            val = val.replace(s, "")
        new.append(val.upper())
    tracker[name] = new
    return tracker


def findLikelyDup(
    tracker,
    name="stripped_cell_line_name",
    signs=["-", "_", ".", " "],
    arxspid="arxspan_id",
    looksub=True,
):
    """
    find cell lines that are likely to be duplicates

    will return ,  as well,
    and

    Args:
    -----
        tracker: dataframe of the sample tracker
        looksub: bool, look if a name if within another name (can flag many derivatives)

    Returns:
    --------
        a list[tuples(str,str)] of likly duplicate names as tuples (rh13, RH-13)
        a list[tuples(str,str)] of associated arxspan ids
        a dict[str:set(str)] of arxspan ids that have multiple cell line names associated
    """
    names = set(tracker[name])
    simi = []
    arxsp = []
    issues = {}
    for i, name1 in enumerate(names):
        h.showcount(i, len(names))
        n1 = name1
        for s in signs:
            name1 = name1.replace(s, "")
        name1 = name1.upper()
        for name2 in names - set([n1]):
            n2 = name2
            for s in signs:
                name2 = name2.replace(s, "")
            name2 = name2.upper()
            if name1 == name2:
                if (
                    looksub
                    and (name1 in name2 or name2 in name1)
                    and abs(len(name1) - len(name2)) < 2
                ) or not looksub:
                    if (n1, n2) not in simi and (n2, n1) not in simi:
                        simi.append((n1, n2))
                        arxsp.append(
                            (
                                tracker[tracker[name] == n1][arxspid][0],
                                tracker[tracker[name] == n2][arxspid][0],
                            )
                        )
    for val in set(tracker[name]):
        v = set(tracker[tracker[name] == val][arxspid])
        if len(v) > 1:
            issues.update({val: v})
    return simi, arxsp, issues


def update(
    table,
    samplesetname,
    failed,
    lowqual,
    newgs="",
    refworkspace=None,
    bamfilepaths=["internal_bam_filepath", "internal_bai_filepath"],
    dry_run=False,
    samplesinset=[],
    todrop=[],
    billing_proj=None,
):
    """updates the sample tracker (or Gumbo omicSequencing table) with the new samples and the QC metrics
    Args:
        table (df): [description]
        selected (list[str]): which samples were selected in the release of the analysis
        samplesetname (str): the name of the sample set or of the current analysis
        samplesinset (list[str]): list of samples in the analysis.
        lowqual (list[str]): list of samples that failed QC
        newgs (str, optional): google storage path where to move the files. Defaults to ''.
        sheetcreds (str, optional): google sheet service account file path. Defaults to SHEETCREDS.
        sheetname (str, optional): google sheet service account file path. Defaults to SHEETNAME.
        refworkspace (str, optional): if provideed will extract workspace values (bam files path, QC,...). Defaults to None.
        bamfilepaths (list, optional): Terra columns containing the bam filepath for which to change the location. Defaults to ['internal_bam_filepath', 'internal_bai_filepath'].
    """
    # updating locations of bam files and extracting infos
    if newgs and refworkspace is not None:
        if not samplesinset:
            samplesinset = [
                i["entityName"]
                for i in dm.WorkspaceManager(refworkspace)
                .get_entities("sample_set")
                .loc[samplesetname]
                .samples
            ]
        res, _, _ = terra.changeGSlocation(
            refworkspace,
            newgs=newgs,
            bamfilepaths=bamfilepaths,
            entity="sample",
            keeppath=False,
            dry_run=dry_run,
            onlysamples=samplesinset,
            workspaceto=refworkspace,
        )
        table.loc[res.index.tolist(), ["bam_filepath", "bai_filepath"]] = res[
            bamfilepaths[:2]
        ].values
        table.loc[res.index.tolist(), "bam_size"] = [
            gcp.extractSize(i)[1]
            for i in gcp.lsFiles(
                res[bamfilepaths[0]].tolist(),
                "-l",
                billing_proj=billing_proj,
            )
        ]
        table.loc[res.index.tolist(), "bam_crc32c_hash"] = [
            gcp.extractHash(i)
            for i in gcp.lsFiles(
                res[bamfilepaths[0]].tolist(),
                "-L",
                billing_proj=billing_proj,
            )
        ]

    table.loc[samplesinset, ["low_quality", "blacklist", "prioritized"]] = False
    table.loc[samplesinset, "processed_sequence"] = True
    table.loc[lowqual, "low_quality"] = True
    failed_not_dropped = list(set(failed) - set(todrop))
    table.loc[failed_not_dropped, "blacklist"] = True
    if dry_run:
        return table
    else:
        mytracker = SampleTracker()
        mytracker.write_seq_table(table)
        mytracker.close_gumbo_client()
        print("updated gumbo")
        return None


def getQC(workspace, only=[], qcname=[], match=""):
    """
    Copied from depmap_omics.
    Will get from a workspace, the QC data for each samples

    Args:
    -----
        workspace: the workspace name
        only: do it only for this set of samples
        qcname: col name where the QC is in the workspace samples
        match: for example'.Log.final.out' get only that QC if you have a list of QCs in you qcname col

    Returns:
    --------
        a dict(sample_id:list[QC_filepaths])
    """
    if type(qcname) is str:
        qcname = [qcname]
    res = {}
    wm = dm.WorkspaceManager(workspace)
    sam = wm.get_samples()
    if len(only) > 0:
        sam = sam[sam.index.isin(only)]
    for k, val in sam.iterrows():
        res[k] = []
        for i in val[qcname]:
            if type(i) is list:
                if match:
                    res[k].extend([e for e in i if match in e])
                else:
                    res[k].extend(i)
            else:
                res[k].append(i)
    return res


def updateTrackerRNA(
    failed,
    lowqual,
    tracker,
    samplesetname,
    refworkspace=config["rnaworkspace"],
    bamfilepaths=config["starbamcolterra"],
    newgs=config["rna_hg38_path"],
    dry_run=False,
    qcname="star_logs",
    match=".Log.final.out",
    samplesinset=[],
    starlogs={},
    todrop=[],
    billing_proj=None,
):
    """updates the sample tracker with the new rna samples and the QC metrics

    Args:
        tracker (dataframe[datatype, prioritized, arxspan_id, index, ($newname)]): the sample tracker containing necessary info to compute which duplicates to keep
        selected (list[str]): which samples were selected in the release of the analysis
        samplesetname (str): the name of the sample set or of the current analysis
        samplesinset (list[str]): list of samples in the analysis.
        lowqual (list[str]): list of samples that failed QC
        newgs (str, optional): google storage path where to move the files. Defaults to ''.
        sheetcreds (str, optional): google sheet service account file path. Defaults to SHEETCREDS.
        sheetname (str, optional): google sheet service account file path. Defaults to SHEETNAME.
        qcname (str, optional): Terra column containing QC files. Defaults to "star_logs".
        refworkspace (str, optional): if provideed will extract workspace values (bam files path, QC,...). Defaults to None.
        bamfilepaths (list, optional): Terra columns containing the bam filepath for which to change the location. Defaults to STARBAMCOLTERRA.
        todrop (list, optional): list of samples to be dropped. Defaults to []
        samplesinset (list[str], optional): list of samples in set when refworkspace is None (bypass interacting with terra)
        starlogs (dict(str:list[str]), optional): dict of samples' star qc log locations when refworkspace is None (bypass interacting with terra)
    """
    refwm = dm.WorkspaceManager(refworkspace)
    if samplesinset == []:
        samplesinset = [
            i["entityName"]
            for i in refwm.get_entities("sample_set").loc[samplesetname].samples
        ]
    if starlogs == {}:
        starlogs = getQC(
            workspace=refworkspace, only=samplesinset, qcname=qcname, match=match
        )
    for k, v in starlogs.items():
        if k == "nan":
            continue
        a = tracker.loc[k, "processing_qc"]
        a = "" if a is None else a
        tracker.loc[k, "processing_qc"] = str(v) + "," + a
        if tracker.loc[k, "bam_qc"] != v[0]:
            tracker.loc[k, "bam_qc"] = v[0]
    return update(
        tracker,
        samplesetname,
        failed,
        lowqual,
        newgs=newgs,
        refworkspace=refworkspace,
        bamfilepaths=bamfilepaths,
        dry_run=dry_run,
        todrop=todrop,
        billing_proj=billing_proj,
    )


def updateTrackerWGS(
    tracker,
    samplesetname,
    lowqual,
    datatype,
    newgs=config["wgs_hg38_bam_path"],
    samplesinset=[],
    procqc=[],
    bamqc=[],
    refworkspace=None,
    bamfilepaths=["internal_bam_filepath", "internal_bai_filepath"],
    dry_run=False,
    billing_proj=None,
):
    """updates the sample tracker with the new wgs samples and the QC metrics

    Args:
        tracker (dataframe[datatype, prioritized, arxspan_id, index, ($newname)]): the sample tracker containing necessary info to compute which duplicates to keep
        selected (list[str]): which samples were selected in the release of the analysis
        samplesetname (str): the name of the sample set or of the current analysis
        samplesinset (list[str]): list of samples in the analysis.
        lowqual (list[str]): list of samples that failed QC
        newgs (str, optional): google storage path where to move the files. Defaults to ''.
        sheetcreds (str, optional): google sheet service account file path. Defaults to SHEETCREDS.
        sheetname (str, optional): google sheet service account file path. Defaults to SHEETNAME.
        procqc (list, optional): list of Terra columns containing QC files. Defaults to [].
        bamqc (list, optional): list of Terra columns containing bam QC files. Defaults to [].
        refworkspace (str, optional): if provideed will extract workspace values (bam files path, QC,...). Defaults to None.
        bamfilepaths (list, optional): Terra columns containing the bam filepath for which to change the location. Defaults to ['internal_bam_filepath', 'internal_bai_filepath'].
    """
    # computing QC
    print("looking for QC..")
    if refworkspace is not None:
        if not samplesinset:
            samplesinset = [
                i["entityName"]
                for i in dm.WorkspaceManager(refworkspace)
                .get_entities("sample_set")
                .loc[samplesetname]
                .samples
            ]
        dataProc = (
            {}
            if procqc == []
            else getQC(workspace=refworkspace, only=samplesinset, qcname=procqc)
        )
        dataBam = (
            {}
            if bamqc == []
            else getQC(workspace=refworkspace, only=samplesinset, qcname=bamqc)
        )
        for k, v in dataProc.items():
            if k == "nan":
                continue
            a = tracker.loc[k, "processing_qc"]
            a = "" if a is None else a
            tracker.loc[k, "processing_qc"] = str(v) + "," + a
        for k, v in dataBam.items():
            if k == "nan":
                continue
            a = tracker.loc[k, "bam_qc"]
            a = "" if a is None else a
            tracker.loc[k, "bam_qc"] = str(v) + "," + a
    if type(datatype) is str:
        datatype = [datatype]
    return update(
        tracker,
        samplesetname,
        lowqual,
        lowqual,
        newgs=newgs,
        refworkspace=refworkspace,
        bamfilepaths=bamfilepaths,
        dry_run=dry_run,
        samplesinset=samplesinset,
        billing_proj=billing_proj,
    )
