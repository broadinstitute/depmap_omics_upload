from depmap_omics_upload.mgenepy.google import gcp
from depmap_omics_upload.mgenepy.utils import helper as h
import subprocess
import pandas as pd
import dalmatian as dm


def changeToBucket(
    samples,
    gsfolderto,
    billing_proj=None,
    name_col=None,
    values=["bam", "bai"],
    filetypes=None,
    catchdup=False,
    dryrun=True,
):
    """
    moves all bam/bai files in a sampleList from Terra, to another gs bucket and rename them in the sample list

    will prevent erasing a duplicate sample by adding a random string or by flagging them and not copying them

    Args:
    ----
      samples: pandas.dataframe with columns to move
      gsfolderto: the bucket path to move the data to
      values: list of the cols in the dataframe containing the gs object path to be moved
      filetypes: list[str] of size values for each columns, give a suffix (.txt, .bam, ...)
      catchdup: if false will prepend a random string to the names before moving them, else will flag duplicate names
      dryrun: only shows the output but does not move the files

    Returns:
    --------
      the updated sample pandas.dataframe
    """
    # to do the download to the new dataspace
    cmds = []
    for i, val in samples.iterrows():
        ran = h.randomString(6, "underscore", withdigits=False)
        for j, ntype in enumerate(values):
            # TODO try:catch
            filetype = (
                ".".join(val[ntype].split("/")[-1].split(".")[1:])
                if filetypes is None
                else filetypes[j]
            )
            if name_col is None:
                name = val[ntype].split("/")[-1].split(".")[0]
            elif name_col == "index":
                name = val.name
            else:
                name = val[name_col]
            name = (
                name + "." + filetype if catchdup else name + "_" + ran + "." + filetype
            )
            if (
                not gcp.exists(gsfolderto + name, billing_proj=billing_proj)
                or not catchdup
            ):
                cmd = "gsutil cp " + val[ntype] + " " + gsfolderto + name
                if billing_proj is not None:
                    cmd = (
                        "gsutil -u "
                        + billing_proj
                        + " cp "
                        + val[ntype]
                        + " "
                        + gsfolderto
                        + name
                    )
                if dryrun:
                    print(cmd)
                    cmds.append(cmd)
                else:
                    cmds.append(cmd)
                    res = subprocess.run(cmd, shell=True, capture_output=True)
                    if res.returncode != 0:
                        raise ValueError(str(res.stderr))
                    samples.loc[i, ntype] = gsfolderto + name
            else:
                print(name + " already exists in the folder: " + gsfolderto)
                print(gcp.lsFiles([gsfolderto + name], "-la"))
                samples.loc[i, ntype] = gsfolderto + name
    return samples, cmds


def changeGSlocation(
    workspacefrom,
    newgs,
    workspaceto=None,
    prevgslist=[],
    index_func=None,
    flag_non_matching=False,
    onlysamples=[],
    bamfilepaths=[],
    entity="samples",
    droplists=True,
    keeppath=True,
    dry_run=True,
    par=20,
):
    """
  Function to move data around from one workspace to a bucket or to another workspace.

  can also work on dataframes containing lists of paths

  Args:
  -----
    workspacefrom: the workspace name where the data is
    newgs: the newgs bucket where to copy the data in
    workspaceto: if we should have these new samples and columns added to another workspace instead \
    of just updating the same one (usefull to copy one workspace to another)
    prevgslist: if providded, will only move files that are in the set of google bucket listed here
    index_func: *WIP* unused
    flag_non_matching: if set to true and prevgslist is set to some value, will return a list of samples that were not
    matched to anything in the prevgslist
    bamfilepaths: do this only on a subset of columns in terra workspace
    entity: the entity in the terra workspace on which to do this
    droplists: if set to true remove all columns containing list of paths (list of path are not uploaded well in terra)
    keeppath: if set to true, will keep the full object path and just change the bucket
    dry_run: if set to true will not update anything on Terra but just return the result
    par: on how many processor do the gs copy commands.

  Returns:
  -------
    torename: the pandas.df containing the new paths
    flaglist: the samples that were non matching (if flag_non_matching is set to true)
  """
    flaglist = []
    wmfrom = dm.WorkspaceManager(workspacefrom)
    a = wmfrom.get_entities(entity)
    if len(onlysamples) > 0:
        a = a[a.index.isin(onlysamples)]
    print("using the data from " + workspacefrom + " " + entity + " list")
    if len(a) == 0:
        raise ValueError("no " + entity)
    if bamfilepaths:
        a = a[bamfilepaths]
    todrop = set()
    torename = {}
    print(
        'this should only contains gs:// paths otherwise precise columns using "bamfilepaths"'
    )
    cmd = []
    for col in a.columns.tolist():
        val = []
        for k, prev in a[col].iteritems():
            if type(prev) is str:
                new = prev
                if newgs not in new:
                    if len(prevgslist) > 0:
                        for prevgs in prevgslist:
                            new = new.replace(prevgs, newgs)
                        if flag_non_matching:
                            if new == prev:
                                flaglist.append(prev)
                    if not keeppath:
                        new = newgs + new.split("/")[-1]
                    else:
                        new = newgs + "/".join(new.split("/")[3:])
                else:
                    print("sample " + str(k) + " was already in the new gs")
                val.append(new)
            # IN CASE WE HAVE A LIST
            if type(prev) is list:
                if droplists:
                    todrop.add(k)
                    continue
                ind = []
                for prevname in prev:
                    newname = prevname
                    if newgs not in newname:
                        if len(prevgslist) > 0:
                            for prevgs in prevgslist:
                                new = new.replace(prevgs, newgs)
                            if flag_non_matching:
                                if new == prev:
                                    flaglist.append(prev)
                        if not keeppath:
                            new = newgs + new.split("/")[-1]
                        else:
                            new = newgs + "/".join(new.split("/")[3:])
                    else:
                        print("sample " + str(k) + " was already in the new gs")
                    ind.append(newname)
                val.append(ind)
        torename.update({col: val})
        if not dry_run:
            if keeppath:
                h.parrun(
                    [
                        "gsutil mv " + a.iloc[i][col] + " " + v
                        for i, v in enumerate(val)
                    ],
                    cores=20,
                )
            else:
                gcp.mvFiles(a[col].tolist(), newgs)
        else:
            if keeppath:
                cmd = [
                    "gsutil mv " + a.iloc[i][col] + " " + v for i, v in enumerate(val)
                ]
            else:
                cmd = "mv " + str(a[col].tolist()) + " " + newgs
    torename = pd.DataFrame(
        data=torename, index=[i for i in a.index.tolist() if i != "nan"]
    )
    if workspaceto is not None:
        wmto = dm.WorkspaceManager(workspaceto)
        if not dry_run:
            wmto.disable_hound().update_entity_attributes(entity, torename)
    return torename, flaglist, cmd

async def deleteHeavyFiles(workspaceid, unusedOnly=True, ma=100, whitelist=[]):
    """
    deletes all files above a certain size in a workspace (that are used or unused)

    Args:
    ----
      workspaceid: str the name off the workspace
      unusedOnly: bool whether to delete used files as well (files that appear in one of the sample/samplesets/pairs data tables)
    """
    print("cleaning up, threshold is ", ma, "MB")
    wm = dm.WorkspaceManager(workspaceid)
    bucket = wm.get_bucket_id()
    sizes = gcp.get_all_sizes("gs://" + bucket + "/")
    print("we got " + str(len(sizes)) + " files")
    a = list(sizes.keys())
    a.sort()
    torm = []
    tot = 0
    print("whitelisting jobs: ", whitelist)
    for i in a[::-1]:
        if i > 1000000 * ma:
            for val in sizes[i]:
                if len(whitelist) > 0:
                    if not val.startswith(
                        tuple(["gs://" + bucket + "/" + wl for wl in whitelist])
                    ):
                        tot += i
                        torm.append(val)
                else:
                    tot += i
                    torm.append(val)
    print("we might remove more than " + str(tot / 1000000000) + "GB")
    if unusedOnly:
        sam = pd.concat([wm.get_samples(), wm.get_pairs(), wm.get_sample_sets()])
        tokeep = set(
            [
                val
                for val in sam.values.ravel()
                if type(val) is str and val[:5] == "gs://"
            ]
        )
        torm = set(torm) - tokeep
    return torm


def deleteJob(workspaceid, subid, taskid, deleteCurrent=False, dryrun=True):
    """
    removes files generated by a job on Terra

    Args:
    -----
      workspaceid: str wokspace name
      subid: str the name of the job
      taskid: str the name of the task in this job
      DeleteCurrent: bool whether or not to delete files if they appear in one of the sample/samplesets/pairs data tables
      dryrun: bool just plot the commands but don't execute them
    """
    wm = dm.WorkspaceManager(workspaceid)
    bucket = wm.get_bucket_id()
    data = []
    if deleteCurrent:
        if dryrun:
            print("gsutil -m rm gs://" + bucket + "/" + subid + "/*/" + taskid + "/**")
        else:
            res = subprocess.run(
                "gsutil -m rm gs://" + bucket + "/" + subid + "/*/" + taskid + "/**",
                shell=True,
                capture_output=True,
            )
            if res.returncode != 0:
                raise ValueError(str(res.stderr))
    else:
        res = subprocess.run(
            "gsutil -m ls gs://" + bucket + "/" + subid + "/*/" + taskid + "/**",
            shell=True,
            capture_output=True,
        )
        if res.returncode != 0 or len(str(res.stdout)) < 4:
            raise ValueError(str(res.stderr))
        data += str(res.stdout)[2:-1].split("\\n")[:-1]
        if "TOTAL:" in data[-1]:
            data = data[:-1]
        sam = pd.concat([wm.get_samples(), wm.get_pairs(), wm.get_sample_sets()])
        tokeep = set(
            [
                val
                for val in sam.values.ravel()
                if type(val) is str and val[:5] == "gs://"
            ]
        )
        torm = set(data) - tokeep
        if dryrun:
            print(torm)
        else:
            h.parrun(["gsutil rm " + i for i in torm], cores=12)


def removeFromFailedWorkflows(
    workspaceid, maxtime="2020-06-10", everythingFor=[], dryrun=False
):
    """
    Lists all files from all jobs that have failed and deletes them.

    Can be very long

    Args:
    -----
      workspaceid: str the workspace name
      maxtime: str date format (eg. 2020-06-10) does not delete files generated past this date
      everythingFor: list[str] removes from these workflows even if not failed
      dryrun: bool whether or not to execute or just print commands
    """
    wm = dm.WorkspaceManager(workspaceid)
    for k, val in wm.get_submission_status(filter_active=False).iterrows():
        if (
            val.Failed > 0 or val.configuration in everythingFor
        ) and val.date.date() > pd.to_datetime(maxtime):
            for w in wm.get_submission(val.submission_id)["workflows"]:
                if w["status"] == "Failed" or val.configuration in everythingFor:
                    try:
                        a = w["workflowId"]
                    # else it was not even run
                    except:
                        continue
                    deleteJob(workspaceid, val.submission_id, a, dryrun=dryrun)