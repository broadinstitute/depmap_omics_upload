# GCPFunction.py

import os
import subprocess
import re
from depmap_omics_upload.mgenepy.utils import helper as h


def lsFiles(files, add="", group=50, billing_proj=None):
    """
    list a set of files in parallel (when the set is huge)

    Args:
    ----
        files: gs paths
        add: additional params to add
        group: files to do in parallel
    """
    print("listing files in gs")
    by = len(files) if len(files) < group else group
    res = []
    for sfiles in h.grouped(files, by):
        a = ""
        for val in sfiles:
            a += val + " "
        data = None
        if billing_proj is None:
            data = subprocess.run(
                "gsutil -m ls " + add + " " + a,
                capture_output=True,
                shell=True,
            )
        else:
            data = subprocess.run(
                "gsutil -u " + billing_proj + " -m ls " + add + " " + a,
                capture_output=True,
                shell=True,
            )
        if data.returncode != 0:
            if "One or more URLs matched no objects" not in str(data.stderr):
                raise ValueError("issue with the command: " + str(data.stderr))
        if len(str(data.stdout)) < 4:
            return []
        res += (
            str(data.stdout)[2:-1].split("\\n")[:-1]
            if "L" not in add
            else ["gs://" + i for i in str(data.stdout).split("\\ngs://")]
        )
        if "TOTAL:" in res[-1] and "L" not in add:
            res = res[:-1]
    return res


def cpFiles(files, location, group=50):
    """
    copy a set of files in parallel (when the set is huge)

    Args:
    ----
        files: gs paths
        location to copy
        group: files to do in parallel
    """
    by = len(files) if len(files) < group else group
    for sfiles in h.grouped(files, by):
        a = ""
        for val in sfiles:
            a += val + " "
        code = os.system("gsutil -m cp " + a + location)
        if code != 0:
            print("pressed ctrl+c or command failed")
            break


def exists(val, billing_proj=None):
    """
    tells if a gcp path exists
    """
    if type(val) is str:
        if billing_proj is None:
            return os.popen("gsutil ls " + val).read().split("\n")[0] == val
        else:
            return (
                os.popen("gsutil -u " + billing_proj + " ls " + val)
                .read()
                .split("\n")[0]
                == val
            )
    elif type(val) is list:
        rest = set(val) - set(lsFiles(val))
        return len(rest) == 0, rest


def extractSize(val):
    """
    extract the size from the string returned by an ls -l|a command
    """
    return "gs://" + val.split("gs://")[1].split("#")[0], int(
        re.split("\d{4}-\d{2}-\d{2}", val)[0]
    )


def extractTime(val):
    """
    extract the size from the string returned by an ls -l|a command
    """
    return val.split("  ")[1].split("T")[0]


def extractPath(val):
    """
    extract the path from the string returned by an ls -l|a command
    """
    return "gs://" + val.split("gs://")[1].split("#")[0]


def extractHash(val, typ="crc32c"):
    """
    extract the crc32 from the string returned by an ls -L command

    Args:
    ----
        type: flag ['crc32c','md5']
    """
    if "    Hash (crc32c):" in val and typ == "crc32c":
        return (
            val.split("    Hash (crc32c):          ")[-1]
            .split("\\\\n")[0]
            .split("\\n")[0]
        )
    elif "    Hash (md5):" in val and typ == "md5":
        return (
            val.split("    Hash (md5):          ")[-1].split("\\\\n")[0].split("\\n")[0]
        )
    else:
        return None
