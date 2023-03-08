import os
import re
import numpy as np
import signal
from tqdm import tqdm


def getBamDate(bams, split="-", order="des", unknown="U"):
    """
    from bam files (could be in a google bucket) returns their likely sequencing date if available in the header

    Args:
    -----
      bams: the bams file|bucket paths
      split: the splitter in the output date
      unknown: maybe the some dates can't be found the program will output unknown for them
      order: if 'asc', do d,m,y else do y,m,d

    Returns:
    -------
      a list of likely dates or [unknown]s
    """
    DTs = []
    for i, bam in enumerate(tqdm(bams)):
        data = os.popen(
            "export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`\
    && samtools view -H "
            + bam
            + ' | grep "^@RG"'
        )
        if data == signal.SIGINT:
            print("Awakened")
            break
        else:
            res = data.read()
            dt = re.findall("(?<=\tDT:).+?\t", res)
        if len(dt) > 1:
            arr = np.array(dt[0].split("T")[0].split(split)).astype(int)
            for val in dt[1:]:
                arr = np.vstack(
                    (arr, np.array(val.split("T")[0].split(split)).astype(int))
                )
            arr = arr.T
            i = (
                arr[0] * 365 + arr[1] * 31 + arr[2]
                if order == "asc"
                else arr[2] * 365 + arr[1] * 31 + arr[0]
            )
            DTs.append(dt[np.argsort(i)[0]].split("T")[0])
        elif len(dt) == 1:
            DTs.append(dt[0].split("T")[0])
        else:
            DTs.append(unknown)
    return DTs
