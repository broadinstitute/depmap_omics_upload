import sys
import logging
import json
import pandas as pd
import datetime

from mgenepy import terra
from depmap_omics_upload import tracker as track


def addSamplesToGumbo(
    gumbo_env,
    samples,
    stype,
    bucket,
    name_col="index",
    values=["hg19_bam_filepath", "hg19_bai_filepath"],
    filetypes=["bam", "bai"],
    billing_proj=None,
    dryrun=False,
):
    """update the samples on gumbo's sequencing sheet

    Args:
        samples ([type]): [description]
        stype ([type]): [description]
        bucket ([type]): [description]
        name_col (str, optional): [description]. Defaults to "index".
        values (list, optional): [description]. Defaults to ['legacy_bam_filepath', 'legacy_bai_filepath'].
        filetypes (list, optional): [description]. Defaults to ['bam', 'bai'].
    """
    # uploading to our bucket (now a new function)
    assert set(values).issubset(set(samples.columns))
    logging.info("copying files to depmap omics bucket")
    samples, cmds = terra.changeToBucket(
        samples,
        bucket,
        billing_proj=billing_proj,
        name_col=name_col,
        values=values,
        filetypes=filetypes,
        catchdup=True,
        dryrun=dryrun,
    )
    logging.info("finished copying files to depmap omics bucket")

    mytracker = track.SampleTracker(gumbo_env=gumbo_env)
    ccle_refsamples = mytracker.read_seq_table()

    names = []
    subccle_refsamples = ccle_refsamples[ccle_refsamples["expected_type"] == stype]
    for k, val in samples.iterrows():
        val = val["ProfileID"]
        if val != "":
            names.append(val)
            samples.loc[k, "version"] = len(
                subccle_refsamples[subccle_refsamples["ProfileID"] == val]
            ) + names.count(val)
    samples["version"] = samples["version"].astype(int)

    if not dryrun:
        logging.info("inserting new samples to gumbo sequencing table")
        mytracker.insert_to_seq_table(samples)
    mytracker.client.close()

    return samples, cmds


if __name__ == "__main__":
    config_filename = sys.argv[1]
    with open(config_filename, "r") as f:
        config = json.load(f)

    today = (
        str(datetime.datetime.now().replace(microsecond=0))
        .replace(" ", "_")
        .replace(":", "-")
    )

    logging.basicConfig(
        filename=config["loading_workingdir"] + today + "_add-to-gumbo.log",
        level=logging.INFO,
    )

    for fn in config["new_wgs_sample_tables"]:
        sampletable = pd.read_csv(fn, index_col=0)
        addSamplesToGumbo(
            config["gumbo_env"],
            sampletable,
            "wgs",
            config["wgs_hg38_cram_path"],
            values=["hg38_cram_filepath", "hg38_crai_filepath"],
            filetypes=["cram", "crai"],
            billing_proj=config["gcp_billing_proj"],
            dryrun=config["dryrun"],
        )

    for fn in config["new_rna_sample_tables"]:
        sampletable = pd.read_csv(fn, index_col=0)
        addSamplesToGumbo(
            config["gumbo_env"],
            sampletable,
            "rna",
            config["rna_hg19_path"],
            values=["hg19_bam_filepath", "hg19_bai_filepath"],
            filetypes=["bam", "bai"],
            billing_proj=config["gcp_billing_proj"],
            dryrun=config["dryrun"],
        )

    for fn in config["new_wes_sample_tables"]:
        sampletable = pd.read_csv(fn, index_col=0)
        addSamplesToGumbo(
            config["gumbo_env"],
            sampletable,
            "wes",
            config["wes_hg19_bam_path"],
            values=["hg19_bam_filepath", "hg19_bai_filepath"],
            filetypes=["bam", "bai"],
            billing_proj=config["gcp_billing_proj"],
            dryrun=config["dryrun"],
        )
