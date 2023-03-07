# -*- coding: utf-8 -*-
# Jérémie Kalfon
# for BroadInsitute
# in 2019

import dalmatian as dm
from depmap_omics_upload import tracker as track


def addSamplesToDepMapWorkspace(
    stype,
    refworkspace,
    samplesetname="",
    add_to_samplesets=[],
    model_cols_to_add=["PatientID", "ModelID", "StrippedCellLineName"],
):
    """update the samples on a depmapomics terra processing workspace

    Args:
        stype (str): data type
        refworkspace (str): terra processing workspace to import data to
        add_to_samplesets (list, optional): add new samples to additional sample_sets on terra. Defaults to []
    """
    # TODO: make this work for other gumbo tables (improrts to DEPMAP_OMICS)
    mytracker = track.SampleTracker()
    refwm = dm.WorkspaceManager(refworkspace).disable_hound()

    terra_samples = refwm.get_samples()
    seq_table = mytracker.add_model_cols_to_seqtable(cols=model_cols_to_add)

    # terra requires a participant_id column
    seq_table = seq_table.rename(
        columns={
            "PatientID": "participant_id",
            "ModelID": "arxspan_id",
            "StrippedCellLineName": "stripped_cell_line_name",
        }
    )
    seq_table["participant_id"] = seq_table["participant_id"].fillna("nan")

    # check which lines are new and need to be imported to terra
    samples_to_add = seq_table[
        (~seq_table.index.isin(terra_samples.index))
        & (seq_table.expected_type == stype)
    ]
    print("found " + str(len(samples_to_add)) + " new samples to import!")

    # uploading new samples
    samples_to_add.index.name = "sample_id"

    refwm.upload_samples(samples_to_add)

    refwm.update_sample_set(
        sample_set_id="all",
        sample_ids=[i for i in refwm.get_samples().index.tolist() if i != "nan"],
    )

    # creating a sample set, if samplesetname is provided
    if samplesetname != "":
        print("adding new samples to sampleset: " + samplesetname)
        refwm.update_sample_set(
            sample_set_id=samplesetname, sample_ids=samples_to_add.index
        )

    # add new samples to additional existing sample_sets if needed
    for sname in add_to_samplesets:
        samples_in_sname = refwm.get_sample_sets().loc[sname, "samples"]
        new_samples = samples_to_add.index.tolist()

        refwm.update_sample_set(
            sample_set_id=sname, sample_ids=list(set(samples_in_sname + new_samples))
        )
