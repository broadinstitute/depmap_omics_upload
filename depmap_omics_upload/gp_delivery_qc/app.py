import logging

import dalmatian
import functions_framework
import numpy as np
from cloudevents.http import CloudEvent

from depmap_omics_upload.gp_delivery_qc.gcp_utils import get_objects_metadata
from depmap_omics_upload.gp_delivery_qc.gumbo_client_wrapper import GumboClientWrapper
from depmap_omics_upload.gp_delivery_qc.pydantic_models import GPDeliveryQCConfig

@functions_framework.cloud_event
def entrypoint(cloud_event: CloudEvent) -> str:
    config = GPDeliveryQCConfig(**cloud_event.data)

    logging.info(
        f"Checking {config.data_type}"
        f"samples in {config.workspace.name}"
        f"against Gumbo {config.gumbo_client_config.env} environment"
    )

    gcw = GumboClientWrapper(config.gumbo_client_config)
    seq_table = gcw.read_seq_table()
    pr_table = gcw.read_pr_table()

    wm = dalmatian.WorkspaceManager(config.workspace.name).disable_hound()
    samples_in_ws = wm.get_samples().replace(np.nan, "", regex=True).reset_index()

    objects_metadata = get_objects_metadata(samples_in_ws[config.bam_col_name])

    return "OK"
