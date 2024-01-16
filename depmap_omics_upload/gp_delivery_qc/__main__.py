import pandas as pd
from cloudevents.http import CloudEvent

from depmap_omics_upload.gp_delivery_qc.app import entrypoint

if __name__ == "__main__":
    pd.set_option("display.max_columns", None)
    pd.set_option("display.max_colwidth", 50)
    pd.set_option("display.max_info_columns", 999)
    pd.set_option("display.max_info_rows", 999)
    pd.set_option("display.max_rows", 20)
    pd.set_option("display.max_seq_items", None)
    pd.set_option("display.width", 150)
    pd.set_option("expand_frame_repr", True)
    pd.set_option("mode.chained_assignment", "raise")

    payload = {
        "workspace": {
            "source": "DEPMAP",
            "name": "broad-genomics-data/DepMap_WGS",
            "file_type": "cram",
        },
        "id_col_name": "SmId",
        "profile_id_col_name": "ProfileID",
        "data_type": "wgs",
        "bam_col_name": "cram_path",
        "col_name_map": {
            "name": "sample_alias",
            "bai": "crai_or_bai_path",
            "bam": "cram_or_bam_path",
            "cram": "cram_path",
            "crai": "crai_path",
            "ref_bam": "hg19_bam_filepath",
            "ref_type": "Datatype",
            "ref_bai": "hg19_bai_filepath",
            "ref_cram": "hg38_cram_filepath",
            "ref_crai": "hg38_crai_filepath",
            "version": "version",
            "primary_disease": "primary_disease",
            "ref_arxspan_id": "arxspan_id",
            "ref_name": "stripped_cell_line_name",
            "source": "source",
            "size": "size",
            "legacy_size": "legacy_size",
            "from_arxspan_id": "individual_alias",
            "ref_id": "sample_id",
            "PDO_id_terra": "PDO",
            "PDO_id_gumbo": "PdoId",
            "update_time": "update_time",
            "from_patient_id": "individual_alias",
            "patient_id": "participant_id",
            "ref_date": "date_sequenced",
            "hs_hs_library_size": "hs_hs_library_size",
            "hs_het_snp_sensitivity": "hs_het_snp_sensitivity",
            "hs_mean_bait_coverage": "hs_mean_bait_coverage",
            "hs_mean_target_coverage": "hs_mean_target_coverage",
            "hs_on_target_bases": "hs_on_target_bases",
            "total_reads": "total_reads",
            "release_date": "sequencing_date",
            "hash": "crc32c_hash",
            "legacy_hash": "legacy_crc32c_hash",
            "mean_depth": "mean_depth",
            "root_sample_id": "",
            "sm_id": "SmId",
            "profile_id": "ProfileID",
            "expected_type": "expected_type",
            "blacklist": "blacklist",
            "issue": "issue",
        },
        "min_file_size": 50000000000,
        "gumbo_client_config": {
            "env": "staging",
            "username": "gumbo-staging",
            "pr_table_name": "omics_profile",
            "pr_table_index": "ProfileID",
            "seq_table_name": "omics_sequencing",
            "seq_table_index": "SequencingID",
        },
    }

    cloud_event = CloudEvent(
        attributes={
            "specversion": "1.0",
            "type": "com.github.pull_request.opened",
            "source": "https://github.com/cloudevents/spec/pull",
            "id": "A234-1234-1234",
            "time": "2018-04-05T17:31:00Z",
        },
        data=payload,
    )

    entrypoint(cloud_event)
