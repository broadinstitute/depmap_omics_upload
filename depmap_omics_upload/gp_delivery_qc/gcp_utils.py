from collections.abc import Iterable

import pandas as pd
from google.cloud import storage


def get_objects_metadata(uris: Iterable[str]) -> pd.DataFrame:
    storage_client = storage.Client()

    blobs = {}

    with storage_client.batch():
        for uri in uris:
            blob = storage.Blob.from_string(uri, client=storage_client)
            bucket = storage_client.bucket(blob.bucket.name)
            blob = bucket.get_blob(blob.name)
            blobs[uri] = blob

    metadata = [
        {"uri": k, "crc32c": v.crc32c, "size": v.size, "updated": v.updated}
        for k, v in blobs.items()
    ]

    return pd.DataFrame(metadata)
