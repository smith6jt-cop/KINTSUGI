from typing import Dict, List

EXTENSION_TYPES: Dict[str, List[str]] = {
    "IMAGE": ["tiff", "tif", "png", "jpg", "jpeg", "ome.tiff"],
    "ARCHIVE": ["tar", "gz", "zip"],
    "DATA": ["csv", "feather", "bin", "json"],
}
