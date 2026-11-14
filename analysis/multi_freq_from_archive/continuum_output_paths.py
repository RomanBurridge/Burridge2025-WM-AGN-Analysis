"""Input-bound model paths; do not inherit another source's loop variables."""
from __future__ import annotations

from pathlib import PurePosixPath

from continuum_measurements import portable_release_product
from selected_clean_manifest import normalize_source_name


def fit_output_directory(telescope: str, source: str, image_path: str) -> str:
    telescope = str(telescope).upper()
    normalized_source = normalize_source_name(source)
    if (not normalized_source or normalized_source in {'.', '..'}
            or '/' in normalized_source or '\\' in normalized_source):
        raise ValueError('Invalid source for a fit output directory')
    product = portable_release_product(image_path)
    parts = PurePosixPath(product).parts
    if '..' in parts:
        raise ValueError('Parent traversal in continuum input path')
    if telescope in {'ALMA', 'VLA'}:
        if len(parts) < 3 or parts[0] != telescope:
            raise ValueError('Archive input does not belong to its telescope')
        folder_source = normalize_source_name(parts[1])
        if folder_source.casefold() != normalized_source.casefold():
            raise ValueError('Archive input does not belong to the fitted source')
        return f'{telescope}/{parts[1]}/fitfiles'
    if telescope == 'SMA':
        image_source = parts[-1].removesuffix('_test.image') if parts else ''
        if normalize_source_name(image_source).casefold() != normalized_source.casefold():
            raise ValueError('SMA image does not belong to the fitted source')
        return f'SMA/{normalized_source}/fitfiles'
    raise ValueError(f'Unsupported telescope: {telescope}')
