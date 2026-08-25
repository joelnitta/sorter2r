def match_progenitor_row(seq_id_parts, csv_data):
    """Resolve a sequence ID's progenitor tag to a clade using mapfile rows.

    A mapfile row's `hybrid` value may be either a bare species code
    (applies to every voucher with that code) or a `voucher_species`
    compound key (applies to that one voucher only). Compound matches
    take precedence over species-level matches for the same
    progenitor value. See joelnitta/sorter2r#14.

    Returns (clade, warn_species_key). warn_species_key is set only
    when the match came from a species-level row, so callers can warn
    that the match applied to every voucher sharing that species
    code.
    """
    sample_key = '_'.join(seq_id_parts[:2])
    species_key = seq_id_parts[1]
    target_progenitor = seq_id_parts[2]

    species_match = None
    for row in csv_data:
        if row['progenitor'] != target_progenitor:
            continue
        if row['hybrid'] == sample_key:
            return row['clade'], None
        if row['hybrid'] == species_key and species_match is None:
            species_match = row

    if species_match is not None:
        return species_match['clade'], species_key
    return None, None
