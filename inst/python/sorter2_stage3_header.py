def build_hybrid_header(
    splithits2, own_voucher, own_species, phase_tag, tag_specimen=False
):
    """Build the polyploid-cluster FASTA header for a phased hybrid
    haplotype's alt-lineage haplotype.

    `splithits2` is a usearch_global hit target ID already split on
    '_' (e.g. `L201_cl2_Nitta4483_ogracilis_ph0`.split('_')):
    splithits2[0]/[1] are the locus/cluster prefix, splithits2[2] is
    the matched diploid specimen (voucher), and splithits2[3] is that
    specimen's species code.

    By default (`tag_specimen=False`), the matched specimen
    (splithits2[2]) is dropped, matching existing behavior. When
    `tag_specimen` is True, it is appended as an extra field just
    before `phase_tag`, so a mapfile can target a specific voucher on
    the `progenitor` side (see joelnitta/sorter2r#16), the way #14 lets
    it target a specific voucher on the `hybrid` side.
    """
    parts = [
        splithits2[0], splithits2[1], own_voucher, own_species,
        splithits2[3],
    ]
    if tag_specimen:
        parts.append(splithits2[2])
    parts.append(phase_tag)
    return '>' + '_'.join(parts)
