import re


def parse_readstats(lines, verbose=False):
    """Parse samtools flagstat output plus the two appended samtools
    depth lines into a dict of {stat_label: [int_value]}.

    The two depth lines (mean depth, then coverage percent) are
    matched by content -- a bare float, with no flagstat-style "+"
    separator -- rather than by line position, since flagstat's line
    count varies by samtools version/build. See
    joelnitta/sorter2r#13.
    """
    stats = {}
    for line in lines:
        if '+' in line:
            parts = line.split(" ")
            statlabel = parts[3].strip('\n')
            statint = parts[0].strip('\n')
            if verbose:
                print(statlabel)
                print(statint)
            stats[statlabel] = [int(statint)]
        elif re.fullmatch(r'[\d.]+\n?', line):
            statlabel = 'readdepth' if 'readdepth' not in stats else 'coverage'
            statint = line.strip('\n')
            if verbose:
                print(statlabel + ' = ' + statint)
            stats[statlabel] = [int(float(statint))]
    return stats
