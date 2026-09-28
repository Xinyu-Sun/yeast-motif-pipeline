#!/usr/bin/env python3
"""Summarize shared raw SGD domains across the 48 polymerase TF proteins.

Outputs one row per shared non-MobiDBLite raw domain key with:
- SGD description(s)
- number of proteins sharing the domain
- which proteins share it
- category counts and percentages
- dominant category and its percentage
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_HITS = REPO_ROOT / "data" / "output" / "polymerase_tf_domain_hits_raw.csv"
DEFAULT_OUT = REPO_ROOT / "data" / "output" / "all48_shared_domain_category_summary.csv"


def pct(n: int, d: int) -> str:
    return f'{(100.0 * n / d):.1f}%' if d else '0.0%'


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hits', type=Path, default=DEFAULT_HITS)
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--min-proteins', type=int, default=2)
    args = parser.parse_args()

    grouped = {}
    with args.hits.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            raw_key = (row.get('raw_domain_key') or '').strip()
            source = (row.get('source') or '').strip()
            if not raw_key or source == 'MobiDBLite' or 'MobiDBLite' in raw_key:
                continue
            entry = grouped.setdefault(
                raw_key,
                {
                    'source': source,
                    'raw_domain_id': (row.get('raw_domain_id') or '').strip(),
                    'descriptions': set(),
                    'proteins': set(),
                    'standard_names': set(),
                    'categories_by_y': {},
                },
            )
            desc = (row.get('domain_description') or '').strip()
            if desc and desc != '-':
                entry['descriptions'].add(desc)
            y_name = row['y_name']
            entry['proteins'].add(y_name)
            std = (row.get('standard_name') or '').strip()
            if std:
                entry['standard_names'].add(std)
            entry['categories_by_y'][y_name] = row['category']

    out_rows = []
    for raw_key, entry in grouped.items():
        shared_n = len(entry['proteins'])
        if shared_n < args.min_proteins:
            continue
        category_counts = Counter(entry['categories_by_y'].values())
        dominant_category, dominant_count = sorted(category_counts.items(), key=lambda kv: (-kv[1], kv[0]))[0]
        descriptions = '; '.join(sorted(entry['descriptions']))
        proteins = sorted(entry['proteins'])
        standard_names = sorted(entry['standard_names'])
        out_rows.append(
            {
                'raw_domain_key': raw_key,
                'source': entry['source'],
                'raw_domain_id': entry['raw_domain_id'],
                'domain_description': descriptions,
                'proteins_sharing_count': shared_n,
                'shared_y_names': '; '.join(proteins),
                'shared_standard_names': '; '.join(standard_names),
                'pol_i_count': category_counts.get('Pol I', 0),
                'pol_i_pct': pct(category_counts.get('Pol I', 0), shared_n),
                'pol_ii_count': category_counts.get('Pol II', 0),
                'pol_ii_pct': pct(category_counts.get('Pol II', 0), shared_n),
                'pol_iii_count': category_counts.get('Pol III', 0),
                'pol_iii_pct': pct(category_counts.get('Pol III', 0), shared_n),
                'dominant_category': dominant_category,
                'dominant_category_count': dominant_count,
                'dominant_category_pct': pct(dominant_count, shared_n),
            }
        )

    out_rows.sort(key=lambda row: (-int(row['proteins_sharing_count']), row['raw_domain_key']))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                'raw_domain_key',
                'source',
                'raw_domain_id',
                'domain_description',
                'proteins_sharing_count',
                'shared_y_names',
                'shared_standard_names',
                'pol_i_count',
                'pol_i_pct',
                'pol_ii_count',
                'pol_ii_pct',
                'pol_iii_count',
                'pol_iii_pct',
                'dominant_category',
                'dominant_category_count',
                'dominant_category_pct',
            ],
        )
        writer.writeheader()
        writer.writerows(out_rows)

    print(f'Wrote summary: {args.out}')
    print(f'Shared domains summarized: {len(out_rows)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
