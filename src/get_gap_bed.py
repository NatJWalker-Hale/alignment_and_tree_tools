#!/usr/bin/env python3

"""
Convert FASTA assembly N-regions to BED format.
Identifies stretches of 'N' characters and outputs their positions in BED format.
Claude AI Sonnet 4.5 / NWH
"""

import sys
import argparse


def fasta_n_to_bed(input_fasta, invert=False):
    """
    Process FASTA file and write N-regions (or non-N regions if inverted) to BED format on stdout.
    
    Args:
        input_fasta: Path to input FASTA file
        invert: If True, output non-gap regions instead of gap regions (default: False)
    """
    with open(input_fasta, 'r', encoding='utf-8') as fin:
        current_chrom = None
        position = 0
        region_start = None

        for line in fin:
            line = line.rstrip('\n\r')

            # Handle header lines
            if line.startswith('>'):
                # Write any pending region from previous chromosome
                if region_start is not None:
                    print(f"{current_chrom}\t{region_start}\t{position}")
                    region_start = None

                # Extract chromosome name (everything after '>' until first whitespace)
                current_chrom = line[1:].split()[0]
                position = 0
                continue

            # Process sequence line character by character
            for char in line:
                char_upper = char.upper()
                is_gap = char_upper == 'N'

                # Determine if this character matches what we're tracking
                matches_target = is_gap if not invert else not is_gap

                if matches_target:
                    # Start of target region or continuation
                    if region_start is None:
                        region_start = position
                else:
                    # Non-target character: close any open region
                    if region_start is not None:
                        print(f"{current_chrom}\t{region_start}\t{position}")
                        region_start = None

                position += 1

        # Write final region if file ends in target characters
        if region_start is not None:
            print(f"{current_chrom}\t{region_start}\t{position}")


def main():
    """Main function to parse arguments and call processing function."""
    parser = argparse.ArgumentParser(
        description='Convert N-regions in FASTA assembly to BED format'
    )
    parser.add_argument('input_fasta',help='input FASTA file')
    parser.add_argument('-i', '--invert', action='store_true',
                        help='invert the output to report non-N regions')
    args = parser.parse_args(sys.argv[1:] or ["--help"])
    try:
        fasta_n_to_bed(args.input_fasta, args.invert)
    except FileNotFoundError:
        print(f"Error: Input file '{args.input_fasta}' not found", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
