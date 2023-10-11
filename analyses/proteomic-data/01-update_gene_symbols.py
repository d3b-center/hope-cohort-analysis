#!/usr/bin/env python3

"""Using the HGNC Gene Name Database TSV, update any old gene names in the input TSV

A lightweight program to iterate through the input tsv file and check all requested columns
that contain gene names against the HGNC database. The program does this simply by creating
a dict from the records in the HGNC file. That dict has old gene names as the keys and any
associated new gene names are stored in the value as a list. As the program iterates through
the file it will update the columns as requested using simple key lookups on
the created dict. Outputs are written line by line to an output file that can be compressed
if named GZ.

Typical usage example:
    update_gene_symbols.py -g hgnc_complete_set.txt -f fusion-dgd.tsv.gz -o test.tsv.gz
"""

import argparse
import gzip
import sys

def get_args():
    """Parse the arguments

    Args:
        None

    Returns:
        A namespace object with the arguments for this program
    """
    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    required.add_argument(
            "-g",
            "--hgnc_tsv",
            required=True,
            help="Gene name database TSV file from HGNC: hgnc_complete_set.txt")
    optional.add_argument(
            "-p",
            "--old_symbol",
            default="prev_symbol",
            help="Column name for the old gene symbol(s) in the HGNC TSV. Default: %(default)s")
    optional.add_argument(
            "-n",
            "--new_symbol",
            default="symbol",
            help="Column name for the new gene symbol in the HGNV TSV. Default: %(default)s")
    required.add_argument(
            "-f",
            "--input_tsv",
            required=True,
            help="Input TSV file: fusion-dgd.tsv.gz")
    optional.add_argument(
            "-u",
            "--update_columns",
            nargs='+',
            default=["FusionName","Gene1A","Gene1B"],
            help="Space-separated column names from the Input TSV where to update gene names \
                    (e.g. -u foo bar blah). Default: %(default)s")
    optional.add_argument(
            "-z",
            "--fake_columns",
            nargs='+',
            help="Space-separated column names to use as header for Input TSV. Overrides script header detection (e.g. -z foo bar blah).")
    optional.add_argument(
            "--retain_records",
            action='store_true',
            help="When updating a record with a new gene name, keep the original record.")
    optional.add_argument(
            "--explode_records",
            action='store_true',
            help="Return all available updated names. Will create additional records for each additional new gene name.")
    required.add_argument(
            "-o",
            "--output_filename",
            required=True,
            help="Name for the output TSV file")
    parser._action_groups.append(optional)
    return parser.parse_args()

def hgnc_tsv_to_dict(hgnc_file, old_sym, new_sym):
    """Creates dict for gene name conversion.

    Iterates through the HGNC TSV and creates a dict where old gene symbols
    are the keys and a list of associated new gene symbols are the values.

    Args:
        hgnc_file: An unopened TSV file containing the HGNC gene name information
        old_sym: String name of the column in the TSV that contains the old symbols
        new_sym: String name of the column in the TSV that contains the new symbols

    Returns:
        sym_dict: A dict mapping keys to the old gene names and values to associated
        new gene names
    """
    sym_dict = {}
    with open(hgnc_file, 'rt', encoding="utf-8") as f:
        header = f.readline().strip().split('\t')
        old_sym_index = header.index(old_sym)
        new_sym_index = header.index(new_sym)
        for line in f:
            split_line = line.strip().split('\t')
            # TSV adds quotes around entries where there are multiple old symbols; remove them
            old_sym_info = split_line[old_sym_index].replace('"','')
            new_sym_info = split_line[new_sym_index]
            # If there's no info on old symbols, return to the for loop
            if old_sym_info == '':
                continue
            # Build the dict either by appending the existing list or creating a new one
            for old_symbol in old_sym_info.split('|'):
                if old_symbol in sym_dict:
                    sym_dict[old_symbol].append(new_sym_info)
                else:
                    sym_dict[old_symbol] = [new_sym_info]
    return sym_dict

def update_input_tsv(input_tsv, sym_dict, update_columns, out_file, retain_records, explode_records, fake_columns):
    """Reads lines of the input_tsv, updates the genes, and writes to an output

    Iterates through the input tsv file and feeds each line to the line processor.
    Takes the output of the line processor and writes that to the output file.

    Args:
        input_tsv: An unopened TXT or TXT.GZ file that contains gene information to be updated
        sym_dict: A dict with old symbol keys and corresponding value list of new symbols
        update_columns: A list of columnnames in the input_tsv that must be updated
        out_file: A string filename (files ending in GZ will be compressed) for the output
        retain_records: A boolean that indicates if original records should be kept
        explode_records: A boolean that enables multirecord outputs
        fake_columns: A list of strings representing the header fields
    Returns:
        None
    """
    # output to STDERR old -> new gene changes
    print ('{}\t{}'.format("old gene", "new gene"), file=sys.stderr)
    with (gzip.open if input_tsv.endswith("gz") else open)(input_tsv, "rt", encoding="utf-8") as f:
        with (gzip.open if out_file.endswith("gz") else open)(out_file, "wt", encoding="utf-8") as w:
            if fake_columns is None:
                header = f.readline().strip().split('\t')
                w.write('\t'.join(header) + '\n')
            else:
                header = '\t'.join(fake_columns)
            for line in f:
                if retain_records:
                    w.write(line)
                # Start with one line
                new_lines = [line.strip()]
                # Iteratively recreate and expand the new_lines list
                for column in update_columns:
                    new_lines = process_lines(new_lines, sym_dict, header, column, explode_records)
                # It's possible there are duplicate entires, get rid of them
                uniq_lines = list(set(new_lines))
                for new_line in uniq_lines:
                    # No need to print it again, we already have it
                    if retain_records and new_line == line:
                        continue
                    # STDOUT log for new entries
                    if new_line != line:
                        print(f'New entry added: {new_line.strip()}')
                    w.write(new_line)

def process_lines(input_lines, sym_dict, header, column_name, explode_records):
    """Given a block of lines, process each one with update_tsv_line

    This function is an abstraction to handle the iterative process introduced by the
    explode_records functionality. Basically a single record can become many records depending
    on how many new genes and old gene can have. Here we have a function that iterates over those
    lines and calls the update_tsv_line function to give us a new set of lines that can come
    back through here in future iterations. One key thing is to flatten the array on the way out as
    we'll be creating a nested array as we iterate.

    Args:
        input_lines: A list of line strings to iterate through and process
        sym_dict: A dict with old symbol keys and corresponding value list of new symbols
        column_name: A string of the columnname in the input_tsv that must be updated
        header: A list containing the columnnames from the input_tsv header
        explode_records: A boolean that enables multirecord outputs

    Returns:
        A list of record lines with updated gene names.
    """
    outarr = []
    for line in input_lines:
        outarr.append(update_tsv_line(line, sym_dict, header, column_name, explode_records))
    flatarr = [item + '\n' for sub_list in outarr for item in sub_list]
    return flatarr

def update_tsv_line(input_line, sym_dict, header, column_name, explode_records):
    """Take an input line and update the requested column, return the line

    Given a single line from the input tsv file, this function is tasked with invoking the
    gene name updater appropriately. In most cases, it expects only a single gene name to be present
    in a given column. In this case, it will use the gene name updater to get the new gene name. It
    will also check for fusions (a special use case) by string checking for '--'. If it does identify a fusion it will
    split those gene names and send both off to be processed individually. It will then reanneal
    them with the fusion notation ('--'). Once it has the new gene names it will iterate over them,
    update the column value, and add that line to an array. That array is then returned.

    Args:
        input_line: A string record from the input tsv file. Represents a single line
        sym_dict: A dict with old symbol keys and corresponding value list of new symbols
        column_name: A string of the columnname in the input_tsv that must be updated
        header: A list containing the columnnames from the input_tsv header
        explode_records: A boolean that enables multirecord outputs

    Returns:
        List of all potential newly updated lines.
    """
    split_fuse = input_line.strip().split('\t')
    outlines = []
    col_index = header.index(column_name)
    if '--' in split_fuse[col_index]:
        genea, geneb = split_fuse[col_index].split('--')
        new_geneas = update_gene_name(genea, sym_dict, explode_records)
        new_genebs = update_gene_name(geneb, sym_dict, explode_records)
        new_entries = [i + "--" + j for i in new_geneas for j in new_genebs]
    else:
        gene = split_fuse[col_index]
        new_entries = update_gene_name(gene, sym_dict, explode_records)
    for new_entry in new_entries:
        split_fuse[col_index] = new_entry
        outlines.append('\t'.join(split_fuse))
    return outlines 

def update_gene_name(old_gene, sym_dict, explode_records):
    """Given an gene name and a symbol dict, update the gene name, if possible

    Core component of the updating mechanism: Put simply, if the old gene name is a key
    in the sym_dict, return the corresponding value. WARNING: The HGNC file sometimes
    provides multiple new gene names for an old name. We don't have a solution for this
    at the moment so the program should crash if it encounters this scenario.

    Args:
        old_gene: A string name for the old_gene
        sym_dict: A dict with old symbol keys and corresponding value list of new symbols
        explode_records: A boolean that enables multirecord outputs

    Returns:
        new_genes: A list of strings with new gene names or an empty list if no upgrade.

    Raises:
        SystemExit: When we encounter an old gene name with more than one corresponding
        new gene name and explode_records is not set.
    """
    new_genes = [old_gene]
    if old_gene in sym_dict:
        # print change to stderr log
        print ('{}\t{}'.format(old_gene, ','.join(sym_dict[old_gene])), file=sys.stderr)
        new_genes = sym_dict[old_gene]
        if not explode_records and len(new_genes) > 1:
            sys.exit(f"Error updating gene {old_gene}. HGNC has multiple options: {sym_dict[old_gene]}. To output all potential new genes use --explode_records flag.")
    return new_genes

def main():
    """
    Choo choo
    """
    args = get_args()
    sym_dict = hgnc_tsv_to_dict(args.hgnc_tsv, args.old_symbol, args.new_symbol)
    update_input_tsv(args.input_tsv, sym_dict, args.update_columns, args.output_filename, args.retain_records, args.explode_records, args.fake_columns)

if __name__ == "__main__":
    main()
