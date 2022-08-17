import argparse

def argument_parser():
    parser = argparse.ArgumentParser()    
    parser.add_argument('-o', '--output_dir', required=True, help="Output directory")
    parser.add_argument('-f', '--flux', required=True, help='Input flux file')
    parser.add_argument('-mut', '--mutation', required=True, help='Input mutation file')

    return parser
