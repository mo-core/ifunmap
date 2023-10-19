#!/usr/bin/env python
"""Provide a command line tool to validate configuration files."""
import argparse
import yaml
import tarfile
import logging
import sys
from pathlib import Path

logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='validating config files',
        epilog='Example: python check_config.py -c config.yml',
    )
    parser.add_argument(
        '-c',
        '--config_file',
        required=True,
        type=Path,
        help="config yaml file",
    )
    parser.add_argument(
        '-d',
        '--data_file',
        required=True,
        type=Path,
        help="tar gzipped data matrix files",
    )
    parser.add_argument(
        '-l',
        '--log-level',
        help='The desired log level (default WARNING).',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG'),
        default='WARNING',
    )
    return parser.parse_args(argv)


def is_tar_gz_file(filename):
    try:
        with tarfile.open(filename, "r:gz") as tar:
            return True
    except tarfile.ReadError:
        print('File is not a tar.gz file!')
        return False


def check_inputs(config_file, data_file):
    required_fields = ['data_files', 'data_path']
    try:
        with open(config_file, 'r') as file:
            yaml_obj = yaml.safe_load(file)
            if not all(field in yaml_obj for field in required_fields):
                print(f'configuration file required fields: {required_fields}')
                raise ValueError('Missing required fields in config file!')
            # get the value of 'data_path' field
            data_path = yaml_obj['data_path']
    except yaml.YAMLError:
        return False

    is_tar_gz_file(data_file)

    with open('data_path.txt', 'w') as file:
        file.write(data_path)

    # funmap should check if the files listed under data_files are also in the tar.gz file

    # # check all file listed under data_files are also in the tar.gz file
    # data_files = yaml_obj['data_files']
    # # list all the files in the tar.gz file
    # with tarfile.open(data_file, "r:gz") as tar:
    #     tar_files = tar.getnames()
    #     # get the file names only without the path prefix
    #     tar_files = [Path(file).name for file in tar_files]

    # # check if all files in data_files are in tar_files
    # if not all(file['path'] in tar_files for file in data_files):
    #     print(f'Files listed under data_files are not in the tar.gz file!')
    #     raise ValueError('Files listed under data_files are not in the tar.gz file!')

    with open('config.yml', 'w') as file:
        yaml.dump(yaml_obj, file)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.config_file.is_file():
        logger.error(f'The given input file {args.config_file} was not found!')
        sys.exit(2)
    check_inputs(args.config_file, args.data_file)


if __name__ == '__main__':
    sys.exit(main())
