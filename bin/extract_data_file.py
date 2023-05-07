#!/usr/bin/env python
import argparse
import yaml
import tarfile
import os


def arg_parse():
    parser = argparse.ArgumentParser(description='command line arguments.')
    parser.add_argument('-c', '--config-file', type=str, required=True,
                        help='path to the config file.')
    parser.add_argument('-s', '--dataset-name', type=str, required=True,
                        help='the name of dataset that will be used.')
    parser.add_argument('-d', '--data-file', type=str, required=True,
                        help='path to the gzipped tar data file.')
    parser.add_argument('-o', '--output-file', type=str, required=True,
                        help='path to the output data file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = arg_parse()
    config = yaml.load(open(args.config_file, 'r'), Loader=yaml.FullLoader)
    all_data_set = config['data_files']
    data_set = [x for x in all_data_set if x['name'] == args.dataset_name][0]
    target_data_file = data_set['path']

    tmp_dir = 'tmp_data'
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    os.system(f'tar -xzf {args.data_file} --strip-components=1 -C {tmp_dir}')
    os.system(f'ln -s {tmp_dir}/{target_data_file} {args.output_file}')
