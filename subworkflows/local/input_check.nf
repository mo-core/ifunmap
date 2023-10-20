include { INPUT_FILE_CHECK } from '../../modules/local/input_check'

workflow INPUT_CHECK {
    take:
    config_file
    data_file

    main:
    INPUT_FILE_CHECK (config_file, data_file)

    emit:
    config_file = INPUT_FILE_CHECK.out.config_file  // channel: [ config.yml ]
    data_path_file = INPUT_FILE_CHECK.out.data_path_file // channel: [ data_file_path ]
    versions = INPUT_FILE_CHECK.out.versions // channel: [ versions.yml ]
}
