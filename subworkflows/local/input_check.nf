include { INPUT_FILE_CHECK } from '../../modules/local/input_check'

workflow INPUT_CHECK {
    take:
    config_file
    data_file

    main:
    INPUT_FILE_CHECK (config_file, data_file)
        .config_file
        .set { config_file }

    emit:
    config_file // channel: [ config.yml ]
    versions = INPUT_FILE_CHECK.out.versions // channel: [ versions.yml ]
}
