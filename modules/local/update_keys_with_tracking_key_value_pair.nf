// main.nf
process UPDATE_KEYS {

    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(source_file)
    tuple val(meta), path(target_file)

    output:
    tuple val(meta), path("*.json"), emit: update_keys_with_tracking_key_value_pair

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import os
    import json

    def update_keys(source_file, target_file):
        try:
            with open(source_file, 'r') as file:
                source_data = json.load(file)

            with open(target_file, 'r') as file:
                target_data = json.load(file)

            updated_data = {}
            for key, value in target_data.items():
                if key in source_data:
                    # If the key exists in source_data, update the key
                    new_key = source_data[key]
                    # Use setdefault to append values for the same new_key
                    updated_data.setdefault(new_key, []).extend(value)

            # Generate output JSON file name based on target file
            output_json_name = os.path.basename(target_file).replace('.json', '_updated.json')
            with open(output_json_name, 'w') as file:
                json.dump(updated_data, file, indent=4)

            return "Update successful and saved to new file"
        except Exception as e:
            return f"Error: {e}"

    result = update_keys('${source_file}', '${target_file}')
    print(result)
    """

}

// workflow {
//     // Define the input file paths
//     ch_keyValueFromTracking = Channel.fromPath('/Users/dxu/Documents/compareJoelGTFwithMyGTF/updateGeneID/outputFromGffcompare/keyvalue.json')
//     ch_cmpRefDataAsKey = Channel.fromPath('/Users/dxu/Documents/compareJoelGTFwithMyGTF/updateGeneID/outputFromGffcompare/cmp_ref_data.json')

//     // Process the file
//     updated_json = UPDATE_KEYS(ch_keyValueFromTracking.map { [ [:], it ] }, ch_cmpRefDataAsKey.map { [ [:], it ] })



// }
