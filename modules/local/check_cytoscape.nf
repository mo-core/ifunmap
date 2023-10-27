process CHECK_CYTOSCAPE {
    executor 'local'

    output:
    path 'exit_code.txt', emit: exitcode


    script:
    """
    if docker ps --format '{{.Image}}' | grep -q "cytoscape/cytoscape-desktop"; then
       echo "A container based on cytoscape/cytoscape-desktop is running."
       echo "\$?" > exit_code.txt
    else
       echo "No container based on cytoscape/cytoscape-desktop is running."
       echo "Please start a container based on cytoscape/cytoscape-desktop before running this workflow."
       echo 'docker run  -p 1234:1234 -v /data:/data -it cytoscape/cytoscape-desktop:3.9.1'
       echo "\$?" > exit_code.txt
       exit 1
    fi
    """
}
