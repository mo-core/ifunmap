process CHECK_CYTOSCAPE {
    executor 'local'

    script:
    """
    if docker ps --format '{{.Image}}' | grep -q "cytoscape/cytoscape-desktop"; then
       echo "A container based on cytoscape/cytoscape-desktop is running."
       exit 0
    else
       echo "No container based on cytoscape/cytoscape-desktop is running."
       echo "Please start a container based on cytoscape/cytoscape-desktop before running this workflow."
       echo 'docker run  -p 1234:1234 -v /data:/data -it cytoscape/cytoscape-desktop:3.9.1'
       exit 1
    fi
    """
}
