process CHECK_CYTOSCAPE {
    executor 'local'

    script:
    """
    if docker ps --format '{{.Image}}' | grep -q "cytoscape/cytoscape-desktop"; then
       echo "A container based on cytoscape/cytoscape-desktop is running."
    else
       echo "No container based on cytoscape/cytoscape-desktop is running."
       exit 1
    fi
    """
}
