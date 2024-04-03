/*
General concept for the parallel framework:

1.  Domain decomposition:
    1.1 All algorithms only perform streaming and collision for fluid nodes. 
        When decompositing the domain, store the fluid nodes for each respective subdomain.
    1.2 A layer is full-width, i.e. each layer has the horizontal extend HORIZONTAL_NODES.
        The height is a fixed value. It is assumed that an even decomposition can be applied,
        which means that no subdomains have to be clamped.
    1.3 Apart from the lower-most and upper-most subdomain, all subdomains border two buffers
        which are stripes of the vertical extend 1.
2.  Pre-streaming buffer update:
    2.1 Instream algorithms: For every buffer node: Copy bundles 0 (southern) and 2 (northern) from neighbors
    2.2 Outstream algorithms:
        a.  Two-step algorithm: 
            -   Perform outstream for nodes below (active stream) and above (passive stream) buffer layer
            -   Ignore those lines while performing streams in the corresponding directions
            -   "Streaming" from buffer: copy value to node above (active) or below (passive)
        b.  Swap algorithm: 
            For every buffer node: 
            -   Copy bundles 0 (southern) and 2 (northern) from neighbors
            -   Perform swap step for bundle 2
3.  Continue with corresponding algorithm in subdomains
*/