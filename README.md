# immunoglobulin_design

Here we describe a series of design protocols for immunoglobulin-like domains de novo using the Rosetta modeling software

The build_blueprints python script automatically generates all combinations of blueprint files exploring different strand and loop lengths, as well as specific loop ABEGO types for beta-arch loops and helices. For each blueprint the xml RosettaScripts file builds backbones by fragment insertion and performs sequence design calculations for optimizing the total score and other filters related to loop geometry, and core packing efficiency.
