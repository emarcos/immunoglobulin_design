# immunoglobulin_design

Here we describe design protocols for de novo immunoglobulin-like domains using the Rosetta modeling software as described in our manuscript "De novo design of Immunoglobulin-like domains"

https://www.biorxiv.org/content/10.1101/2021.12.20.472081v1

**Contents:**
* **scripts**
  
  The build_blueprints python script automatically generates all combinations of blueprint files exploring different strand and loop lengths, as well as specific loop ABEGO types for beta-arch loops and helices. This scripts uses the Blueprint class to manage blueprints as defined in Blueprint.py.
  
  For each of the generated blueprints the xml RosettaScripts file builds backbones by fragment insertion and performs sequence design calculations for optimizing the total score and other filters related to loop geometry, and core packing efficiency.
