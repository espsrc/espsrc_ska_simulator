# skasim

skasim creates synthetic radio-interferometric observations from sky descriptions and imaging choices. This glossary defines the project language used by the simulator, CLI, and documentation.

## Language

**Catalogue**:
A named or file-backed collection of sky sources used to build a sky model. Built-in catalogues are selected by name, such as `MIGHTEE`, `GLEAM`, or `SKAMid`; `catalogue` is the canonical term and `catalog` may appear only as a CLI spelling alias.
_Avoid_: Catalog, source database, numeric catalogue ID

**Sky model**:
The collection of sky sources that defines what the observation will simulate. It may come from a file, a catalogue, or generated source intensities.
_Avoid_: Source database, input catalogue

**Sky model source**:
One input used to produce a sky model, such as a file, a catalogue, or generated source intensities. Version 0.2 accepts one sky model source per run; later releases may combine multiple sky model sources into one sky model.
_Avoid_: Model, input

**Source intensity**:
The Stokes-I flux density assigned to a generated sky source. Multiple source intensities define multiple generated sources, not the I/Q/U/V polarization components of one source.
_Avoid_: Flux vector, IQUV list

**Run**:
One execution of the simulator from a resolved configuration to a manifest, outputs, and weblog. A run can complete or fail, but it always leaves inspectable records.
_Avoid_: Job, session

**Weblog**:
The always-on human-readable record of a run. It summarizes the run status, configuration, milestones, errors, and image products.
_Avoid_: Optional report, dashboard

**Image product**:
One imaging result produced from simulated visibilities with a particular imaging configuration. A run may produce one image product in 0.2 and multiple image products in later releases.
_Avoid_: The image, plot

**WSClean command**:
The executable invocation used for cleaned imaging. It defaults to `wsclean`, but may name a local wrapper or container invocation when WSClean is provided by the execution environment.
_Avoid_: WSClean path, cleaner binary

## Example Dialogue

Dev: "Should this run use a catalogue or inline random sources?"

Domain expert: "Use catalogue 1, the MIGHTEE catalogue, and keep `--catalog` as an alias for users who type the American spelling."
Correction: for 0.2, say `--catalogue MIGHTEE`; numeric catalogue IDs are not part of the domain language.

Dev: "How should cleaned imaging find WSClean on this machine?"

Domain expert: "Set the WSClean command to the Singularity invocation for this run; the simulator should still default to plain `wsclean`."

Dev: "Does `--flux-density 1 5 10` describe one polarized source?"

Domain expert: "No, it describes three generated sources with those Stokes-I source intensities."

Dev: "Can one run use a background catalogue and a foreground FITS file?"

Domain expert: "Not in 0.2; that is multiple sky model sources, which should be possible in a later release."

Dev: "Does a run always produce exactly one image?"

Domain expert: "No; 0.2 may produce one image product, but later releases should support several image products with different imaging configurations."

Dev: "Can I disable the weblog for a quick run?"

Domain expert: "No; the weblog is part of the run record and is always produced."
