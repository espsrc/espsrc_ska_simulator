# CONTEXT

## Purpose

This file defines the shared domain language for extending `skasim` beyond catalog/component sky models into image- and cube-based model injection. It is intentionally conceptual. It explains the problem space, the entities involved, and how terms are used in this project.

---

## Project scope

`skasim` is a simulation workflow built around Karabo/OSKAR for generating radio-interferometric datasets, derived images, and a run report/web log. It already supports sky models expressed as catalogs or component lists. The next capability is to accept additional model products that are naturally image-like rather than catalog-like.

The goal is not to generate astrophysical models inside the workflow. The workflow accepts prepared model products and injects them into simulated observations.

---

## Core concepts

### Simulation workflow
A simulation workflow is the end-to-end process that turns a validated configuration and one or more model inputs into a Measurement Set and derived outputs such as dirty/clean images and a weblog.

### Sky model
A sky model is any description of the true or assumed sky brightness distribution used as simulation input. In this project, a sky model can take more than one form.

### Component sky model
A component sky model is a catalog-like description of the sky in terms of discrete entries such as point sources or Gaussian components. This is the existing input mode already supported by Karabo.

### Image model
An image model is a sky description expressed on a spatial pixel grid rather than as a source list. Examples include continuum intensity maps, polarization maps, and spectral cubes.

### Model entry point
A model entry point is an accepted top-level model type that the workflow knows how to validate and inject. Model preparation can happen elsewhere, but the workflow needs stable, well-defined entry points.

### Model preparation
Model preparation is the conversion of external products into one of the accepted workflow entry points. It is outside the conceptual core of `skasim`, even if helper scripts are provided for convenience.

### Injection
Injection is the act of turning a model input into simulated visibility data and merging it into the simulated observation.

### Injection backend
An injection backend is the software mechanism used to predict visibilities from a model product. Different model forms may be injected through different backends, but the workflow should expose a stable user-facing model interface.

### Composite simulation
A composite simulation is a run that combines more than one model contribution. Examples include a wide-field background source catalog plus a bright galaxy image model, or a continuum source plus a spectral-line cube.

---

## Model families

### Continuum intensity map
A continuum intensity map is a 2D spatial image of Stokes I at a reference frequency.

### Spectral index map
A spectral index map is a 2D image that gives the power-law spectral index per pixel. Together, a continuum intensity map and a spectral index map define a spatially varying continuum model.

### Continuum I+alpha model
A continuum `I+alpha` model is the pair consisting of:
- a Stokes I map at a reference frequency `ν0`, and
- a spectral index map `alpha(x,y)`.

This is the preferred continuum image entry point for version 1.

### Polarization map set
A polarization map set is a set of static spatial maps for Stokes I, Q, U, and/or V on a common grid. In version 1, these are spatial maps rather than spectral or Faraday-depth products.

### Spectral cube
A spectral cube is a three-dimensional model with axes corresponding to right ascension, declination, and frequency. This is the default future representation for line-emission models.

### Dynamic model
A dynamic model varies with time and frequency during a single observation. Examples include flares and dynamic spectra. These are future model classes and may require a special simulation mode.

### EoR-like signal
An EoR-like signal is a large spectral-spatial signal model used for cosmological simulation. In this context it is treated as a future image/cube-like model family.

---

## Axes, coordinates, and units

### Spatial grid
A spatial grid is the RA-Dec pixel grid on which an image model is defined.

### Spectral axis
The spectral axis is the frequency axis of a spectral model. In the current plan, frequency is the default spectral coordinate used for workflow inputs.

### Stokes axis
The Stokes axis identifies the polarization component represented by a map or cube. The project vocabulary uses the standard Stokes parameters I, Q, U, and V.

### Reference frequency
The reference frequency is the frequency at which a continuum intensity map is defined. It is required for `I+alpha` continuum models.

### Intrinsic / PB-corrected sky
The intrinsic or PB-corrected sky is the sky brightness before any primary-beam attenuation is applied. This is the default assumption for image-based inputs in the current plan.

### Model units
Model units describe the brightness units carried by an accepted workflow input. For image models in this project, the standard is Jy/pixel.

---

## Data products in the workflow

### Measurement Set (MS)
A Measurement Set is the main radio-interferometric visibility data product used by the workflow.

### Visibility model
A visibility model is the predicted interferometric response corresponding to one or more model inputs.

### Final observed data
Final observed data are the simulated visibilities after all intended model contributions and thermal noise have been combined.

### Derived imaging products
Derived imaging products are outputs such as dirty images or WSClean products created after the simulated visibilities have been produced.

### Weblog
The weblog is the static run report that summarizes configuration, generated files, and validation-relevant information for a simulation run.

---

## Scope boundaries

### In scope for version 1
- Existing component sky models.
- Continuum `I+alpha` image models.
- Static Stokes I/Q/U/V map sets.
- Composite runs that combine catalog-based and image-based contributions.
- Final observed Measurement Set output.

### Planned for later versions
- Spectral cubes.
- Polarization cubes.
- Rotation-measure / Faraday-structure products.
- Time-variable models such as flares or dynamic spectra.
- EoR-oriented signal injection modes.

### Out of scope for the core concept
- Physical generation of astrophysical models.
- Primary-beam, ionospheric, or other direction-dependent corruption modelling in the first version.
- Mosaic-specific behaviour in the first version.

---

## Guiding principles

### Stable workflow entry points
The workflow should accept a small number of clear, high-value model entry points rather than every upstream format directly.

### Separation of concerns
Model creation and model injection are separate concerns. The workflow focuses on validation, injection, and reporting.

### Composite-first design
The workflow should support runs that combine multiple model inputs rather than assuming only one sky description per run.

### Future compatibility
Version 1 should solve the immediate continuum and polarization-map needs without blocking later support for spectral cubes, dynamic models, or EoR-like inputs.
