# Model Input Contracts

This file complements the PRD by listing the accepted model-entry contracts in a compact, implementation-facing form.

## 1. `component_sky_model`

### Purpose
Existing Karabo/OSKAR catalog/component path.

### Scope
Version 1.

### Required content
- path to a supported catalog/component model file.

### Injection route
- Karabo/OSKAR.

### Notes
- Existing behaviour should be preserved.

---

## 2. `continuum_i_alpha`

### Purpose
Inject a spatially extended continuum source with spatially varying spectral index.

### Scope
Version 1.

### Required content
- `stokes_i`: FITS image in Jy/pixel.
- `alpha`: FITS image on the same grid.
- `reference_frequency_hz`.

### Expected assumptions
- Common spatial grid.
- PB-corrected / intrinsic sky.
- `alpha` is dimensionless.

### Injection route
- helper conversion to CASA-ready inputs,
- CASA `ft` prediction into `MODEL_DATA`.

### Validation requirements
- files exist,
- matching spatial dimensions,
- matching WCS,
- required metadata present,
- units valid.

---

## 3. `static_stokes_maps`

### Purpose
Inject static Stokes maps for polarization-oriented simulations.

### Scope
Version 1 schema. Backend injection follows the continuum implementation.

### Required content
- one or more of:
  - `stokes_i`
  - `stokes_q`
  - `stokes_u`
  - `stokes_v`
- all supplied maps must be FITS images in Jy/pixel.

### Expected assumptions
- Common spatial grid.
- PB-corrected / intrinsic sky.
- Static spatial polarization structure only.

### Injection route
- schema accepted in version 1,
- helper conversion and CASA `ft` prediction are planned after the continuum implementation.

### Validation requirements
- at least one Stokes map present,
- all supplied maps exist,
- all supplied maps have matching dimensions and WCS,
- units valid.

---

## Planned future contracts

### 4. `spectral_cube`
- scope: version 2+
- expected form: RA-Dec-frequency cube.

### 5. `polarization_cube`
- scope: version 2+
- expected form: cube containing spectral/polarization structure.

### 6. `rm_or_faraday_model`
- scope: version 2+
- expected form to be defined later.

### 7. `dynamic_model`
- scope: version 2+
- expected to cover time-variable and frequency-variable signals.
