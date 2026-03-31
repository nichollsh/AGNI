# AGNI Second-Pass Audit Report
**Chemistry Scheme & User Interface Analysis**

**Date**: 2026-03-31  
**Focus Areas**: Chemistry configuration, parameter validation, user-facing messages, interface consistency

---

## Executive Summary

This second-pass audit identifies **21 distinct issues** across chemistry scheme implementation and user interfaces, ranging from HIGH to MINOR severity. The chemistry scheme shows functional correctness but suffers from:

1. **Confusing field naming** (gas_ovmr, gas_cvmr, gas_vmr triplet)
2. **Missing parameter validation** (FastChem tolerances, iteration bounds)
3. **Inconsistent user-facing messages** (unit notation, error formatting)
4. **Undocumented behavior** (metallicity format, composition mutual exclusivity)

**Overall Assessment**: ⚠️ **FAIR** — Code is functional but user experience could be significantly improved with better validation, clearer documentation, and consistent messaging.

---

## Table of Contents

1. [Chemistry Scheme Analysis](#1-chemistry-scheme-analysis)
2. [Configuration Parameter Audit](#2-configuration-parameter-audit)
3. [User Interface Messages](#3-user-interface-messages)
4. [Critical Issues Summary](#4-critical-issues-summary)
5. [Recommendations](#5-recommendations)

---

## 1. Chemistry Scheme Analysis

### 1.1 Field Naming Inconsistencies

#### **ISSUE 1: Gas VMR Triplet Confusion** ⚠️ **SEVERITY: HIGH**

**Location**: `src/atmosphere.jl` lines 168-170

Three similar fields store gas mixing ratios:

```julia
gas_vmr::Dict{String, Array{Float64,1}}     # runtime calculated VMRs
gas_cvmr::Dict{String, Array{Float64,1}}    # post-chemistry, pre-rainout VMRs
gas_ovmr::Dict{String, Array{Float64,1}}    # original input VMRs
```

**Problems**:
1. **Non-intuitive naming**: The distinction between these three is only explained in `src/chemistry.jl` module docstring (lines 12-15), not in the struct definition
2. **Risk of misuse**: Functions throughout the code must carefully choose which VMR field to use
3. **No inline documentation**: A developer modifying the code might use the wrong field

**Documentation** (only in chemistry.jl docstring):
> * `gas_ovmr` stores VMRs inputted by the user, which are usually constant in height.
> * `gas_cvmr` stores the VMRs calculated by fastchem, which are used as the starting point for condensation calculations.
> * `gas_vmr` stores the runtime gas volume mixing ratios, after all calculations are performed.

**Evidence of Confusion**:
- `chemistry.jl` line 92-93: `restore_composition!()` copies `gas_ovmr` → `gas_vmr` and `gas_ovmr` → `gas_cvmr`
- `chemistry.jl` line 733-738: `reset_to_chem!()` copies `gas_cvmr` → `gas_vmr`
- `chemistry.jl` line 784-790: `calc_composition!()` conditionally calls `reset_to_chem!()` based on `coldtrap` flag

**Recommendation**:
```julia
# Rename for clarity:
gas_vmr_runtime::Dict{String, Array{Float64,1}}     # Current runtime VMRs (post-all-physics)
gas_vmr_chemistry::Dict{String, Array{Float64,1}}   # VMRs after chemistry, before rainout
gas_vmr_initial::Dict{String, Array{Float64,1}}     # User-specified initial VMRs

# OR add extensive inline documentation:
gas_vmr::Dict{String, Array{Float64,1}}             # RUNTIME: Final VMRs after chemistry+rainout+oceans
gas_cvmr::Dict{String, Array{Float64,1}}            # CHECKPOINT: VMRs after FastChem, before rainout
gas_ovmr::Dict{String, Array{Float64,1}}            # INPUT: User-specified initial VMRs (constant)
```

---

#### **ISSUE 2: Metallicity Units Ambiguity** ⚠️ **SEVERITY: HIGH**

**Location**: `src/AGNI.jl` line 275, `src/atmosphere.jl` lines 420, 848-850

**Problem**: Configuration file comment and code documentation disagree on metallicity units.

**In TOML** (`res/config/struct_grid.toml` line 29):
```toml
metallicities = {"C"=1.0, "O"=10.0, "N"=0.3, "S"=0.1}  # temp mass fractions
```
Comment says **"mass fractions"** but code treats them as **molar ratios relative to hydrogen**!

**In Code** (`atmosphere.jl` line 850):
```julia
atmos.metal_orig[k] = metallicities[k] * phys._get_mmw("H") / phys._get_mmw(k)
```
This converts from **molar ratio (relative to H)** to **molar number** by applying molecular weight correction.

**FastChem Input Format** (`chemistry.jl` lines 530-542):
```julia
write(f, @sprintf("%s    %.3f \n", e, log10(atmos.metal_calc[e]) + 12.0))
```
This writes `log10(N_e/N_H) + 12`, confirming **molar number density relative to hydrogen**.

**Documentation** (`atmosphere.jl` line 420):
```julia
# - `metallicities::Dict`   dictionary of elemental metallicities (mass ratio rel to hydrogen)
```
This says **"mass ratio"** but code implementation treats as **molar ratio**!

**Recommendation**:
1. **Fix TOML comment** in `struct_grid.toml`:
   ```toml
   metallicities = {"C"=1.0, "O"=10.0, "N"=0.3, "S"=0.1}  # molar ratios relative to H
   ```
2. **Fix docstring** in `atmosphere.jl` line 420:
   ```julia
   # - `metallicities::Dict`   dictionary of elemental abundances (molar number ratio rel. to hydrogen)
   ```
3. **Add validation** to check reasonable ranges (e.g., 0.01 to 100 for most elements)

---

### 1.2 Missing Parameter Validation

#### **ISSUE 3: FastChem Parameters Lack Bounds Checking** ⚠️ **SEVERITY: MEDIUM**

**Location**: `src/atmosphere.jl` lines 512-517 (defaults), no validation

**Unvalidated Parameters**:

| Parameter | Default | Should Validate | Current Status |
|-----------|---------|-----------------|----------------|
| `fastchem_floor` | 400.0 K | > 0, < 10000 K | ❌ No check |
| `fastchem_maxiter_chem` | 80000 | > 0, < 1e7 | ❌ No check |
| `fastchem_maxiter_solv` | 40000 | > 0, < 1e7 | ❌ No check |
| `fastchem_xtol_chem` | 1.0e-3 | (0, 1) | ❌ No check |
| `fastchem_xtol_elem` | 1.0e-3 | (0, 1) | ❌ No check |

**Risk**: Users could accidentally set:
- `fastchem_floor = -100` (negative temperature)
- `fastchem_xtol_chem = 2.0` (tolerance > 1, meaningless)
- `fastchem_maxiter_chem = 0` (no iterations allowed)

**Recommendation**:
```julia
# Add to setup!() after line 517:
_check_range("FastChem temperature floor", fastchem_floor; min=100.0, max=10000.0) || return false
_check_range("FastChem chemistry tolerance", fastchem_xtol_chem; min=1e-10, max=1.0) || return false
_check_range("FastChem element tolerance", fastchem_xtol_elem; min=1e-10, max=1.0) || return false
_check_range("FastChem chemistry iterations", fastchem_maxiter_chem; min=100, max=1e7) || return false
_check_range("FastChem solver iterations", fastchem_maxiter_solv; min=100, max=1e7) || return false
```

---

#### **ISSUE 4: Condensates Not Validated Against Disallowed List** ⚠️ **SEVERITY: MEDIUM**

**Location**: `src/atmosphere.jl` lines 1055-1061

**Current Code**:
```julia
const COND_DISALLOWED::Array = ["H2","He"]

# Check that condensates are available (lines 1055-1061)
for c in atmos.condensates
    if !(c in atmos.gas_names)
        @error "Condensate '$c' is not in the list of available gases"
        return false
    end
end
```

**Missing Check**: Code never validates against `COND_DISALLOWED` array!

**Risk**: User could set `condensates = ["H2"]` and code would attempt to condense hydrogen, which is physically unrealistic and could cause numerical issues.

**Recommendation**:
```julia
# Add after line 1061:
for c in atmos.condensates
    if c in COND_DISALLOWED
        @error "Gas '$c' is not allowed to condense (disallowed list: $(join(COND_DISALLOWED, ", ")))"
        return false
    end
end
```

---

### 1.3 Configuration Mutual Exclusivity

#### **ISSUE 5: Incomplete Mutual Exclusivity Checking** ⚠️ **SEVERITY: MEDIUM**

**Location**: `src/AGNI.jl` lines 335-340

**Current Check**:
```julia
if comp_set_by > 1
    @error "Config: provide only one of `vmr_dict`, `vmr_file`, `metallicities`"
    return false
end
```

**Problems**:
1. **Vague error**: Doesn't tell user WHICH multiple options they provided
2. **No information about priority**: If user accidentally provides multiple, which takes precedence?

**Improved Implementation**:
```julia
# Track which were provided
provided = []
haskey(cfg["composition"],"vmr_dict") && push!(provided, "vmr_dict")
haskey(cfg["composition"],"vmr_file") && push!(provided, "vmr_file")
haskey(cfg["composition"],"metallicities") && push!(provided, "metallicities")

if length(provided) > 1
    @error "Config: Multiple composition sources provided: $(join(provided, ", "))"
    @error "        Provide ONLY ONE of: vmr_dict, vmr_file, or metallicities"
    @error "        Priority order: metallicities > vmr_file > vmr_dict"
    return false
end
```

---

#### **ISSUE 6: Undocumented Ocean_ini Limitation** ⚠️ **SEVERITY: LOW**

**Location**: `src/atmosphere.jl` line 1015

**Code**:
```julia
atmos.ocean_ini[g] = 0.0  # TODO: make this input from user
```

**Problem**: 
- Users cannot set initial ocean inventories via configuration
- No error message warns users this feature is unavailable
- TODO comment indicates this is a known limitation

**Impact**: If users expect to model planets with pre-existing oceans, they get no guidance that this isn't supported.

**Recommendation**:
1. **Short-term**: Add to documentation and validation:
   ```julia
   # Initial ocean inventory is set to zero (controlled by solver, not user input)
   # To model pre-existing oceans, run solver first, then use output as initial condition
   atmos.ocean_ini[g] = 0.0
   ```
2. **Long-term**: Implement user-configurable initial ocean inventories via TOML parameter `ocean_initial`

---

### 1.4 Documentation Gaps

#### **ISSUE 7: gas_safe Field Purpose Unclear** ⚠️ **SEVERITY: LOW**

**Location**: `src/atmosphere.jl` line 175

**Field Definition**:
```julia
gas_safe::Dict{String, Bool}  # Whether each gas is reliably 'safe', in terms of the physics modelled
```

**Set Only Once** (`atmosphere.jl` lines 1926-1945):
```julia
# Check if this gas is "safe" to use
gas_safe = true
gas_safe = gas_safe && !atmos.gas_dat[gas].stub
gas_safe = gas_safe && !atmos.gas_dat[gas].no_sat
gas_safe = gas_safe && (gas_soc_idx > 0)
atmos.gas_safe[gas] = gas_safe
```

**Criteria for "Safe"**:
1. Not a stub (has thermodynamic data)
2. Has saturation data
3. Appears in SOCRATES spectral file (has opacity data)

**Never Used Afterwards**: Field is set during allocation but never checked or updated.

**Problem**: Purpose and utility unclear from code. Should this trigger warnings? Should unsafe gases be excluded?

**Recommendation**:
```julia
# Add documentation in struct definition:
gas_safe::Dict{String, Bool}  # Flag: gas is "safe" (has opacity, thermo, & saturation data)
                               # Set during allocation; safe=false gases may have limited accuracy
```

---

## 2. Configuration Parameter Audit

### 2.1 Parameter Name Inconsistencies

#### **ISSUE 8: solution_type vs sol_type** ⚠️ **SEVERITY: HIGH**

**Locations**:
- TOML: `res/config/default.toml` line 47 uses `solution_type`
- TOML Comments: lines 17-20 reference `sol_type`
- Code: `src/AGNI.jl` line 418 reads `solution_type` and stores as local variable `sol_type`

**Inconsistency**:
```toml
# Lines 17-20 in default.toml:
    skin_d          = 0.01              # Used when sol_type=2.
    skin_k          = 2.0               # Used when sol_type=2.
    tmp_magma       = 3000.0            # Used when sol_type=2.
    flux_int        = 0.0               # Used when sol_type=3.

# Line 47 in default.toml:
    solution_type   = 0                 # Solution type (see wiki).
```

**Problem**: Users see `sol_type` in parameter descriptions but must use `solution_type` in TOML.

**Recommendation**: 
**Option A** (Preferred): Rename TOML parameter to match comments
```toml
sol_type = 0  # Solution type (see wiki for options)
```

**Option B**: Update all comments to use `solution_type`
```toml
    skin_d = 0.01   # Used when solution_type=2.
```

---

### 2.2 Unit Notation Inconsistencies

#### **ISSUE 9: Temperature Unit Format** ⚠️ **SEVERITY: HIGH**

**Location**: `res/config/default.toml` lines 8, 19, 49

**Inconsistent Notation**:
```toml
Line 8:  tmp_surf = 2000.0            # Surface temperature [kelvin].
Line 19: tmp_magma = 3000.0           # Magma temperature [K].
Line 49: dx_max = 200.0               # Maximum step size [Kelvin], ...
```

Three different formats: `[kelvin]`, `[K]`, `[Kelvin]`

**Recommendation**: Standardize to SI abbreviation `[K]` throughout:
```toml
tmp_surf = 2000.0    # Surface temperature [K]
tmp_magma = 3000.0   # Magma temperature [K]
dx_max = 200.0       # Maximum step size [K]
```

---

#### **ISSUE 10: Pressure Unit Conversion Documentation** ⚠️ **SEVERITY: MEDIUM**

**Location**: `src/energy.jl` lines 24, 681

**Confusing Comment**:
```julia
# Line 24:
const CONVECT_MIN_PRESSURE::Float64 = 1e-9  # lowest pressure at which convection is allowed [bar]

# Line 681:
if atmos.pl[i] <= CONVECT_MIN_PRESSURE * 1.0e5  # convert bar to Pa
```

**Problem**: 
- Constant is in bar but internally converted to Pa
- The `1.0e5` conversion factor is correct (bar → Pa) but buried in code
- Users reading the constant definition might not expect the conversion

**Recommendation**:
```julia
# Option A: Store in Pa directly
const CONVECT_MIN_PRESSURE_PA::Float64 = 1e-4  # lowest pressure [Pa] (= 1e-9 bar)

# Option B: Clarify conversion in constant comment
const CONVECT_MIN_PRESSURE::Float64 = 1e-9  # lowest pressure [bar] (converted to Pa in code)
```

---

#### **ISSUE 11: Log-Pascal Units Unclear** ⚠️ **SEVERITY: LOW**

**Location**: `res/config/default.toml` line 77

**Current**:
```toml
Pwid = 1.0  # Width of Gaussian in log-Pascal units.
```

**Problem**: "log-Pascal units" is non-standard scientific notation.

**Recommendation**:
```toml
Pwid = 1.0  # Width of Gaussian heating profile [log₁₀(Pa)]
```

---

### 2.3 Missing FastChem TOML Parameters

#### **ISSUE 12: FastChem Settings Not Exposed to Users** ⚠️ **SEVERITY: MEDIUM**

**Location**: `src/AGNI.jl` (no TOML keys), `src/atmosphere.jl` lines 512-517

**Current State**: FastChem parameters are hard-coded in `atmosphere.setup!()` with these defaults:
- `fastchem_floor = 400.0`
- `fastchem_maxiter_chem = 80000`
- `fastchem_maxiter_solv = 40000`
- `fastchem_xtol_chem = 1.0e-3`
- `fastchem_xtol_elem = 1.0e-3`
- `fastchem_wellmixed = false`

**Problem**: Users cannot tune FastChem behavior via configuration files.

**Use Cases**:
- High-precision runs may need tighter tolerances
- Exploratory runs may want faster (lower iteration) settings
- Some atmospheres may need higher `fastchem_floor` (e.g., hot Jupiters)

**Recommendation**: Add to TOML schema:
```toml
[chemistry]
    enabled = true                  # Enable FastChem equilibrium chemistry
    floor_temp = 400.0             # Minimum temperature for chemistry [K]
    max_iter_chemistry = 80000     # Maximum chemistry iterations
    max_iter_solver = 40000        # Maximum internal solver iterations
    tolerance_chemistry = 1.0e-3   # Chemistry convergence tolerance
    tolerance_elements = 1.0e-3    # Element conservation tolerance
    wellmixed = false              # Well-mixed (single-level) vs 1D profile
```

---

## 3. User Interface Messages

### 3.1 Critical Error Message Issues

#### **ISSUE 13: Vague Range Validation Errors** ⚠️ **SEVERITY: HIGH**

**Location**: `src/atmosphere.jl` lines 376-390

**Current Implementation**:
```julia
function _check_range(name, val; min=nothing, max=nothing)::Bool
    if !isnothing(min) && !isnothing(max) && ((val<min) || (val>max))
        @error "$name is out of range"
        @error "    Got: $min < $val < $max"
    # ...
end
```

**Problems**:
1. **Ambiguous formatting**: `"Got: $min < $val < $max"` looks like valid range notation but is actually showing the constraint AND the value mixed together
2. **No units**: Error doesn't show what units are expected
3. **Unclear which bound violated**: When `val > max`, still shows `min < val < max`

**Example Error**:
```
ERROR: Planet surface radius is out of range
ERROR:     Got: 10000.0 < 5000.0 < Inf
```
**Interpretation**: Is `10000.0 < 5000.0 < Inf` the valid range, or does it mean `min=10000, val=5000, max=Inf`? **Extremely confusing!**

**Recommendation**:
```julia
function _check_range(name, val; min=nothing, max=nothing, units="")::Bool
    unit_str = units != "" ? " [$units]" : ""
    if !isnothing(min) && !isnothing(max) && ((val<min) || (val>max))
        @error "$name$unit_str is out of valid range"
        @error "    Valid range: $min ≤ value ≤ $max"
        @error "    Actual value: $val"
        return false
    elseif !isnothing(min) && (val < min)
        @error "$name$unit_str is too small"
        @error "    Minimum allowed: $min"
        @error "    Actual value: $val"
        return false
    elseif !isnothing(max) && (val > max)
        @error "$name$unit_str is too large"
        @error "    Maximum allowed: $max"
        @error "    Actual value: $val"
        return false
    end
    return true
end
```

---

#### **ISSUE 14: Undefined Abbreviations in Error Messages** ⚠️ **SEVERITY: MEDIUM**

**Locations**: `src/atmosphere.jl` lines 600-601, 711

**Examples**:

1. **Line 711**:
   ```julia
   _check_range("Surface CBL thickness", atmos.skin_d; min=SKIN_D_MIN)
   ```
   **Problem**: "CBL" never defined. Users don't know this means "Convective Boundary Layer".

2. **Line 600-601**:
   ```julia
   if atmos.flag_rayleigh
       @error "Scattering not supported by grey-gas RT scheme!"
   ```
   **Problem**: "RT" abbreviation undefined (Radiative Transfer).

**Recommendation**:
```julia
# Line 711: Spell out on first use
_check_range("Surface conductive boundary layer (CBL) thickness", atmos.skin_d; min=SKIN_D_MIN)

# Line 600-601: Expand abbreviation
@error "Scattering not supported by grey-gas radiative transfer (RT) scheme!"
```

---

#### **ISSUE 15: Incomplete Error Context** ⚠️ **SEVERITY: MEDIUM**

**Location**: `src/AGNI.jl` lines 454-455, 466-467, 476-477

**Current Messages**:
```julia
# Line 454-455:
@error "Config: solution type $sol_type selected"
@error "        you must provide `planet.skin_k`, `skin_d`, `tmp_magma`"

# Line 466-467:
@error "Config: solution type $sol_type selected"
@error "        you must provide `planet.flux_int`"

# Line 476-477:
@error "Config: solution type $sol_type selected"
@error "        you must provide `planet.target_olr`"
```

**Problems**:
1. Doesn't explain WHAT each solution type is
2. Doesn't suggest where to find documentation
3. Users don't know why they need these parameters

**Recommendation**:
```julia
# Line 454-455:
@error "Config: solution type $sol_type (conductive skin with subsurface coupling)"
@error "        Required parameters: `planet.skin_k`, `planet.skin_d`, `planet.tmp_magma`"
@error "        See documentation: Solution Types > Type 2"

# Line 466-467:
@error "Config: solution type $sol_type (fixed internal flux)"
@error "        Required parameter: `planet.flux_int` [W m⁻²]"
@error "        See documentation: Solution Types > Type 3"

# Line 476-477:
@error "Config: solution type $sol_type (target outgoing longwave radiation)"
@error "        Required parameter: `planet.target_olr` [W m⁻²]"
@error "        See documentation: Solution Types > Type 4"
```

---

### 3.2 Formatting Inconsistencies

#### **ISSUE 16: Inconsistent Multi-Line Error Formatting** ⚠️ **SEVERITY: MEDIUM**

**Locations**: Throughout `src/AGNI.jl`

**Inconsistent Patterns**:

1. **Lowercase continuation**:
   ```julia
   # Line 336:
   @error "Config: if providing p_surf, must also provide VMRs"
   ```

2. **Indented with lowercase**:
   ```julia
   # Lines 443-444:
   @error "Config: sensible heating included"
   @error "        you must provide `planet.roughness` and `planet.wind_speed`"
   ```

3. **Indented with uppercase**:
   ```julia
   # Lines 383-384:
   @error "Config: RFM calculation enabled (rfm_parfile=$rfm_parfile)"
   @error "        You must also provide `rfm_wn_min` AND `rfm_wn_max`"
   ```

**Recommendation**: Standardize to consistent format:
```julia
@error "Config: [Context description]"
@error "        Required: [what's needed]"
@error "        See: [documentation reference]"
```

Example:
```julia
@error "Config: Sensible heating included"
@error "        Required: `planet.roughness` and `planet.wind_speed`"
@error "        See: Physics Options > Sensible Heat"
```

---

#### **ISSUE 17: Logging Inconsistency (@error vs error())** ⚠️ **SEVERITY: MEDIUM**

**Locations**: Multiple files

**Two Different Approaches**:

1. **Using `@error` macro + `return false`**:
   ```julia
   # Line 231 in AGNI.jl:
   @error "Config: missing required key `$k.$kk`"
   return false
   ```
   **Effect**: Logs error, allows graceful exit

2. **Using `error()` function**:
   ```julia
   # Line 187 in AGNI.jl:
   error("Key $h is missing from configuration file at '$cfg_path'")
   ```
   **Effect**: Throws exception immediately

**Problem**: Inconsistent behavior — some config errors allow cleanup, others abort immediately.

**Recommendation**: Standardize to `@error` + `return false` pattern for all configuration validation:
```julia
# Convert all error() → @error + return false
if !haskey(cfg_dict, h)
    @error "Config: Missing required section '$h' in file: $cfg_path"
    return false
end
```

---

### 3.3 Missing Context and Guidance

#### **ISSUE 18: File Not Found Errors Lack Actionable Guidance** ⚠️ **SEVERITY: LOW**

**Location**: `src/spectrum.jl` lines 49, 201, 337, 341

**Current Messages**:
```julia
@error "Spectral file not found: '$spec_file'"
@error "Stellar spectrum file not found: '$star_file'"
@error "Original spectral file not found: '$orig_file'"
```

**Problem**: Doesn't tell users what to do or check.

**Recommendation**:
```julia
@error "Spectral file not found: '$spec_file'"
@error "    Check that the file exists and the path is correct"
@error "    Hint: Relative paths are relative to working directory: $(pwd())"
@error "    Expected absolute path: $(abspath(spec_file))"
```

---

#### **ISSUE 19: Transspec Term Never Defined** ⚠️ **SEVERITY: LOW**

**Location**: `src/atmosphere.jl` lines 738, 740-742

**Code**:
```julia
# Line 738:
atmos.transspec_p = 2e3  # 20 mbar = 2000 Pa

# Line 740-742:
if atmos.p_toa > atmos.transspec_p
    @error "p_top must be less than transspec_p"
    return false
end
```

**Problem**: "transspec" abbreviation never explained in error message or nearby comments.

**Users see**: `"p_top must be less than transspec_p"` with no idea what `transspec_p` is.

**Recommendation**:
```julia
# Line 738-742:
atmos.transspec_p = 2e3  # Transmission spectroscopy reference pressure: 20 mbar = 2000 Pa

if atmos.p_toa > atmos.transspec_p
    @error "Top-of-atmosphere pressure (p_top) must be less than transmission spectroscopy pressure"
    @error "    p_top = $(atmos.p_toa/1e5) bar"
    @error "    transspec_p = $(atmos.transspec_p/1e5) bar (20 mbar reference level)"
    return false
end
```

---

#### **ISSUE 20: MLT Criterion Error Doesn't Explain Options** ⚠️ **SEVERITY: LOW**

**Location**: `src/atmosphere.jl` lines 703-705

**Current**:
```julia
if !(atmos.mlt_criterion in ['s','l'])
    @error "Invalid choice for mlt_criterion: $(atmos.mlt_criterion)"
    @error "    Must be: 's' or 'l' only"
    return false
end
```

**Problem**: Doesn't tell user what 's' and 'l' mean.

**Recommendation**:
```julia
if !(atmos.mlt_criterion in ['s','l'])
    @error "Invalid choice for mlt_criterion: '$(atmos.mlt_criterion)'"
    @error "    Valid options:"
    @error "        's' = Schwarzschild criterion (neglects composition gradients)"
    @error "        'l' = Ledoux criterion (accounts for mean molecular weight gradients)"
    return false
end
```

---

### 3.4 Dead Code and TODOs

#### **ISSUE 21: Commented-Out Code Should Be Removed** ⚠️ **SEVERITY: LOW**

**Location**: `src/atmosphere.jl` lines 700-701

**Code**:
```julia
@warn "    (Will use Ledoux criterion anyway)"
# @warn "    Switching criterion to Schwarzschild, neglecting MMW gradients"
# atmos.mlt_criterion = 's'
```

**Problem**: Dead code creates confusion. Is this planned for future? Should users expect this behavior?

**Recommendation**: Either:
1. **Remove entirely** if not planned
2. **Add issue tracker reference** if planned:
   ```julia
   # TODO(Issue #123): Implement Schwarzschild fallback for real gas
   # @warn "    Switching criterion to Schwarzschild, neglecting MMW gradients"
   # atmos.mlt_criterion = 's'
   ```

---

## 4. Critical Issues Summary

### By Severity

| Severity | Count | Issues |
|----------|-------|--------|
| **HIGH** | 4 | #1 (VMR naming), #2 (metallicity units), #8 (solution_type), #9 (temperature units) |
| **MEDIUM** | 9 | #3 (FastChem validation), #4 (condensates), #5 (mutual exclusivity), #10 (pressure units), #12 (FastChem TOML), #14 (abbreviations), #15 (error context), #16 (formatting), #17 (logging) |
| **LOW** | 8 | #6 (ocean_ini TODO), #7 (gas_safe), #11 (log-Pascal), #13 (range errors), #18 (file errors), #19 (transspec), #20 (MLT error), #21 (dead code) |

### By Category

| Category | Count | Issues |
|----------|-------|--------|
| **Chemistry Scheme** | 7 | #1, #2, #3, #4, #5, #6, #7 |
| **Configuration** | 5 | #8, #9, #10, #11, #12 |
| **Error Messages** | 9 | #13, #14, #15, #16, #17, #18, #19, #20, #21 |

---

## 5. Recommendations

### Immediate Actions (HIGH Priority)

1. **Fix gas VMR field naming** (#1)
   - Add extensive inline documentation to `Atmos_t` struct
   - Consider renaming to `gas_vmr_runtime`, `gas_vmr_chemistry`, `gas_vmr_initial`

2. **Clarify metallicity units** (#2)
   - Fix TOML comment: "molar ratios relative to H" (not "mass fractions")
   - Fix docstring in `atmosphere.jl`
   - Add validation for reasonable ranges

3. **Standardize solution_type vs sol_type** (#8)
   - Rename TOML parameter to `sol_type` OR
   - Update all comments to reference `solution_type`

4. **Standardize temperature units** (#9)
   - Use `[K]` consistently across all TOML files
   - Update documentation to match

### Short-Term Actions (MEDIUM Priority)

5. **Add FastChem parameter validation** (#3)
   - Validate `fastchem_floor`, iteration counts, tolerances
   - Add bounds checks with clear error messages

6. **Validate condensates against disallowed list** (#4)
   - Check against `COND_DISALLOWED` array
   - Provide clear error when H2 or He condensation attempted

7. **Improve mutual exclusivity checking** (#5)
   - Show which multiple options were provided
   - Document priority order

8. **Expose FastChem settings to TOML** (#12)
   - Add `[chemistry]` section with all FastChem parameters
   - Allow per-configuration tuning

9. **Improve error message clarity** (#13-15)
   - Fix range validation formatting
   - Define all abbreviations
   - Add context for solution type requirements

10. **Standardize error message formatting** (#16-17)
    - Use consistent multi-line format
    - Switch all config errors to `@error` + `return false`

### Long-Term Actions (LOW Priority)

11. **Implement ocean_ini user input** (#6)
    - Add TOML parameter for initial ocean inventories
    - Update validation and documentation

12. **Document gas_safe criteria** (#7)
    - Clarify purpose in comments
    - Consider adding warnings for unsafe gases

13. **Improve all file-not-found errors** (#18)
    - Add actionable guidance
    - Show absolute paths
    - Suggest common fixes

14. **Remove dead code** (#21)
    - Delete commented code or link to issue tracker
    - Clean up TODOs with action items

---

## Appendix A: Complete Field Reference

### Gas Mixing Ratio Fields

| Field | Type | Purpose | Updated By |
|-------|------|---------|------------|
| `gas_ovmr` | Dict{String, Array{Float64,1}} | Original user input | `setup!()`, never modified |
| `gas_cvmr` | Dict{String, Array{Float64,1}} | Post-chemistry VMRs | `_chem_gas!()`, `restore_composition!()` |
| `gas_vmr` | Dict{String, Array{Float64,1}} | Runtime VMRs | `_sat_surf!()`, `_sat_aloft!()`, chemistry, etc. |

### Metallicity Fields

| Field | Type | Purpose | Units |
|-------|------|---------|-------|
| `metal_orig` | Dict{String, Float64} | User-provided elem ratios | Molar number rel. to H |
| `metal_calc` | Dict{String, Float64} | Calculated elem ratios | Molar number rel. to H |

### FastChem Configuration Fields

| Field | Type | Default | Units | Validated? |
|-------|------|---------|-------|------------|
| `fastchem_floor` | Float64 | 400.0 | K | ❌ No |
| `fastchem_maxiter_chem` | Int | 80000 | iterations | ❌ No |
| `fastchem_maxiter_solv` | Int | 40000 | iterations | ❌ No |
| `fastchem_xtol_chem` | Float64 | 1.0e-3 | dimensionless | ❌ No |
| `fastchem_xtol_elem` | Float64 | 1.0e-3 | dimensionless | ❌ No |
| `fastchem_wellmixed` | Bool | false | - | ✓ Bool |

---

## Appendix B: Error Message Pattern Guide

### Recommended Pattern

```julia
# Configuration Errors:
@error "Config: [Brief context]"
@error "        Required: [specific requirements]"
@error "        See: [documentation reference]"
return false

# Validation Errors:
@error "[Parameter name] [issue]"
@error "    Valid range: $min ≤ value ≤ $max [$units]"
@error "    Actual value: $val [$units]"
return false

# File Errors:
@error "[File type] not found: $path"
@error "    Check that file exists and path is correct"
@error "    Working directory: $(pwd())"
@error "    Absolute path: $(abspath(path))"
return false
```

### Examples

```julia
# Good: Clear, actionable, consistent
@error "Config: Sensible heating requires boundary layer parameters"
@error "        Required: `planet.roughness` [m] and `planet.wind_speed` [m s⁻¹]"
@error "        See: Physics Options > Sensible Heat Transport"
return false

# Bad: Vague, no context, no units
@error "Config: sensible heating included"
@error "        you must provide roughness and wind speed"
return false
```

---

**End of Second-Pass Audit Report**
