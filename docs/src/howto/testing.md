# Testing in AGNI

This document describes the testing infrastructure and how to run tests for AGNI.

## Running Tests

### Prerequisites
Before running tests, ensure that the `RAD_DIR` environment variable is set to your SOCRATES installation:

```bash
export RAD_DIR=/path/to/SOCRATES
```

### Run All Tests
From the `test/` directory:

```bash
cd test/
julia --project=.. runtests.jl
```

### Run Fast Tests Only
The fast test suite excludes expensive integration tests:

```bash
cd test/
julia --project=.. runtests.jl fast
```

### Run a Single Test File
To run an individual test file:

```bash
cd test/
julia --project=.. -e 'include("test_phys.jl")'
```

## Test Coverage

### Generating Coverage Reports

#### Run Tests with Coverage
```bash
cd test/
julia --project=.. --code-coverage runtests.jl fast
```

#### Process Coverage Data
From the AGNI root directory:

```bash
julia --project=. test/get_coverage.jl
```

This generates:
- `coverage.info` - LCOV format coverage data
- `coverage.json` - Codecov JSON format
- `coverage.total` - Single number with overall coverage percentage

#### Generate Markdown Report
To create a human-readable coverage report:

```bash
julia --project=. test/coverage_to_markdown.jl
```

This produces `coverage.md` with:
- Overall coverage percentage
- Per-file coverage breakdown with color coding (🔴 <50%, 🟡 50-80%, 🟢 >80%)
- Files needing attention (<50% coverage)
- Quick wins (small files with 0% coverage)
- Lists of uncovered line numbers

## Test Suite Structure

The test suite is organized into separate test files, each focusing on a specific module:

| Test File | Module Tested | Description |
|-----------|---------------|-------------|
| `test_consts.jl` | `consts` | Physical constants, molecular weights, atom counts, standard species lists |
| `test_phys.jl` | `phys` | Physical utility functions: EOS, thermodynamics, gas properties, molecular formulas |
| `test_blake.jl` | `phys.blake` | BLAKE2b file hashing and integrity checking |
| `test_spectrum.jl` | `spectrum` | Spectral file operations, gas counting, SOCRATES version |
| `test_setpt.jl` | `setpt` | Temperature profile initialization functions |
| `test_guillot.jl` | `setpt.guillot` | Analytical Guillot (2010) T(p) profiles |
| `test_ocean.jl` | `ocean` | Surface liquid distribution and ocean geometry |
| `test_deep_heating.jl` | (various) | Deep atmospheric heating profiles |
| `test_kzz.jl` | (various) | Eddy diffusion coefficient (Kzz) profiles |
| `test_integration.jl` | (full model) | End-to-end integration tests (slow) |

### Test Design Principles

1. **Fast by default**: Most tests avoid expensive operations like full atmosphere allocation
2. **Unit-focused**: Each test targets specific functions or small groups of related functions
3. **Self-contained**: Tests don't depend on external files when possible (except SOCRATES spectral files)
4. **Comprehensive**: Tests cover:
   - Normal operation (happy path)
   - Edge cases (boundary conditions, empty inputs)
   - Error conditions (missing files, invalid inputs)
   - Type stability and return value validation

### Example Test Pattern

```julia
using Test
using AGNI

@testset "module_name" begin
    @testset "function_name" begin
        # Test normal case
        result = AGNI.module.function(valid_input)
        @test result > 0
        @test typeof(result) == Float64

        # Test edge case
        result_edge = AGNI.module.function(edge_case_input)
        @test isfinite(result_edge)

        # Test error handling
        @test_throws ErrorException AGNI.module.function(invalid_input)
    end
end
```

## Adding New Tests

When adding new tests:

1. Create a new test file in `test/` following the naming convention `test_<module>.jl`
2. Add it to `test/runtests.jl` in the appropriate location
3. Follow the existing test patterns (see above)
4. Avoid expensive operations:
   - Don't call `allocate!`, `setup!`, or `deallocate!` unless absolutely necessary
   - Use minimal atmosphere configurations when needed
   - Keep test execution time under 5 seconds per file
5. Clean up temporary files created during tests
6. Test both success and failure paths
7. Use descriptive test names in `@testset` blocks

### Coverage Goals

- **Target**: >80% line coverage for utility modules (`phys`, `consts`, `ocean`, `guillot`, `blake`)
- **Acceptable**: >50% for core simulation modules (`atmosphere`, `energy`, `solver`)
- **Document**: Any intentionally untested code (e.g., platform-specific branches)

## Continuous Integration

GitHub Actions automatically:
1. Runs the fast test suite on every push
2. Generates coverage reports
3. Updates the coverage badge in the README

See `.github/workflows/install_and_test.yml` for CI configuration.

## Troubleshooting

### Tests fail with "RAD_DIR not set"
Ensure `RAD_DIR` points to a valid SOCRATES installation:
```bash
export RAD_DIR=/path/to/SOCRATES
```

### Tests fail with "Spectral file not found"
The test suite expects spectral files in `res/spectral_files/`. Ensure AGNI's resource directory is complete.

### Coverage is 0% or unexpectedly low
- Make sure you ran tests with `--code-coverage` flag
- Run `get_coverage.jl` from the AGNI root directory (not `test/`)
- Check that `.cov` files were generated in `src/` and `test/`

## Test Maintenance

- Run tests before committing changes
- Update tests when modifying function signatures or behavior
- Keep test execution time reasonable (fast suite <30 seconds)
- Regularly review coverage reports to identify untested code paths
