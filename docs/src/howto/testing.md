# Testing in AGNI

This document describes the testing infrastructure and how to run tests for AGNI.

## Running Tests

### Prerequisites
Before running tests, ensure that the `RAD_DIR` environment variable is set to your SOCRATES installation:

```bash
export RAD_DIR=/path/to/SOCRATES
```

### Run All Tests
From the AGNI root directory:

```bash
julia --project=. test/runtests.jl
```

### Run Fast Tests Only
The fast test suite excludes expensive integration tests:

```bash
julia --project=. test/runtests.jl fast
```

## Test Coverage

### Generating Coverage Reports

#### Run Tests with Coverage
From the AGNI root directory:

```bash
julia --project=. --code-coverage test/runtests.jl
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
- `coverage.md` - Markdown formatted file with report on coverage

The `coverage.md` file outlines:
- Overall coverage percentage
- Per-file coverage breakdown with color coding (🔴 <50%, 🟡 50-80%, 🟢 >80%)
- Files needing attention (<50% coverage)
- Quick wins (small files with 0% coverage)
- Lists of uncovered line numbers

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
2. Any `test_*.jl` file in `test/` will be picked up automatically by `test/runtests.jl`, so no manual edit to `test/runtests.jl` is needed
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

## Troubleshooting the tests

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
