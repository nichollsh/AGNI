
include("../src/ParseFortran.jl")

SOCRATES_DIR = "../../"

strctrl = ParseFortran.parse_type(joinpath(SOCRATES_DIR, "src/radiance_core/def_control.F90"), "StrCtrl")

int_names = [f[:name] for f in strctrl if f[:type]=="integer"]
println("int_names:\n", int_names)

bool_names = [f[:name] for f in strctrl if f[:type]=="logical"]
println("bool_names:\n", bool_names)

real_names = [f[:name] for f in strctrl if contains(f[:type], "real")]
println("real_names:\n", real_names)

alloc_names = [f[:name] for f in strctrl if f[:allocatable]]
println("alloc_names:\n", alloc_names)

vec_names = [f[:name] for f in strctrl if !isempty(f[:dims])]
println("vec_names:\n", vec_names)

string_names = [f[:name] for f in strctrl if contains(f[:type], "char")]
println("string_names:\n", string_names)

strspectrum = ParseFortran.parse_type(joinpath(SOCRATES_DIR, "src/radiance_core/def_spectrum.F90"), "StrSpecData")

pars = ParseFortran.parse_module_parameters(joinpath(SOCRATES_DIR, "src/radiance_core/rad_pcf.F90"), "rad_pcf")
