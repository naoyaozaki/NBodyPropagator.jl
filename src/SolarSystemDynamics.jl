""" Solar System Dynamics Library in Julia
"""

# Export Function
export get_path_of_genker, init_spice_kernels
export SolarSystemDynamics

# Import Modules
import SPICE
using Downloads

"""
Solar System Dynamics Parameter Sets (Struct)

Attributes:
    AU (Float64)               : Astronomical unit, km
    ID (Dict{String,Integer})  : ID of planetary bodies where the key is the body name.
    GM (Dict{String,Float64})  : Gravity constant of planetary bodies where the key is the body name, km3/s2.
    RE (Dict{String,Float64})  : Equitorial radius of planetary bodies where the key is the body name, km.
    NAME (Dict{String,String}) : NAME of planetary bodies where the key is the body id.
    OE (Dict{String,Dict})     : Orbital elements of planets.
"""
struct SolarSystemDynamics
    AU::Float64
    ID::Dict{String,Integer}
    GM::Dict{String,Float64}
    RE::Dict{String,Float64}
    NAME::Dict{Integer,String}
    OE::Dict{String,Dict}

    """ Constructor of the SolarSystemDynamics struct
    """
    function SolarSystemDynamics()

        # 1) Furnish SPICE Kernels
        init_spice_kernels()

        # 2) List of Planetary Bodies
        list_bodies = ["SOLAR_SYSTEM_BARYCENTER",
            "MERCURY BARYCENTER",
            "VENUS BARYCENTER",
            "EARTH BARYCENTER",
            "MARS BARYCENTER",
            "JUPITER BARYCENTER",
            "SATURN BARYCENTER",
            "URANUS BARYCENTER",
            "NEPTUNE BARYCENTER",
            "PLUTO BARYCENTER",
            "SUN",
            "MERCURY",
            "VENUS",
            "EARTH",
            "MOON",
            "MARS",
            "PHOBOS",
            "DEIMOS",
            "JUPITER",
            "IO",
            "EUROPA",
            "GANYMEDE",
            "CALLISTO",
            "SATURN",
            "MIMAS",
            "ENCELADUS",
            "TETHYS",
            "DIONE",
            "RHEA",
            "TITAN",
            "URANUS",
            "NEPTUNE",
            "PLUTO"]

        _AU = 149597870.7  # km

        # 3) Load SPICE Kernels
        _ID, _GM, _RE, _NAME = Dict(), Dict(), Dict(), Dict()
        for body in list_bodies
            id_temp = SPICE.bodn2c(body)
            if !isnothing(id_temp)
                _ID[body] = id_temp
                if SPICE.bodfnd(_ID[body], "GM")
                    _GM[body] = SPICE.bodvcd(_ID[body], "GM", 1)[1]
                end
                if SPICE.bodfnd(_ID[body], "RADII")
                    _RE[body] = SPICE.bodvcd(_ID[body], "RADII", 3)[1]
                end
                _NAME = Dict(value => key for (key, value) in _ID)
            end
        end


        # 4) Mean Orbital Elements
        _OE = Dict{String,Dict}()
        list_oe = ["MERCURY", "VENUS", "EARTH", "MARS", "JUPITER", "SATURN",
            "URANUS", "NEPTUNE", "PLUTO"]

        current_directory = dirname(abspath(@__FILE__))
        open(current_directory * "/../data/lib/p_elem_t1.txt", "r") do f
            list = readlines(f)
            for id in 1:9
                _elem = parse.(Float64, split(list[15+2*id][8:end]))
                _OE[list_oe[id]] = Dict("sma" => _elem[1] * _AU, #= km =#
                    "ecc" => _elem[2], #= non-dim =#
                    "inc" => _elem[3] * (π / 180), #= rad =#
                    "meanlong" => _elem[4] * (π / 180), #= rad =#
                    "argp" => _elem[5] * (π / 180), #= rad =#
                    "lnode" => _elem[6] * (π / 180)) #= rad =#
            end
        end

        new(_AU, _ID, _GM, _RE, _NAME, _OE)
    end
end


"""Get the path of specific generic kernels.

Args:
    kernel_name (String): Name of the generic kernels

Returns:
    String: Full path of the generic kernels

Note:
    The downloaded generic kernels are supposed to 
    be stored in <current_directory> * "/../data/lib"
"""
function get_path_of_genker(kernel_name::String)::String
    # Get the spice directory
    current_directory = dirname(abspath(@__FILE__))
    spice_directory = current_directory * "/../data/lib"

    # Split the input string
    kernel_dir = split(kernel_name, "/")

    # Find naif0012.tls in the directory
    is_downloaded = false
    genker_path = ""
    for (root, dirs, files) in walkdir(spice_directory)
        is_downloaded = any(occursin.(kernel_dir[end], files))
        if is_downloaded
            genker_path = root * "/" * kernel_dir[end]
            break
        end
    end
    if !is_downloaded
        mkpath(dirname(spice_directory * "/" * kernel_name))
        genker_path = Downloads.download("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/" * kernel_name, spice_directory * "/" * kernel_name)
    end

    return genker_path
end

"""Initialize SPICE library

Note:
    The downloaded generic kernels are supposed to 
    be stored in <current_directory> * "/../data/lib"
"""
function init_spice_kernels()
    # Initialize SPICE Kernels
    SPICE.furnsh(get_path_of_genker("pck/earth_200101_990825_predict.bpc")) # Earth-moon system kernel
    SPICE.furnsh(get_path_of_genker("pck/earth_fixed.tf")) # Earth fixed system
    SPICE.furnsh(get_path_of_genker("lsk/naif0012.tls")) # Leap seconds kernel
    SPICE.furnsh(get_path_of_genker("pck/gm_de440.tpc")) # Gravity Constant
    SPICE.furnsh(get_path_of_genker("pck/pck00011.tpc")) # P-Constant
    SPICE.furnsh(get_path_of_genker("spk/satellites/nep102.bsp")) # Neptune system kernel
    SPICE.furnsh(get_path_of_genker("spk/satellites/sat452.bsp")) # Saturn system kernel
    SPICE.furnsh(get_path_of_genker("spk/satellites/mar097.bsp")) # Mars system kernel
    SPICE.furnsh(get_path_of_genker("spk/satellites/ura116.bsp")) # Uranus system kernel
    SPICE.furnsh(get_path_of_genker("spk/satellites/jup365.bsp")) # Jupiter system kernel
    SPICE.furnsh(get_path_of_genker("spk/planets/de440.bsp")) # Planetary ephemeris kernel
end
