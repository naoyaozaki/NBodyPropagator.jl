""" Solar System Dynamics Library in Julia
"""

# Export Function
export get_path_of_genker
export ssd

# Import Modules
import SPICE
using Downloads

"""Solar System Dynamics Parameter Sets (Struct)

Attributes:
    AU                : Astronomical unit, km
    ID(dict) : ID of planetary bodies where the key is the body name.
    NAME(dict) : NAME of planetary bodies where the key is the body id.
    GM(dict) : Gravity constant of planetary bodies where the key is the body name, km3/s2.
    RE(dict) : Equitorial radius of planetary bodies where the key is the body name, km.
"""
struct SolarSystemDynamics
    AU::Real
    ID::Dict{String,Integer}
    GM::Dict{String,Real}
    RE::Dict{String,Real}
    NAME::Dict{Integer,String}

    """ Constructor of the SolarSystemDynamics struct
    """
    function SolarSystemDynamics()

        # 1) Furnish SPICE Kernels
        # SPICE.furnsh(get_path_of_genker("spk/satellites/jup310.bsp")) # Jupiter system kernel
        # SPICE.furnsh(get_path_of_genker("spk/satellites/ura090_1.bsp")) # Uranus system kernel
        # SPICE.furnsh(get_path_of_genker("spk/satellites/mar097.bsp")) # Mars system kernel
        # SPICE.furnsh(get_path_of_genker("spk/satellites/sat317.bsp")) # Saturn system kernel
        # SPICE.furnsh(get_path_of_genker("spk/satellites/nep077.bsp")) # Neptune system kernel
        # SPICE.furnsh(get_path_of_genker("spk/satellites/nep081.bsp")) # Neptune system kernel
        SPICE.furnsh(get_path_of_genker("pck/earth_fixed.tf")) # Earth fixed system
        SPICE.furnsh(get_path_of_genker("pck/earth_200101_990628_predict.bpc")) # Earth-moon system kernel
        SPICE.furnsh(get_path_of_genker("lsk/naif0012.tls")) # Leap seconds kernel
        SPICE.furnsh(get_path_of_genker("pck/gm_de431.tpc")) # Gravity Constant
        SPICE.furnsh(get_path_of_genker("pck/pck00010.tpc")) # P-Constant
        SPICE.furnsh(get_path_of_genker("spk/planets/de430.bsp")) # Planetary ephemeris kernel

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

        new(_AU, _ID, _GM, _RE, _NAME)
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


ssd = SolarSystemDynamics()