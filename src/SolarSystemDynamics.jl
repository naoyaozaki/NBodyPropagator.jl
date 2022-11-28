""" Solar System Dynamics Library in Julia
"""

# Export Function
export get_path_of_genker

# Import Modules
using Downloads


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